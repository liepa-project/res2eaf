#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import string
import argparse

import webvtt
from pympi import Eaf
from pathlib import Path
import statistics as stats
from itertools import pairwise
from prettytable import PrettyTable
from datetime import timedelta, datetime, timezone

__author__ = "Laimonas Vėbra"
__copyright__ = "Copyright 2024"
__license__ = "BSD"
__version__ = "0.1a"


parser = argparse.ArgumentParser(
    prog='res2eaf',
    description='Converts Semantika transcription results to EAF',)

parser.add_argument('-v', '--verbose', action='store_true',
                    help='print verbose information')
parser.add_argument('-d', '--debug', action='store_true',
                    help='print debug information')
parser.add_argument('--version', action='version',
                    version='%(prog)s {version}'.format(version=__version__))


group = parser.add_argument_group('Input/Output files')
group.add_argument('-l', '--lattice', metavar='F_LAT', required=True,
                   help='Sematika recognizer result (lattice file)')
group.add_argument('-w', '--webvtt', metavar='F_VTT', required=False,
                   help='webvtt file; if given it will be read & parsed for '
                   'segments along with F_LAT for speaker info, otherwise '
                   'only F_LAT is sufficient')
group.add_argument('-o', '--outfile', metavar='F_EAF', required=True,
                   help='eaf file to save the output')



group = parser.add_argument_group('General')
group.add_argument('--skip-overlaping', action='store_true',
                   help='Skip/omit overlaping segments from all tiers')

group = parser.add_argument_group('Joining segments')
group.add_argument('-j', '--join-segments', action='store_true', 
                   help='join consequitive segments. Note: segments up to '
                   '--max-length ar joined if they are separated by gaps '
                   '<= --max-gap and then, if --allow-overlength, up to '
                   '--ultimate-length if gaps <= --max-overlength-gap')

group.add_argument('--max-gap', metavar='INT_MS', required=False,
                   type=int, default=100,
                   help='max interval (in ms) between segments to be joined '
                   '(default: %(default)d ms)')

group.add_argument('--max-length', metavar='LEN_MS', required=False,
                   type=int, default=5000,
                   help='max length (in ms) of the combined segment if '
                   'overlength is not allowed (default: %(default)d ms)')

group.add_argument('--allow-overlength', action='store_true', 
                   help='Allow combining segments over --max-length up to '
                   '--ultimate-length')

group.add_argument('--ultimate-length', metavar='LEN_MS', required=False,
                   type=int, default=7000,
                   help='Ultimate length (in ms) of the combined segment if '
                   'overlength is allowed (default: %(default)d ms)')

group.add_argument('--max-overlength-gap', metavar='INT_MS', required=False,
                   type=int, default=50,
                   help='max interval (in ms) between segments over '
                   '--max-length to be joined (default: %(default)d ms)')

group = parser.add_argument_group('Text processing')
group.add_argument('--strip-punctuation', action='store_true', 
                   help='remove punctuation characters from segments text')

group = parser.add_argument_group('EAF creation')
group.add_argument('--author',  required=False,
                   default='{0}:{1}'.format(Path(__file__).stem, __version__),
                   help='Set author of the EAF (default: %(default)s)')
                   
group.add_argument('--link-media', metavar='F_WAV', required=False,
                   help='media/wav file to link to in EAF')
group.add_argument('--orig-media', required=False,
                   help='Original media/wav file name (--link-media may link '
                   'to renamed file). This name is stored in EXTRACTED_FROM '
                   'att. and is preserved when EAF is saved in ELAN')
group.add_argument('--annotator', help='Add annotator to all tiers')
group.add_argument('--prefill-meta', metavar='INFO',
                   help='Prefill all tiers with meta info in participant field')
group.add_argument('--overlap-tier', action='store_true', 
                   help='Add overlap tier with overlaping speech intervals')




args = parser.parse_args()


if not Path(args.lattice).is_file():
    print("The specified lattice file '{0}' does not exist"
          .format(args.lattice))
    sys.exit(1)

if (args.webvtt and not Path(args.webvtt).is_file()):
    print("The specified webvtt file '{0}' does not exist"
          .format(args.webvtt))
    sys.exit(1)
    

def to_ms(ts):
     return int(float(ts) * 1000)

def ms_to_ts(ms):
    """ Converts milliseconds to time string in [HH:]MM:SS.mmm format """
    dt = datetime.fromtimestamp(ms/1000, tz=timezone.utc)

    # NOTE: up to 24 hours
    if dt.hour >= 1:
        return dt.strftime('%H:%M:%S.%f')[:-3]
    else:
        return dt.strftime('%M:%S.%f')[:-3]
    

    
speech = {}
speech_blocks = {}
overlaps = []


with open(args.lattice, 'r', encoding='utf-8') as lat_file:
    header = re.compile(
        r"^#\s+(?P<blk>\d+)\s+"
        r"(?P<sid>.+)$")

    segment = re.compile(
        r"^(?P<hyp>1)\s+"
        r"(?P<beg>0|\d+(\.\d{1,2})?)\s+"
        r"(?P<end>\d+(\.\d{1,2})?)\s+"
        r"(?P<val>.+)$")

    sid = None
    
    last_blk = 0
    last_sid = ''
    last_end = 0
    lineno = 0
    overlap = False
        
    for line in lat_file:
        lineno += 1
        
        if line.strip() == '':
            continue

        elif line.startswith('#'):
            if (m := header.match(line)):
                
                # New speach block
                sid  = m.group('sid')
                blk = int(m.group('blk'))

                # Skip fix.lattice.time inserted silence (TYLA) blocks
                if sid == 'TYLA':
                    last_sid, last_blk = sid, blk
                    continue
                
                if sid not in speech:
                    speech[sid] = []

                speech_blocks[blk] = { 'sid': sid, 'segs': []}
                

                if (blk != last_blk + 1):
                    print("WARN: non-consequitive speech block number\n"
                          "    prev. sid: {0}, seq: {1}\n"
                          "    curr. sid: {2}, seq: {3}".format(last_sid, last_blk, sid, blk))

                if (last_sid and sid == last_sid):
                    print("WARN: same speaker '{0}' "
                          "consequitive blocks {1} and {2}".format(last_blk, blk))
                    
                last_sid, last_blk = sid, blk

            else:
                print("WARN: Line:{0} '{1}' doesn't match header format".format(lineno, line))
        
        elif (m := segment.match(line)):
            assert (sid is not None) # should never happen: speech line without (prior) header

            # Skip fix.lattice.time inserted silence (TYLA) blocks/segments
            # they are empty anyway, but https://github.com/airenas/list/issues/1
            if sid == 'TYLA':
                continue

            seg = {
                'hyp': m.group('hyp'),
                'beg': to_ms(m.group('beg')),
                'end': to_ms(m.group('end')),
                'val': m.group('val')
            }

            # Skip silence/noise segments
            if seg['val'].strip() == '<eps>':
                continue
                
            speech[sid].append(seg)
            speech_blocks[blk]['segs'].append(list(seg.values()))

            # Overlaping segments. XXX: based on asumption that speech
            # blocks/segments are in chronological order without gaps and
            # subsequent blocks/segments those beg < farthest_segment_end
            # read so far are overlaping
            if seg['beg'] < last_end:

                if not overlap:
                    # new overlap interval
                    overlap = True
                    overlap_beg  = seg['beg']

                # There may be gaps between overlaping segments;
                # (separate overlaping intervals then)
                # TODO: maybe increase gap size (larger and less intermittent overlaps)
                elif abs(seg['beg'] - overlap_end) > 200:
                    overlaps.append((overlap_beg, overlap_end))
                    if args.debug:
                        print("      * distinct overlap {0} - {1}"
                              .format(overlap_beg, overlap_end))
                    overlap_beg  = seg['beg']

                overlap_end = seg['end']

                if args.verbose:
                    print("INFO: overlaping segment {0:.2f} - {1:.2f}"
                        .format(seg['beg']/1000, seg['end']/1000))


            else:
                last_end = seg['end']
                if overlap:
                    # end of overlaping
                    overlap = False
                    overlaps.append((overlap_beg, overlap_end))
                    if args.debug:
                        print("-------------------------------------------")
                        print("      End of overlaping; {0} - {1}\n"
                              .format(overlap_beg, overlap_end))

        else:
            print("WARN: Line '{0}' doesn't match segment format".format(lineno, line))


    
class Stats:
    segs = {}
    gaps = {}

    total_len = 0
    total_seg = 0
    total_seg_single = 0
    total_seg_joined = 0
    total_seg_combof = 0

    @classmethod
    def add(cls, tier, segments):
        if tier not in cls.segs:
            cls.segs[tier] = {
                'data': [],

                'single_cnt': 0,
                'joined_cnt': 0,
                'combof_cnt': 0,
                
                'min_len': 0,
                'max_len': 0,
                'avg_len': 0,
                'sum_len': 0,
                'len_quarts': []
            }


        if tier not in cls.gaps:
            cls.gaps[tier] = {
                'data': [],
                'min_len': 0,
                'max_len': 0,
                'avg_len': 0,
                'sum_len': 0,
                'len_quarts': []
            }

        cls.segs[tier]['data'].append(segments)


     
    @classmethod
    def collect(cls):
        
        for tier in cls.segs:
            seg_lens = []
            gap_lens = cls.gaps[tier]['data']
            
            for segs in cls.segs[tier]['data']:
                if len(segs) == 1:
                    cls.segs[tier]['single_cnt'] += 1
                elif len(segs) > 1:
                    cls.segs[tier]['joined_cnt'] += 1
                    cls.segs[tier]['combof_cnt'] += len(segs)
                
                beg = segs[0][0]  # first beg
                end = segs[-1][1] # last end
                seg_lens.append(end - beg)

                for (a, b) in pairwise(segs):
                    this_end = a[1]
                    next_beg = b[0]
                    gap_len = (next_beg - this_end)
                    if gap_len > 0:
                        gap_lens.append(gap_len)

            if seg_lens:    
                cls.segs[tier]['min_len'] = min(seg_lens)
                cls.segs[tier]['max_len'] = max(seg_lens)
                cls.segs[tier]['avg_len'] = stats.mean(seg_lens)
                cls.segs[tier]['sum_len'] = sum(seg_lens)
                if len(seg_lens) > 4:
                    cls.segs[tier]['len_quarts'] = stats.quantiles(
                        seg_lens, method='inclusive')

            if gap_lens:
                cls.gaps[tier]['min_len'] = min(gap_lens)
                cls.gaps[tier]['max_len'] = max(gap_lens)
                cls.gaps[tier]['avg_len'] = stats.mean(gap_lens)
                cls.gaps[tier]['sum_len'] = sum(gap_lens)
                if len(gap_lens) > 4:
                    cls.gaps[tier]['len_quarts'] = stats.quantiles(
                        gap_lens, method='inclusive')

            cls.total_len += cls.segs[tier]['sum_len']
            cls.total_seg += len(cls.segs[tier]['data'])
            cls.total_seg_single += cls.segs[tier]['single_cnt']
            cls.total_seg_joined += cls.segs[tier]['joined_cnt']
            cls.total_seg_combof += cls.segs[tier]['combof_cnt']

                
    
class Segment:
            
    def __init__(self, beg, end, text, sid=None):
        self.beg = beg
        self.text =  ''

        if sid:
            self.tier = sid
        else:
            self.tier = self.get_tier_name(beg, end)
        self.segments = []

        if not args.allow_overlength:
            self.max_length = args.max_length
        else:
            self.max_length = args.ultimate_length
        
        self.append(beg, end, text)


    @classmethod
    def overlaping(cls, beg, end):
        for (o_beg, o_end) in overlaps:
            if ((beg >= o_beg and beg <= o_end) or
                (end >= o_beg and end <= o_end)):
                return True

        return False

    def append(self, beg, end, text):

        # note: preserving original text/value
        self.segments.append((beg, end, text))
        
        self.end = end
        self.length = (self.end - self.beg)

        if args.debug:
            print("        {1} - {2}  seg. append: '{3}', comb. len.: {4}"
                  .format(self.tier, beg, end, text, self.length))


        text = self.process(text)
        
        if self.text == '':
            self.text =  text
        elif text:
            self.text += ' ' + text
        

    def can_join(self, beg, end, text, sid=None):

        consequent = ((beg - self.end) >= 0)
        combined_length = (end - self.beg)
        overlength = combined_length > args.max_length
        
        if not overlength:
            right_gap = ((beg - self.end) <= args.max_gap)
        else:
            right_gap = ((beg - self.end) <= args.max_overlength_gap)

        if sid:
            same_tier = (self.tier == sid)
        else:
            same_tier = (self.tier == self.get_tier_name(beg, end))


                        
        if (args.join_segments and
            same_tier and consequent and right_gap and
            combined_length > self.max_length):
            if args.verbose:
                print("NOTE: seg. {0:.2f} - {1:.2f} can still be joined, "
                      "but combined length {2:.2f} > {3:.2f} s"
                      .format(beg/1000, end/1000,
                              combined_length/1000, self.max_length/1000))
            
        return (args.join_segments and
                same_tier and consequent and right_gap and
                combined_length <= self.max_length)


    def process(self, text):
        text = text.strip()

        if args.strip_punctuation:
            text = text.translate(str.maketrans(
                {c: None for c in string.punctuation + '–„“'}))
            return re.sub(r'\s+', ' ', text).strip()

        # strip whitespace before punct. at the end ('-' - special case :-))
        text = re.sub(r'\s+(?=[^\s-]$)', '', text)

        # (some) combinations are separated by _
        text = text.replace('_', ' ')

        # webvtt segments may have newlines?
        text = text.replace("\n", ' ').replace("\r", '')

        # if (m := re.match(r'^\D*?(?P<num>\d+)\D*$', val)):
        #     print(text)

        return text

    
    def get_tier_name(self, beg, end):
        for blk in speech_blocks:
            # the first and the last segment of the speech block
            (_, b_beg, _, _) = speech_blocks[blk]['segs'][0]
            (_, _, b_end, _) = speech_blocks[blk]['segs'][-1]
            if beg >= b_beg and end <= b_end:
                return speech_blocks[blk]['sid']
             
        raise Exception("Can't find tier for time interval: {0} to {1}"
                        .format(beg, end))

    def added(self):
        Stats.add(self.tier, self.segments)

        if args.debug:
            print("{0}: ({2} - {3}) segment of {1} added, "
                  "length: {4}, text: '{5}'\n".format(self.tier,
                            len(self.segments), self.beg, self.end,
                            self.length, self.text))



def create_eaf():
    eaf = Eaf(author=args.author)
    
    if args.link_media:
        if args.orig_media:
            eaf.add_linked_file(args.link_media, ex_from=args.orig_media)
        else:
            eaf.add_linked_file(args.link_media)
        

    for sid in speech:
        eaf.add_tier(sid)

        if args.annotator:
            eaf.get_parameters_for_tier(sid)['ANNOTATOR'] = args.annotator

        if args.prefill_meta:
            eaf.get_parameters_for_tier(sid)['PARTICIPANT'] = args.prefill_meta

    return eaf

def last_setup(eaf):

    # Eaf() adds it
    eaf.remove_tier('default')
    
    if args.overlap_tier:
        eaf.add_tier('overlap')
        for (beg, end) in overlaps:
            if args.debug:
                print("Adding interval {0} - {1} to overlap tier"
                      .format(beg, end))
            eaf.add_annotation('overlap', beg, end)

    # TODO: CV for noise tier
    eaf.add_tier('noise')
    


def convert_lattice_to_eaf():

    eaf = create_eaf()
    
    segment = None
    
    for (sid, segs) in [(blk['sid'], blk['segs'])
                        for blk in speech_blocks.values()]:
        for (hyp, beg, end, val) in segs:

            seg = (beg, end, val)

            if args.skip_overlaping:
                if Segment.overlaping(beg, end):
                    # TODO: stats of overlaping segments (count, length)
                    if args.verbose:
                        print("INFO: {0} skipping overlaping segment "
                              "{1} - {2}".format(sid, beg, end))
                    continue
                
                        
            if segment is None:
                segment = Segment(*seg, sid)
            elif segment.can_join(*seg, sid):
                segment.append(*seg)
            else:
                # Segment can no longer be joined; add it to eaf
                eaf.add_annotation(segment.tier, segment.beg, segment.end, segment.text)
                segment.added()
                
                # New segment
                segment = Segment(*seg, sid)

        # Last segment
        if segment:
            eaf.add_annotation(segment.tier, segment.beg, segment.end, segment.text)
            segment.added()
            segment = None

    last_setup(eaf)    
    eaf.to_file(args.outfile)
        

def convert_webvtt_to_eaf(filename):

    eaf = create_eaf()

    segment = None
    
    for caption in webvtt.read(filename):
        beg = (caption.start_time.in_seconds() *
            1000 + caption.start_time.milliseconds)

        end = (caption.end_time.in_seconds() *
            1000 + caption.end_time.milliseconds)

        text = caption.text

        seg = (beg, end, text)
        
        # XXX: webvtt segments are already joined, so skipping whole chunk due to
        # overlap interval hit may be overkill
        if args.skip_overlaping:
            if Segment.overlaping(beg, end):
                # TODO: stats of overlaping segments (count, length)
                if args.verbose:
                    print("INFO: {0} skipping overlaping segment "
                          "{1} - {2}".format(sid, beg, end))
                continue
            

        if segment is None:
            segment = Segment(*seg)
        elif segment.can_join(*seg):
            segment.append(*seg)
        else:
            # Segment can no longer be joined; add it to eaf
            eaf.add_annotation(segment.tier, segment.beg, segment.end, segment.text)
            segment.added()
                
            # New segment
            segment = Segment(*seg)

    # Last segment
    if segment:
        eaf.add_annotation(segment.tier, segment.beg, segment.end, segment.text)
        segment.added()
        segment = None

    last_setup(eaf)
    eaf.to_file(args.outfile)
    

if args.webvtt:
    convert_webvtt_to_eaf(args.webvtt)
else:
    convert_lattice_to_eaf()

    
def quarts_in_seconds(quarts):
    """ Returns length quartiles in seconds """
    
    if not quarts:
        # (reasonable) data is not available 
        return ['-'] * 3
    else:
        return list(map(lambda q: q/1000, quarts))
    
    

Stats.collect()

table_totals = PrettyTable()
table_totals.align = "r"
table_totals.float_format = '0.2'

table_totals.field_names = [
    "Tiers", "Segs", "Single", "Joined", "Comb.of", "Duration"]
table_totals.add_row([
    len(Stats.segs), Stats.total_seg, Stats.total_seg_single,
    Stats.total_seg_joined, Stats.total_seg_combof,
    ms_to_ts(Stats.total_len)])

print("Total:")
print(table_totals.get_string())


table_by_tier = PrettyTable()
table_by_tier.field_names = [
    "Tier", "Segs", "Single", "Joined", "Comb.of", "Duration",
    "Min", "Avg", "Max", "Q25%", "Q50%", "Q75%"]

table_by_tier.align = "r"
table_by_tier.float_format = '0.2'
table_by_tier.sortby = "Duration"

for tier in Stats.segs:
    seg = Stats.segs[tier]
    table_by_tier.add_row([
        tier, len(seg['data']), seg['single_cnt'], seg['joined_cnt'],
        seg['combof_cnt'], ms_to_ts(seg['sum_len']), seg['min_len']/1000,
        seg['avg_len']/1000, seg['max_len']/1000,
        *quarts_in_seconds(seg['len_quarts'])])

print("By Tier:")
print(table_by_tier.get_string())

# print("\n\t{0} intervals/gaps between segments, total length: {1:.1f} s"
#       "\n\tLengths (ms): min: {2}, max: {3}, avg: {4:.0f}, "
#       "quart.: 25%: {5:.0f}, 50%: {6:.0f}, 75%: {7:.0f}"
#       .format(
#           len(Stats.gaps), Stats.gap_sum_len/1000,
#           Stats.gap_min_len, Stats.gap_max_len, Stats.gap_avg_len,
#           Stats.gap_quarts[0], Stats.gap_quarts[1], Stats.gap_quarts[2]))

