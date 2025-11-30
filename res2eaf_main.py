#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import string
import argparse

# from res2eaf_lib import Segment, Stats
import res2eaf_lib
import webvtt
from pympi import Eaf
from pathlib import Path
from prettytable import PrettyTable
from datetime import datetime, timezone

__author__ = "Laimonas Vėbra"
__copyright__ = "Copyright 2024"
__license__ = "BSD"
__version__ = "0.1d"


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
group.add_argument('--noise-tier', action='store_true',
                   help='Add noise tier')
group.add_argument('--noise-tier-cv', metavar='F_CSV', default='./noise_cv.csv',
                   help='CV (Controlled Vocabulary) file (tab(s) separated) '
                   'for noise tier')


args = parser.parse_args()
args_dict = vars(args)
print(f"**Parsed Arguments (Dict):** {args_dict}")


config = res2eaf_lib.Res2EafConfig.from_kwargs(**args_dict)

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
                          "    curr. sid: {2}, seq: {3}"
                          .format(last_sid, last_blk, sid, blk))

                if (last_sid and sid == last_sid):
                    print("WARN: same speaker '{0}' "
                          "consequitive blocks {1} and {2}"
                          .format(last_blk, blk))

                last_sid, last_blk = sid, blk

            else:
                print("WARN: Line:{0} '{1}' doesn't match header format"
                      .format(lineno, line))

        elif (m := segment.match(line)):
            # should never happen: speech line without (prior) header
            assert (sid is not None)

            # Skip fix.lattice.time inserted silence (TYLA) blocks/segments
            # they are empty anyway, but:
            # https://github.com/airenas/list/issues/1
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
                # TODO: maybe increase gap size (larger and less intermittent
                # overlaps)
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
            print("WARN: Line '{0}' doesn't match segment format"
                  .format(lineno, line))


if args.debug:
    print("Overlaps ({0}) before cleanup:".format(len(overlaps)))
    for overlap in overlaps:
        print(overlap)
    print()

# fix overlaps (remove inclusions, merge overlaping)
i = 0; _len = len(overlaps)
while i < _len:
    (beg1, end1) = overlaps[i]

    j = 0
    while j < _len:
        (beg2, end2) = overlaps[j]
        if i != j:
            if (beg1 >= beg2) and (end1 <= end2):
                if args.debug:
                    print("inclusive overlap: {0} in {1}; removing {0}"
                          .format(overlaps[i], overlaps[j]))
                overlaps.pop(i)
                i -= 1; _len -= 1
                break

            elif (beg1 >= beg2) and (beg1 <= end2) and (end1 > end2):
                if args.debug:
                    print("extending overlap: {0} by {1}; "
                          "removing {1}, extending: {0} -> {2}"
                          .format(overlaps[j], overlaps[i], (beg2, end1)))
                overlaps[j] = (beg2, end1)
                overlaps.pop(i)
                i -= 1; _len -= 1
                break
        j += 1
    i += 1

if args.debug:
    print("\nOverlaps ({0}) after cleanup:".format(len(overlaps)))
    for overlap in overlaps:
        print(overlap)
    print()






def read_noise_cv(noise_cv_path):
    noise_cv = []
    with open(noise_cv_path, encoding='utf-8') as f:
        header = f.readline()
        for line in f:
            if line.strip().startswith('#'):
                continue
            fields = re.split("\t+", line)
            assert len(fields) >= 2

            noise_cv.append({
                "name": fields[0].strip(),
                "desc": fields[1].strip()
            })
    return noise_cv


def add_overlap_tier(eaf):
    eaf.add_tier('overlap')
    for (beg, end) in overlaps:
        if args.debug:
            print("Adding interval {0} - {1} to overlap tier"
                  .format(beg, end))
        eaf.add_annotation('overlap', beg, end)

def add_noise_tier(eaf):
    eaf.add_linguistic_type('noise-lt')

    if args.noise_tier_cv:
        if not Path(args.noise_tier_cv).is_file():
            print("WARN: the specified noise CV file '{0}' does not exist",
                  args.noise_tier_cv)
        else:
            noise_cv = read_noise_cv(args.noise_tier_cv)
            if noise_cv:
                eaf.add_controlled_vocabulary('noise-cv')
                eaf.add_cv_description('noise-cv', 'lit',
                                       description='triukšmų žodynas')

                for i, cve in enumerate(noise_cv):
                    eaf.add_cv_entry('noise-cv', ('cveid' + str(i)),
                                     [(cve['name'], 'lit', cve['desc'])])

                    (eaf.get_parameters_for_linguistic_type('noise-lt')
                     ['CONTROLLED_VOCABULARY_REF']) = 'noise-cv'

    eaf.add_tier('noise', ling='noise-lt', language='und')


def create_eaf():
    eaf = Eaf(author=args.author)

    eaf.add_language('und',
                     'http://cdb.iso.org/lg/CDB-00130975-001',
                     'undetermined (und)')

    eaf.add_language('lit',
                     'http://cdb.iso.org/lg/CDB-00138562-001',
                     'Lithuanian (lit)')

    if args.link_media:
        if args.orig_media:
            eaf.add_linked_file(args.link_media, relpath=args.link_media,
                                ex_from=args.orig_media)
        else:
            eaf.add_linked_file(args.link_media, relpath=args.link_media)


    for sid in speech:
        eaf.add_tier(sid, language='lit')

        if args.annotator:
            eaf.get_parameters_for_tier(sid)['ANNOTATOR'] = args.annotator

        if args.prefill_meta:
            eaf.get_parameters_for_tier(sid)['PARTICIPANT'] = args.prefill_meta

    return eaf


def last_setup(eaf):
    eaf.remove_tier('default') # Eaf() adds it

    # remove empty tiers
    for tier in list(eaf.get_tier_names()):
        if (len(eaf.tiers[tier][0]) == 0):
            if args.verbose:
                print("INFO: removing empty tier '{0}'".format(tier))
            eaf.remove_tier(tier)

    if args.overlap_tier:
        add_overlap_tier(eaf)

    if args.noise_tier:
        add_noise_tier(eaf)



def convert_lattice_to_eaf(speech_blocks):

    eaf = create_eaf()

    segment = None

    for (sid, segs) in [(blk['sid'], blk['segs'])
                        for blk in speech_blocks.values()]:
        for (hyp, beg, end, val) in segs:

            seg = (beg, end, val)

            if args.skip_overlaping:
                if res2eaf_lib.Segment.overlaping(beg=beg, end=end, overlaps=overlaps):
                    # TODO: stats of overlaping segments (count, length)
                    if args.verbose:
                        print("INFO: {0} skipping overlaping segment "
                              "{1} - {2}".format(sid, beg, end))
                    continue

            if segment is None:
                segment = res2eaf_lib.Segment(beg=beg, end=end, text=val, config=config, speech_blocks =speech_blocks,sid=sid)
            elif segment.can_join(*seg, sid):
                segment.append(*seg)
            else:
                # Segment can no longer be joined; add it to eaf
                eaf.add_annotation(segment.tier, segment.beg, segment.end,
                                   segment.text)
                segment.added()

                # New segment
                segment = res2eaf_lib.Segment(beg=beg, end=end, text=val, config=config, speech_blocks =speech_blocks, sid=sid)

        # Last segment
        if segment:
            eaf.add_annotation(segment.tier, segment.beg, segment.end,
                               segment.text)
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

        # seg = (beg, end, text)

        # XXX: webvtt segments are already joined, so skipping whole chunk due
        # to overlap interval hit may be overkill
        if args.skip_overlaping:
            if res2eaf_lib.Segment.overlaping(beg=beg, end=end):
                # TODO: stats of overlaping segments (count, length)
                if args.verbose:
                    print("INFO: {0} skipping overlaping segment "
                          "{1} - {2}".format(sid, beg, end))
                continue

        if segment is None:
            segment = res2eaf_lib.Segment(beg=beg, end=end, text=text,speech_blocks =speech_blocks)
        elif segment.can_join(beg=beg, end=end, text=text):
            segment.append(beg=beg, end=end, text=text)
        else:
            # Segment can no longer be joined; add it to eaf
            eaf.add_annotation(segment.tier, segment.beg, segment.end,
                               segment.text)
            segment.added()

            # New segment
            segment = res2eaf_lib.Segment(beg=beg, end=end, text=text, speech_blocks =speech_blocks)

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
    convert_lattice_to_eaf(speech_blocks)


def quarts_in_seconds(quarts):
    """ Returns length quartiles in seconds """

    if not quarts:
        # (reasonable) data is not available
        return ['-'] * 3
    else:
        return list(map(lambda q: q/1000, quarts))



res2eaf_lib.Stats.collect()

table_totals = PrettyTable()
table_totals.align = "r"
table_totals.float_format = '0.2'

table_totals.field_names = [
    "Tiers", "Segs", "Single", "Joined", "Comb.of", "Duration"]
table_totals.add_row([
    len(res2eaf_lib.Stats.segs), res2eaf_lib.Stats.total_seg, res2eaf_lib.Stats.total_seg_single,
    res2eaf_lib.Stats.total_seg_joined, res2eaf_lib.Stats.total_seg_combof,
    ms_to_ts(res2eaf_lib.Stats.total_len)])

print("Total:")
print(table_totals.get_string())


table_by_tier = PrettyTable()
table_by_tier.field_names = [
    "Tier", "Segs", "Single", "Joined", "Comb.of", "Duration",
    "Min", "Avg", "Max", "Q25%", "Q50%", "Q75%"]

table_by_tier.align = "r"
table_by_tier.float_format = '0.2'
table_by_tier.sortby = "Duration"

for tier in res2eaf_lib.Stats.segs:
    seg = res2eaf_lib.Stats.segs[tier]
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

