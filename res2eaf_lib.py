import statistics as stats
from itertools import pairwise
from datetime import datetime, timezone
# from datetime import timedelta
# from typing import Dict, List, Tuple, Any, Optional
from typing import  Optional
from dataclasses import dataclass, fields
import re

import webvtt
from pympi import Eaf
from pathlib import Path
from io import StringIO

@dataclass
class Res2EafConfig:
    """Holds configuration parameters for LAT to EAF conversion."""
    
    # Input/Output (now handled by caller, but kept for context if needed)
    lattice_path: Optional[str] = None # The path is needed for media linking/naming
    webvtt_path: Optional[str] = None 
    
    # General
    verbose: bool = False
    debug: bool = False
    skip_overlaping: bool = False

    # Joining segments
    join_segments: bool = False
    max_gap: int = 100
    max_length: int = 5000
    allow_overlength: bool = False
    ultimate_length: int = 7000
    max_overlength_gap: int = 50

    # Text processing
    strip_punctuation: bool = False

    # EAF creation
    author: str = 'res2eaf_refactored:0.1d' # Default from your original script
    link_media: Optional[str] = None
    orig_media: Optional[str] = None
    annotator: Optional[str] = None
    prefill_meta: Optional[str] = None
    overlap_tier: bool = False
    noise_tier: bool = False
    noise_tier_cv: str = './noise_cv.csv'
    outfile: Optional[str] = None


    @classmethod
    def from_kwargs(cls, **kwargs):
        names = {f.name for f in fields(cls)}
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in names}
        return cls(**filtered_kwargs)
    
    


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

    def __init__(self, beg, end, text, config: Res2EafConfig, speech_blocks, sid=None):
        self.beg = beg
        self.text =  ''
        self.config = config
        self.speech_blocks =speech_blocks
        if sid:
            self.tier = sid
        else:
            self.tier = self.get_tier_name(beg, end)
        self.segments = []

        if not self.config.allow_overlength:
            self.max_length = self.config.max_length
        else:
            self.max_length = self.config.ultimate_length

        self.append(beg, end, text)


    @classmethod
    def overlaping(cls, beg, end, overlaps):
        for (o_beg, o_end) in overlaps:
            if (
                    # hit: beg or end in overlap interval
                    (beg >= o_beg and beg <= o_end) or
                    (end >= o_beg and end <= o_end)
                    or
                    # segment includes overlap interval
                    (beg < o_beg and end > o_end)
            ):
                return True

        return False

    def append(self, beg, end, text):

        # note: preserving original text/value
        self.segments.append((beg, end, text))

        self.end = end
        self.length = (self.end - self.beg)

        if self.config.debug:
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
        overlength = combined_length > self.config.max_length

        if not overlength:
            right_gap = ((beg - self.end) <= self.config.max_gap)
        else:
            right_gap = ((beg - self.end) <= self.config.max_overlength_gap)

        if sid:
            same_tier = (self.tier == sid)
        else:
            same_tier = (self.tier == self.get_tier_name(beg, end))



        if (self.config.join_segments and
            same_tier and consequent and right_gap and
            combined_length > self.max_length):
            if self.config.verbose:
                print("NOTE: seg. {0:.2f} - {1:.2f} can still be joined, "
                      "but combined length {2:.2f} > {3:.2f} s"
                      .format(beg/1000, end/1000,
                              combined_length/1000, self.max_length/1000))

        return (self.config.join_segments and
                same_tier and consequent and right_gap and
                combined_length <= self.max_length)


    def process(self, text):
        text = text.strip()

        if self.config.strip_punctuation:
            text = text.translate(str.maketrans(
                {c: None for c in self.config.punctuation + '–„“'}))
            return re.sub(r'\s+', ' ', text).strip()

        # strip whitespace before punct. at the end ('-' - special case :-))
        text = re.sub(r'\s+(?=[^\s-]$)', '', text)

        # (some) combinations are separated by _
        text = text.replace('_', ' ')

        # webvtt segments may have newlines?
        text = text.replace("\n", ' ').replace("\r", '')

        return text


    def get_tier_name(self, beg, end):
        for blk in self.speech_blocks:
            # the first and the last segment of the speech block
            (_, b_beg, _, _) = self.speech_blocks[blk]['segs'][0]
            (_, _, b_end, _) = self.speech_blocks[blk]['segs'][-1]
            if beg >= b_beg and end <= b_end:
                return self.speech_blocks[blk]['sid']

        raise Exception("Can't find tier for time interval: {0} to {1}"
                        .format(beg, end))

    def added(self):
        Stats.add(self.tier, self.segments)

        if self.config.debug:
            print("{0}: ({2} - {3}) segment of {1} added, "
                  "length: {4}, text: '{5}'\n".format(self.tier,
                            len(self.segments), self.beg, self.end,
                            self.length, self.text))


def convert_webvtt_to_eaf(filename, sid, config: Res2EafConfig,  speech_blocks):

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
        if config.skip_overlaping:
            if Segment.overlaping(beg=beg, end=end):
                # TODO: stats of overlaping segments (count, length)
                if config.verbose:
                    print("INFO: {0} skipping overlaping segment "
                          "{1} - {2}".format(sid, beg, end))
                continue

        if segment is None:
            segment = Segment(beg=beg, end=end, text=text, config=config, speech_blocks=speech_blocks)
        elif segment.can_join(beg=beg, end=end, text=text):
            segment.append(beg=beg, end=end, text=text)
        else:
            # Segment can no longer be joined; add it to eaf
            eaf.add_annotation(segment.tier, segment.beg, segment.end,
                               segment.text)
            segment.added()

            # New segment
            segment = Segment(beg=beg, end=end, text=text, speech_blocks=speech_blocks)

    # Last segment
    if segment:
        eaf.add_annotation(segment.tier, segment.beg, segment.end, segment.text)
        segment.added()
        segment = None

    last_setup(eaf)
    eaf.to_file(config.outfile)


    

def create_eaf(config: Res2EafConfig, speech):
    eaf = Eaf(author=config.author)

    eaf.add_language('und',
                     'http://cdb.iso.org/lg/CDB-00130975-001',
                     'undetermined (und)')

    eaf.add_language('lit',
                     'http://cdb.iso.org/lg/CDB-00138562-001',
                     'Lithuanian (lit)')

    if config.link_media:
        if config.orig_media:
            eaf.add_linked_file(config.link_media, relpath=config.link_media,
                                ex_from=config.orig_media)
        else:
            eaf.add_linked_file(config.link_media, relpath=config.link_media)


    for sid in speech:
        eaf.add_tier(sid, language='lit')

        if config.annotator:
            eaf.get_parameters_for_tier(sid)['ANNOTATOR'] = config.annotator

        if config.prefill_meta:
            eaf.get_parameters_for_tier(sid)['PARTICIPANT'] = config.prefill_meta

    return eaf


def convert_lattice_to_eaf(speech_blocks, config: Res2EafConfig, speech, overlaps):
    """
    ***************************************************************
    """

    eaf = create_eaf(config=config,speech=speech)

    segment = None

    for (sid, segs) in [(blk['sid'], blk['segs'])
                        for blk in speech_blocks.values()]:
        for (hyp, beg, end, val) in segs:

            seg = (beg, end, val)

            if config.skip_overlaping:
                if Segment.overlaping(beg=beg, end=end, overlaps=overlaps):
                    # TODO: stats of overlaping segments (count, length)
                    if config.verbose:
                        print("INFO: {0} skipping overlaping segment "
                              "{1} - {2}".format(sid, beg, end))
                    continue

            if segment is None:
                segment = Segment(beg=beg, end=end, text=val, config=config, speech_blocks =speech_blocks,sid=sid)
            elif segment.can_join(*seg, sid):
                segment.append(*seg)
            else:
                # Segment can no longer be joined; add it to eaf
                eaf.add_annotation(segment.tier, segment.beg, segment.end,
                                   segment.text)
                segment.added()

                # New segment
                segment = Segment(beg=beg, end=end, text=val, config=config, speech_blocks =speech_blocks, sid=sid)

        # Last segment
        if segment:
            eaf.add_annotation(segment.tier, segment.beg, segment.end,
                               segment.text)
            segment.added()
            segment = None

    last_setup(eaf, config=config, overlaps=overlaps)
    eaf.to_file(config.outfile)


def last_setup(eaf, config: Res2EafConfig, overlaps):
    """
    ****************************************
    """
    eaf.remove_tier('default') # Eaf() adds it

    # remove empty tiers
    for tier in list(eaf.get_tier_names()):
        if (len(eaf.tiers[tier][0]) == 0):
            if config.verbose:
                print("INFO: removing empty tier '{0}'".format(tier))
            eaf.remove_tier(tier)

    if config.overlap_tier:
        add_overlap_tier(eaf, config=config, overlaps=overlaps)

    if config.noise_tier:
        add_noise_tier(eaf, config=config)


def add_overlap_tier(eaf, config: Res2EafConfig, overlaps):
    eaf.add_tier('overlap')
    for (beg, end) in overlaps:
        if config.debug:
            print("Adding interval {0} - {1} to overlap tier"
                  .format(beg, end))
        eaf.add_annotation('overlap', beg, end)


def add_noise_tier(eaf,config: Res2EafConfig):
    """
    *****************************************
    """
    eaf.add_linguistic_type('noise-lt')

    if config.noise_tier_cv:
        if not Path(config.noise_tier_cv).is_file():
            print("WARN: the specified noise CV file '{0}' does not exist",
                  config.noise_tier_cv)
        else:
            noise_cv = read_noise_cv(config.noise_tier_cv)
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


def read_noise_cv(noise_cv_path):
    """
    *************************************
    """
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
    
def quarts_in_seconds(quarts):
    """ Returns length quartiles in seconds """

    if not quarts:
        # (reasonable) data is not available
        return ['-'] * 3
    else:
        return list(map(lambda q: q/1000, quarts))

def parse_lat_content_initiate(lat_content:StringIO, config: Res2EafConfig):
    """
    ****************************************
    """
    speech = {}
    speech_blocks = {}
    overlaps = []


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

    for line in lat_content:
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
                    if config.debug:
                        print("      * distinct overlap {0} - {1}"
                              .format(overlap_beg, overlap_end))
                    overlap_beg  = seg['beg']

                overlap_end = seg['end']

                if config.verbose:
                    print("INFO: overlaping segment {0:.2f} - {1:.2f}"
                        .format(seg['beg']/1000, seg['end']/1000))


            else:
                last_end = seg['end']
                if overlap:
                    # end of overlaping
                    overlap = False
                    overlaps.append((overlap_beg, overlap_end))
                    if config.debug:
                        print("-------------------------------------------")
                        print("      End of overlaping; {0} - {1}\n"
                              .format(overlap_beg, overlap_end))

        else:
            print("WARN: Line '{0}' doesn't match segment format"
                  .format(lineno, line))
    if config.debug:
        print("Overlaps ({0}) before cleanup:".format(len(overlaps)))
        for overlap in overlaps:
            print(overlap)
        print()
    return (speech, speech_blocks, overlaps, sid)

def fix_overlaps(overlaps, config:Res2EafConfig):
    # fix overlaps (remove inclusions, merge overlaping)
    i = 0; _len = len(overlaps)
    while i < _len:
        (beg1, end1) = overlaps[i]

        j = 0
        while j < _len:
            (beg2, end2) = overlaps[j]
            if i != j:
                if (beg1 >= beg2) and (end1 <= end2):
                    if config.debug:
                        print("inclusive overlap: {0} in {1}; removing {0}"
                            .format(overlaps[i], overlaps[j]))
                    overlaps.pop(i)
                    i -= 1; _len -= 1
                    break

                elif (beg1 >= beg2) and (beg1 <= end2) and (end1 > end2):
                    if config.debug:
                        print("extending overlap: {0} by {1}; "
                            "removing {1}, extending: {0} -> {2}"
                            .format(overlaps[j], overlaps[i], (beg2, end1)))
                    overlaps[j] = (beg2, end1)
                    overlaps.pop(i)
                    i -= 1; _len -= 1
                    break
            j += 1
        i += 1

    if config.debug:
        print("\nOverlaps ({0}) after cleanup:".format(len(overlaps)))
        for overlap in overlaps:
            print(overlap)
        print()
    return overlaps

def parse_lat_content(lat_content:StringIO, config: Res2EafConfig):
    (speech, speech_blocks, overlaps, sid)=parse_lat_content_initiate(lat_content, config=config)
    overlaps=fix_overlaps(overlaps=overlaps, config=config)
    return (speech, speech_blocks, overlaps, sid)