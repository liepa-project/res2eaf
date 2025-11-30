import statistics as stats
from itertools import pairwise
# from datetime import timedelta
# from typing import Dict, List, Tuple, Any, Optional
from typing import  Optional
from dataclasses import dataclass, fields
import re

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
