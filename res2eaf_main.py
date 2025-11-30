#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import string
import argparse

# from res2eaf_lib import Segment, Stats
import res2eaf_lib
# import webvtt
from pympi import Eaf
from pathlib import Path
from prettytable import PrettyTable


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



with open(args.lattice, 'r', encoding='utf-8') as lat_file:
    (speech, speech_blocks, overlaps, sid)=res2eaf_lib.parse_lat_content(lat_file, config=config)



if args.webvtt:
    res2eaf_lib.convert_webvtt_to_eaf(args.webvtt, sid=sid, config=config, speech_blocks=speech_blocks)
else:
    res2eaf_lib.convert_lattice_to_eaf(speech_blocks, config=config, speech=speech,overlaps=overlaps)


res2eaf_lib.Stats.collect()

table_totals = PrettyTable()
table_totals.align = "r"
table_totals.float_format = '0.2'

table_totals.field_names = [
    "Tiers", "Segs", "Single", "Joined", "Comb.of", "Duration"]
table_totals.add_row([
    len(res2eaf_lib.Stats.segs), res2eaf_lib.Stats.total_seg, res2eaf_lib.Stats.total_seg_single,
    res2eaf_lib.Stats.total_seg_joined, res2eaf_lib.Stats.total_seg_combof,
    res2eaf_lib.ms_to_ts(res2eaf_lib.Stats.total_len)])

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
        seg['combof_cnt'], res2eaf_lib.ms_to_ts(seg['sum_len']), seg['min_len']/1000,
        seg['avg_len']/1000, seg['max_len']/1000,
        *res2eaf_lib.quarts_in_seconds(seg['len_quarts'])])

print("By Tier:")
print(table_by_tier.get_string())

# print("\n\t{0} intervals/gaps between segments, total length: {1:.1f} s"
#       "\n\tLengths (ms): min: {2}, max: {3}, avg: {4:.0f}, "
#       "quart.: 25%: {5:.0f}, 50%: {6:.0f}, 75%: {7:.0f}"
#       .format(
#           len(Stats.gaps), Stats.gap_sum_len/1000,
#           Stats.gap_min_len, Stats.gap_max_len, Stats.gap_avg_len,
#           Stats.gap_quarts[0], Stats.gap_quarts[1], Stats.gap_quarts[2]))

