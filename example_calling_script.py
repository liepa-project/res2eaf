# In your calling_script.py
from io import StringIO
import res2eaf_lib 

config=res2eaf_lib.Res2EafConfig()
config.join_segments=True
output_file_path = 'test/test_example.eaf'

# lat_file_path = 'test/test.lat'
# with open(lat_file_path, 'r', encoding='utf-8') as lat_file:
    # (speech, speech_blocks, overlaps, sid)=res2eaf_lib.parse_lat_content(lat_file, config=config)
# eaf_obj=res2eaf_lib.convert_lattice_to_eaf_obj(speech_blocks, config=config, speech=speech,overlaps=overlaps)
# eaf_obj.to_file(output_file_path)


str_lat="""
# 1 S0001
1 0.03  0.83 skvernelis
1 0.85  1.19 atskleidžia
1 1.21  1.55 planą
1 1.61  2.10 chuliganą
1 2.12  2.86 baltarusijai
1 3.04  3.52 astravo
1 3.54  3.94 atominę
1 4.02  4.62 pakeisti
1 4.78  5.20 dujine.

# 2 S0002
1 5.36  6.04 atominės
1 6.06  7.29 elektrinės
1 7.33  8.39 konversija
1 8.41  8.43 į
1 8.45  8.57 dujinę
1 9.05  9.59 praktiškai
1 9.61  10.93 neįmanoma.

# 3 S0001
1 10.97  11.05 apie
1 11.07  11.65 pirmadienį
1 11.68  11.82 daugiau
1 11.84  12.48 įprastu
1 12.50  12.94 panoramos
1 12.96  22.06 laiku.

# 4 UNKOWN
1 22.08  22.16 pusę

# 5 S0004
1 22.18  22.34 devintos
1 22.38  23.00 antradienio
1 23.02  23.36 ryto
1 23.56  23.94 allegro
1 24.16  24.72 įsijungti
1 24.74  24.76 į
1 24.78  25.32 gyvenimą
1 25.34  26.00 naujuoju
1 26.09  27.25 albumu
1 27.39  28.11 kviečiantys
1 28.21  28.71 subtilu
1 28.73  29.31 zė
"""

(speech, speech_blocks, overlaps, sid)=res2eaf_lib.parse_lat_content(StringIO(str_lat), config=config)
eaf_obj=res2eaf_lib.convert_lattice_to_eaf_obj(speech_blocks, config=config, speech=speech,overlaps=overlaps)
eaf_obj.to_file(output_file_path)