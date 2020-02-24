# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Aim: Cherry-pick ChEBML 24_1 Scaffold-Key ordered scaffolds for analysis.
#
#

import pandas as pd

sk_ordered_scaffolds_fname = '../../data/hc_space.tab'
fnameout = '../../data/cherrypicked_scaffolds.tab'

df = pd.read_csv (sk_ordered_scaffolds_fname, sep = '\t')

major_markers = [5000, 15000, 16000, 25000, 26000, 35000, 44000, 45000, 55000]

selection = []
df_res = pd.DataFrame()
first = True

color_idx = 1
for marker in major_markers:
    selection = []
    for i in range (marker - 50, marker):
        selection.append(i)
    for i in range (marker, marker + 51):
        selection.append(i)
    df_cp = df[df['order'].isin(selection)].copy()
    df_cp['color'] = color_idx
    color_idx += 1
    if first:
        df_res = df_cp
        first = False
    else:
        df_res = df_res.append(df_cp, ignore_index = True)

df_res = df_res.reset_index (drop = True)

df_res.to_csv(fnameout, sep ='\t', index = False)

print('[Done.]')

