# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
# Ref: https://stackoverflow.com/questions/37790429/seaborn-heatmap-using-pandas-dataframe
# Ref: https://python-graph-gallery.com/4-add-title-and-axis-label/
# Ref: https://python-graph-gallery.com/92-control-color-in-seaborn-heatmaps/
# Ref: https://stackoverflow.com/questions/34232073/seaborn-heatmap-y-axis-reverse-order


import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
import math

import sys

print ('[SYNTAX] python aggregate_embedding.py <input: embedding> <output: heatmap filename>')


fnamein = sys.argv[1]
fnameout = sys.argv[2]

df = pd.read_csv (fnamein, sep ='\t')
df_agg = df.groupby(['Dim_1', 'Dim_2'], as_index = False).agg({'ID': 'count'})
df_agg['log_n'] = df_agg.apply(lambda x: math.log10(x['ID']), axis = 1)

# Note that Dim_1 and Dim_2 seem to be reversed in the next line, the reason is to get a heatmap which better reflect the oritnation of HCASE embedding plots. -- GZK
df_pivot = df_agg.pivot(index='Dim_2', columns='Dim_1', values='log_n')
df_pivot = df_pivot.reset_index(drop = True)
df_pivot = df_pivot.fillna(0)


plt.figure()


ax = sns.heatmap(df_pivot, annot=False, cmap='Blues')
plt.xlabel('Dim 1')
plt.ylabel('Dim 2')
ax.invert_yaxis()

plt.savefig (fnameout, dpi=300)


#plt.show()
