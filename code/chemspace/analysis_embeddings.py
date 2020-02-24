# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
#
#
# Ref: https://github.com/matplotlib/matplotlib/issues/3466/
# Ref: https://matplotlib.org/3.1.1/gallery/shapes_and_collections/scatter.html
# Ref: https://stackoverflow.com/questions/17682216/scatter-plot-and-color-mapping-in-python
# Ref: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
# Ref: https://python-graph-gallery.com/106-seaborn-style-on-matplotlib-plot/
# Ref: https://python-graph-gallery.com/122-multiple-lines-chart/
# Ref: https://python-graph-gallery.com/125-small-multiples-for-line-chart/
# Ref: https://seaborn.pydata.org/examples/many_facets.html
# Ref: https://matplotlib.org/3.1.1/api/markers_api.html#module-matplotlib.markers
# Ref: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.scatter.html
# Ref: https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/subplots_demo.html
# Ref: https://matplotlib.org/3.1.1/tutorials/colors/colors.html
# Ref: https://stackoverflow.com/questions/14827650/pyplot-scatter-plot-marker-size
# Ref: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.title.html
# Ref: https://stackoverflow.com/questions/31686530/matplotlib-generate-a-new-graph-in-a-new-window-for-subsequent-program-runs/31686783
# Ref: https://stackoverflow.com/questions/39870642/matplotlib-how-to-plot-a-high-resolution-graph
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html
# Ref: https://iwatobipen.wordpress.com/2017/11/03/draw-high-quality-molecular-image-in-rdkit-rdkit/
# Ref: https://stackoverflow.com/questions/14432557/matplotlib-scatter-plot-with-different-text-at-each-data-point
# Ref: https://python-graph-gallery.com/197-available-color-palettes-with-matplotlib/
# Ref: https://stackoverflow.com/questions/57568311/matplotlib-scatter-issue-with-python-3-x
# Ref: https://www.science-emergence.com/Articles/How-to-create-a-scatter-plot-with-several-colors-in-matplotlib-/
# Ref: https://www.pluralsight.com/guides/choosing-color-palettes
# Ref: https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf
# Ref: https://htmlcolorcodes.com/color-picker/



# #1179B0
# #F58C30
# #74BB5A
# #BC412C
# #795B9A
# #764A0C
# #D37DB5
# #7A7A7A
# #B8C449



import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

import rdkit
from rdkit import Chem
#from rdkit.Chem import Draw
#from rdkit.Chem.Draw import DrawingOptions

import pandas as pd


def create_plot (df, outputfile, title):
    print(df.columns)

    plt.figure()

    plt.scatter(df['Dim_1'].values, df['Dim_2'].values, alpha = 0.3, marker = 'o', s = 10)

    plt.savefig (outputfile, dpi=300)

def create_dual_plot (df1, df2, outputfile, title):

    plt.figure()

    plt.scatter(df1['Dim_1'].values, df1['Dim_2'].values, c = 'blue', alpha = 0.3, marker = 'o', s = 10)
    plt.scatter(df2['Dim_1'].values, df2['Dim_2'].values, c = 'yellow', alpha = 0.3, marker = 'o', s = 10)


    plt.savefig (outputfile, dpi=300)




df_2 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_5_dim_2.tab', sep = '\t')
df_6 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_6_dim_2.tab', sep = '\t')
df_7 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_7_dim_2.tab', sep = '\t')
df_8 = pd.read_csv ('../../data/canvass_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')

canvass_chembl_bms_dfs = []

canvass_chembl_bms_dfs.append(df_2)
canvass_chembl_bms_dfs.append(df_3)
canvass_chembl_bms_dfs.append(df_4)
canvass_chembl_bms_dfs.append(df_5)
canvass_chembl_bms_dfs.append(df_6)
canvass_chembl_bms_dfs.append(df_7)
canvass_chembl_bms_dfs.append(df_8)








create_plot (df_2, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_2_dim_2.png', 'CANVASS\nPHC order 2, ChEMBL 24_1 BMSs')
create_plot (df_3, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_3_dim_2.png', 'CANVASS\nPHC order 3, ChEMBL 24_1 BMSs')
create_plot (df_4, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_4_dim_2.png', 'CANVASS\nPHC order 4, ChEMBL 24_1 BMSs')
create_plot (df_5, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_5_dim_2.png', 'CANVASS\nPHC order 5, ChEMBL 24_1 BMSs')
create_plot (df_6, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_6_dim_2.png', 'CANVASS\nPHC order 6, ChEMBL 24_1 BMSs')
create_plot (df_7, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_7_dim_2.png', 'CANVASS\nPHC order 7, ChEMBL 24_1 BMSs')
create_plot (df_8, '../../plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_8_dim_2.png', 'CANVASS\nPHC order 8, ChEMBL 24_1 BMSs')


df_2['phc_order'] = 2
df_3['phc_order'] = 3
df_4['phc_order'] = 4
df_5['phc_order'] = 5
df_6['phc_order'] = 6
df_7['phc_order'] = 7
df_8['phc_order'] = 8


df_coords = df_2.append(df_3, ignore_index = True)
df_coords = df_coords.append(df_4, ignore_index = True)
df_coords = df_coords.append(df_5, ignore_index = True)
df_coords = df_coords.append(df_6, ignore_index = True)
df_coords = df_coords.append(df_7, ignore_index = True)
df_coords = df_coords.append(df_8, ignore_index = True)

df_coords.to_csv ('../../data/canvass_all_embedding_coords_chembl_24_1_bms.tab', sep = '\t', index = False)



########## Approved Drugs DrugBank into NatProd Space

df_2 = pd.read_csv ('../../data/app_drugbank_into_hc_natprod_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/app_drugbank_into_hc_natprod_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/app_drugbank_into_hc_natprod_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/app_drugbank_into_hc_natprod_bms_ord_5_dim_2.tab', sep = '\t')
#df_6 = pd.read_csv ('../../data/app_drugbank_into_hc_natprod_bms_ord_6_dim_2.tab', sep = '\t')

drugbank_natprod_dfs = []

drugbank_natprod_dfs.append (df_2)
drugbank_natprod_dfs.append (df_3)
drugbank_natprod_dfs.append (df_4)
drugbank_natprod_dfs.append (df_5)
#drugbank_natprod_dfs.append (df_6)





create_plot (df_2, '../../plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_2_dim_2.png', 'Approved Drugs DugBank\nPHC order 2, NatProd BMSs')
create_plot (df_3, '../../plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_3_dim_2.png', 'Approved Drugs DugBank\nPHC order 3, NatProd BMSs')
create_plot (df_4, '../../plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_4_dim_2.png', 'Approved Drugs DugBank\nPHC order 4, NatProd BMSs')
create_plot (df_5, '../../plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_5_dim_2.png', 'Approved Drugs DugBank\nPHC order 5, NatProd BMSs')
#create_plot (df_6, '../../plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_6_dim_2.png', 'Approved Drugs DugBank\nPHC order 6, NatProd BMSs')



df_2['phc_order'] = 2
df_3['phc_order'] = 3
df_4['phc_order'] = 4
df_5['phc_order'] = 5
#df_6['phc_order'] = 6


df_coords = df_2.append(df_3, ignore_index = True)
df_coords = df_coords.append(df_4, ignore_index = True)
df_coords = df_coords.append(df_5, ignore_index = True)
#df_coords = df_coords.append(df_6, ignore_index = True)

df_coords.to_csv ('../../data/app_drugbank_all_embedding_coords_natprod_bms.tab', sep = '\t', index = False)

########## CANVASS into NatProd Space

df_2 = pd.read_csv ('../../data/canvass_into_hc_natprod_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/canvass_into_hc_natprod_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/canvass_into_hc_natprod_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/canvass_into_hc_natprod_bms_ord_5_dim_2.tab', sep = '\t')
#df_6 = pd.read_csv ('../../data/canvass_into_hc_natprod_bms_ord_6_dim_2.tab', sep = '\t')

canvass_natprod_dfs = []

canvass_natprod_dfs.append(df_2)
canvass_natprod_dfs.append(df_3)
canvass_natprod_dfs.append(df_4)
canvass_natprod_dfs.append(df_5)
#canvass_natprod_dfs.append(df_6)





create_plot (df_2, '../../plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_2_dim_2.png', 'CANVASS\nPHC order 2, NatProd BMSs')
create_plot (df_3, '../../plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_3_dim_2.png', 'CANVASS\nPHC order 3, NatProd BMSs')
create_plot (df_4, '../../plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_4_dim_2.png', 'CANVASS\nPHC order 4, NatProd BMSs')
create_plot (df_5, '../../plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_5_dim_2.png', 'CANVASS\nPHC order 5, NatProd BMSs')
#create_plot (df_6, '../../plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_6_dim_2.png', 'CANVASS\nPHC order 6, NatProd BMSs')



df_2['phc_order'] = 2
df_3['phc_order'] = 3
df_4['phc_order'] = 4
df_5['phc_order'] = 5
#df_6['phc_order'] = 6


df_coords = df_2.append(df_3, ignore_index = True)
df_coords = df_coords.append(df_4, ignore_index = True)
df_coords = df_coords.append(df_5, ignore_index = True)
#df_coords = df_coords.append(df_6, ignore_index = True)

df_coords.to_csv ('../../data/canvass_all_embedding_coords_natprod_bms.tab', sep = '\t', index = False)



######## Overlaid plots ##########

####### Importing Approved Drugs DrugBank embedded into ChEMBL 24.1 BMS HC-space data

df_2 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2.tab', sep = '\t')
df_6 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2.tab', sep = '\t')
df_7 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2.tab', sep = '\t')
df_8 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')

app_drugbank_chembl_bms_dfs = []

app_drugbank_chembl_bms_dfs.append(df_2)
app_drugbank_chembl_bms_dfs.append(df_3)
app_drugbank_chembl_bms_dfs.append(df_4)
app_drugbank_chembl_bms_dfs.append(df_5)
app_drugbank_chembl_bms_dfs.append(df_6)
app_drugbank_chembl_bms_dfs.append(df_7)
app_drugbank_chembl_bms_dfs.append(df_8)







outfile = ''
title = ''
for i in range(7):
    outfile = '../../plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms'
    outfile += '_ord_' + str(i + 2) + '_dim_2.png'
    title = 'CANVASS and Approved Drugs DrugBank\nin ChEMBL 24.1 BMSs HC Space, order: ' + str(i + 2) 
    df1 = canvass_chembl_bms_dfs[i]
    df2 = app_drugbank_chembl_bms_dfs[i]
    create_dual_plot (df1, df2, outfile, title)
 






outfile = ''
title = ''
for i in range(5):
    outfile = '../../plots/comparative_embeddings/dual_canvass_drugbank_natrpod'
    outfile += '_ord_' + str(i + 2) + '_dim_2.png'
    title = 'CANVASS and Approved Drugs DrugBank\nin NatProd BMSs HC Space, order: ' + str(i + 2) 
    df1 = canvass_natprod_dfs[i]
    df2 = drugbank_natprod_dfs[i]
    create_dual_plot (df1, df2, outfile, title)
    
 



print ('[Done.]')

