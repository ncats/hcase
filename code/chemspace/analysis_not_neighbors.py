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



import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

import rdkit
from rdkit import Chem
#from rdkit.Chem import Draw
#from rdkit.Chem.Draw import DrawingOptions

import knn
from knn import get_mol, get_fingerprint, get_Tanimoto
import pandas as pd

"""
def plot_knn (df, outputfile, title):
    #output_plot_fname = '../../plots/' + outfile_stem + '.png'
    
    plt.scatter(df['Dim_1'].values, df['Dim_2'].values, c = df['knn_color'], alpha = 0.7, cmap = 'Paired')
    plt.savefig (outputfile)
"""

#data_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y']
data_labels = {}

data_labels['DB11994'] = 'A'
data_labels['DB00593'] = 'B'
data_labels['DB00613'] = 'C'
data_labels['DB00172'] = 'D'
data_labels['DB00704'] = 'E'
data_labels['DB00204'] = 'F'
data_labels['DB01104'] = 'G'
data_labels['DB01114'] = 'H'
data_labels['DB00860'] = 'I'
data_labels['DB00770'] = 'J'
data_labels['DB08924'] = 'K'
data_labels['DB06282'] = 'L'
data_labels['DB06736'] = 'M'
data_labels['DB11184'] = 'N'
data_labels['DB00496'] = 'O'
data_labels['DB08877'] = 'P'
data_labels['DB00947'] = 'Q'
data_labels['DB08943'] = 'R'
data_labels['DB01599'] = 'S'
data_labels['DB02546'] = 'T'
data_labels['DB01138'] = 'U'
data_labels['DB00237'] = 'V'
data_labels['DB01325'] = 'W'
data_labels['DB01124'] = 'X'
data_labels['DB00177'] = 'Y'


def plot_multi (df, outputfile, title):

    plt.figure()
 
    # Comment this out if you want image titles.
    #plt.title (title)


    df_query = df[df['not_nn_color'] != 'gray']
    df_not_nn = df[df['not_nn_color'] == 'gray']



    plt.scatter(df_query['Dim_1'].values, df_query['Dim_2'].values, c = df_query['not_nn_color'], alpha = 0.7, marker = 'x', s = 200)


    plt.scatter(df_not_nn['Dim_1'].values, df_not_nn['Dim_2'].values, c = df_not_nn['not_nn_color'], alpha = 0.7, marker = '+', s = 100)

    for i, txt in enumerate(list(df_not_nn['not_nn_label'])):
        plt.annotate('  ' + txt, ((df_not_nn['Dim_1'].values[i], df_not_nn['Dim_2'].values[i])))



    plt.savefig (outputfile, dpi=300)




   

def get_letter (id, labels):
    return(labels[id])

   

def compute_sim (df, str1_col, str2_col, id1_col, id2_col, radius, fplength):

    df['not_nn_mol_query'] = df.apply (lambda x: get_mol (x[str1_col]), axis = 1)
    df['not_nn_mol_target'] = df.apply (lambda x: get_mol (x[str2_col]), axis = 1)
    df['not_nn_fp_query'] = df.apply (lambda x: get_fingerprint (x['not_nn_mol_query'], radius, fplength), axis = 1)
    df['not_nn_fp_target'] = df.apply (lambda x: get_fingerprint (x['not_nn_mol_target'], radius, fplength), axis = 1)
    df['sim'] = df.apply (lambda x: get_Tanimoto(x['not_nn_fp_query'], x['not_nn_fp_target']), axis = 1)
    
    return (df)


def get_query_color (id, colors):
    return (colors[id])

def make_plots (df, query_ids, not_nn_ids, outfile, title):
    df_rnd = df[df['ID'].isin(query_ids)].copy()
    df_not_nn = df[df['ID'].isin(not_nn_ids)].copy()

    # df_rnd: randomly selected 5  query molecules (same that were randomly selected for KNN analysis)
    df_rnd['not_nn_color'] = df_rnd.apply (lambda x: get_query_color (x['ID'], colors), axis = 1)

    # df_not_nn: randomly selected 25 molecules, that are not neightbors of the query molecules, at least no similarity information was used to select them.
    df_not_nn['not_nn_color'] = 'gray'

    df_not_nn = df_not_nn.reset_index (drop = True)
    df_not_nn['idx'] = df_not_nn.index
    df_rnd['idx'] = -1


    df_not_nn['not_nn_label'] = df_not_nn.apply(lambda x: get_letter(x['ID'], data_labels), axis = 1)
    df_rnd['not_nn_label'] = ''

    df_not_nn_final = df_rnd.append(df_not_nn, ignore_index = True)

    #print (df_not_nn_final)

    plot_multi (df_not_nn_final, outfile, title)



# It's sufficient to perform this step once, from any of the datasets created at the space-embedding step, as the structural information is constant, and the structures with nonsense/empty 
# Bemis-Murcko scaffodls were already eliminated.

df = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')

df_valid = knn.filter_valid_mols (df, 'Structure')


df_rnd = df_valid.sample (5, random_state = 55555)

# Remove randomly selected 5 molecules from original (valid) set
df_valid = df_valid[~df_valid['ID'].isin(df_rnd['ID'])]

df_not_nn = df_valid.sample(n = 25, random_state = 11111)

df_rnd['key'] = 1
df_not_nn['key'] = 1

df_paired = df_rnd.merge(df_not_nn, on = 'key', how = 'inner')

df_sim = compute_sim (df_paired, 'Structure_x', 'Structure_y', 'ID_x', 'ID_y', 3, 2048)

#print (df_sim)
#df_all_knns = knn.morgan_knn (df_valid, df_valid, str1_col = 'Structure', str2_col = 'Structure', id1_col = 'ID', id2_col = 'ID', radius = 3, fplength = 2048, k = 5)

df_sim.to_csv ('../../data/similarity_app_drugs_drugbank_all_not_nn_5.tab', sep = '\t', index = False)

colors = {}

colors['DB00006'] = 'green'
colors['DB00849'] = 'orange'
colors['DB00977'] = 'purple'
colors['DB01362'] = 'aqua'
colors['DB04837'] = 'blue'


df_rnd['not_nn_color'] = df_rnd.apply (lambda x: get_query_color (x['ID'], colors), axis = 1)
df_not_nn['not_nn_color'] = 'gray'

df_not_nn = df_not_nn.reset_index (drop = True)
df_not_nn['idx'] = df_not_nn.index
df_rnd['idx'] = -1


df_not_nn['not_nn_label'] = df_not_nn.apply(lambda x: get_letter(x['ID'], data_labels), axis = 1)
df_rnd['not_nn_label'] = ''

df_not_nn_final = df_rnd.append(df_not_nn, ignore_index = True)
df_not_nn_final.to_csv ('../../data/app_drugs_drugbank_all_not_nn_5.tab', sep = '\t', index = False)

query_ids = df_rnd['ID'].values
not_nn_ids = df_not_nn['ID'].values



df_2 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2.tab', sep = '\t')
df_6 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2.tab', sep = '\t')
df_7 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2.tab', sep = '\t')
df_8 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')








make_plots (df_2, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_2_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 2, ChEMBL 24_1 BMSs')
make_plots (df_3, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_3_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 3, ChEMBL 24_1 BMSs')
make_plots (df_4, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_4_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 4, ChEMBL 24_1 BMSs')
make_plots (df_5, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_5_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 5, ChEMBL 24_1 BMSs')
make_plots (df_6, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_6_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 6, ChEMBL 24_1 BMSs')
make_plots (df_7, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_7_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 7, ChEMBL 24_1 BMSs')
make_plots (df_8, query_ids, not_nn_ids, '../../plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_8_dim_2.png', 'Not NN, Approved Drugs DrugBank\nPHC order 8, ChEMBL 24_1 BMSs')








