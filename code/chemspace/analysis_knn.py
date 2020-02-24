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
import pandas as pd


def separate_query_from_target_mols (df_knn):
    df_target = df_knn[['knn_target_id', 'knn_target_structure', 'knn_color', 'data_label']].copy()
    df_target['knn_type'] = 'nn'
    df_target = df_target.rename(columns = {
        'knn_target_id': 'ID',
        'knn_target_structure': 'smiles'
    #    'knn_color': 'color'
    })

    df_query = df_knn[['knn_query_id', 'knn_query_structure', 'knn_color', 'data_label']].copy()
    df_query['knn_type'] = 'parent'
    df_query = df_query.rename(columns = {
        'knn_query_id': 'ID',
        'knn_query_structure': 'smiles'
    #    'knn_color': 'color'
    })

    df_query = df_query.groupby(['ID'], as_index = False).agg({
        'smiles': 'first',
         'knn_color': 'first',
         'knn_type': 'first',
         'data_label': 'first'
    })

    df_query = df_query.reset_index(drop = True)

    df_knn = df_query.append(df_target, ignore_index = True)

    return (df_knn)


#data_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y']

data_labels = {}

data_labels['DB04931'] = 'A'
data_labels['DB01284'] = 'B'
data_labels['DB00050'] = 'C'
data_labels['DB09067'] = 'D'
data_labels['DB06825'] = 'E'
data_labels['DB01174'] = 'F'
data_labels['DB00794'] = 'G'
data_labels['DB05246'] = 'H'
data_labels['DB01437'] = 'I'
data_labels['DB00252'] = 'J'
data_labels['DB01357'] = 'K'
data_labels['DB04575'] = 'L'
data_labels['DB00655'] = 'M'
data_labels['DB00783'] = 'N'
data_labels['DB04573'] = 'O'
data_labels['DB01249'] = 'P'
data_labels['DB09135'] = 'Q'
data_labels['DB09134'] = 'R'
data_labels['DB09313'] = 'S'
data_labels['DB01578'] = 'T'
data_labels['DB11609'] = 'U'
data_labels['DB00257'] = 'V'
data_labels['DB00333'] = 'W'
data_labels['DB01231'] = 'X'
data_labels['DB08944'] = 'Y'

def plot_multi (df, df_knn, outputfile, title):
    #df = df [df[df['hc_order'] == hc_order]
    df_other = df[df['knn_type'] == 'other']
    
    
    df_1 = df[df['knn_color'] == 1]
    df_1_parent = df_1[df_1['knn_type'] == 'parent']
    df_1_nn = df_1[df_1['knn_type'] == 'nn']

    #print (df.head)
    #print (df_knn.head)
    df_1_nn = df_1_nn.drop(columns = ['knn_color'])
    df_1_nn = df_1_nn.merge(df_knn, on = 'ID', how = 'inner')
    
    #print (df_1_parent.head)
    #print (df_1_nn.head)


    df_2 = df[df['knn_color'] == 2]
    df_2_parent = df_2[df_2['knn_type'] == 'parent']
    df_2_nn = df_2[df_2['knn_type'] == 'nn']

    df_2_nn = df_2_nn.drop(columns = ['knn_color'])
    df_2_nn = df_2_nn.merge(df_knn, on = 'ID', how = 'inner')



    df_3 = df[df['knn_color'] == 3]
    df_3_parent = df_3[df_3['knn_type'] == 'parent']
    df_3_nn = df_3[df_3['knn_type'] == 'nn']

    df_3_nn = df_3_nn.drop(columns = ['knn_color'])
    df_3_nn = df_3_nn.merge(df_knn, on = 'ID', how = 'inner')


    df_4 = df[df['knn_color'] == 4]
    df_4_parent = df_4[df_4['knn_type'] == 'parent']
    df_4_nn = df_4[df_4['knn_type'] == 'nn']

    df_4_nn = df_4_nn.drop(columns = ['knn_color'])
    df_4_nn = df_4_nn.merge(df_knn, on = 'ID', how = 'inner')


    df_5 = df[df['knn_color'] == 5]
    df_5_parent = df_5[df_5['knn_type'] == 'parent']
    df_5_nn = df_5[df_5['knn_type'] == 'nn']

    df_5_nn = df_5_nn.drop(columns = ['knn_color'])
    df_5_nn = df_5_nn.merge(df_knn, on = 'ID', how = 'inner')


    plt.figure()
 
    # Comment this out if you want image titles.
    #plt.title (title)

    plt.scatter(df_other['Dim_1'].values, df_other['Dim_2'].values, c = 'gray', alpha = 0.3, marker = '.')

    plt.scatter(df_1_parent['Dim_1'].values, df_1_parent['Dim_2'].values, c = 'blue', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_1_nn['Dim_1'].values, df_1_nn['Dim_2'].values, c = 'blue', alpha = 0.7, marker = '+', s = 100)

    for i, txt in enumerate(list(df_1_nn['data_label'])):
        plt.annotate('  ' + txt, ((df_1_nn['Dim_1'].values[i], df_1_nn['Dim_2'].values[i])))

    plt.scatter(df_2_parent['Dim_1'].values, df_2_parent['Dim_2'].values, c = 'green', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_2_nn['Dim_1'].values, df_2_nn['Dim_2'].values, c = 'green', alpha = 0.7, marker = '+', s = 100)

    for i, txt in enumerate(list(df_2_nn['data_label'])):
        plt.annotate('  ' + txt, ((df_2_nn['Dim_1'].values[i], df_2_nn['Dim_2'].values[i])))



    plt.scatter(df_3_parent['Dim_1'].values, df_3_parent['Dim_2'].values, c = 'purple', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_3_nn['Dim_1'].values, df_3_nn['Dim_2'].values, c = 'purple', alpha = 0.7, marker = '+', s = 100)

    for i, txt in enumerate(list(df_3_nn['data_label'])):
        plt.annotate('  ' + txt, ((df_3_nn['Dim_1'].values[i], df_3_nn['Dim_2'].values[i])))



    plt.scatter(df_4_parent['Dim_1'].values, df_4_parent['Dim_2'].values, c = 'orangered', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_4_nn['Dim_1'].values, df_4_nn['Dim_2'].values, c = 'orangered', alpha = 0.7, marker = '+', s = 100)


    for i, txt in enumerate(list(df_4_nn['data_label'])):
        plt.annotate('  ' + txt, ((df_4_nn['Dim_1'].values[i], df_4_nn['Dim_2'].values[i])))



    plt.scatter(df_5_parent['Dim_1'].values, df_5_parent['Dim_2'].values, c = 'aqua', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_5_nn['Dim_1'].values, df_5_nn['Dim_2'].values, c = 'aqua', alpha = 0.7, marker = '+', s = 100)

    for i, txt in enumerate(list(df_5_nn['data_label'])):
        plt.annotate('  ' + txt, ((df_5_nn['Dim_1'].values[i], df_5_nn['Dim_2'].values[i])))



    plt.savefig (outputfile, dpi=300)





def get_color (id, df_knn):
    # This only works if Test 1 and Test 2 below are OK.
    c = 0
    if id in list(df_knn['ID']):
        x = df_knn[df_knn['ID'] == id]
        c = list(x['knn_color'])[0]
    else:
        c = 0
   
    return (c)
    
 
def get_knn_type (id, df_knn):
    t = 'other'
    parents = df_knn[df_knn['knn_type'] == 'parent']
    parent_ids = list(parents['ID'])
    nns = df_knn[df_knn['knn_type'] == 'nn']
    nn_ids = list(nns['ID'])


    if id in parent_ids:
        t = 'parent'
    elif id in nn_ids:
        t = 'nn'
    else:
        t = 'other'
   
    return (t) 

    
def color_knns (df, df_knn):
    df['knn_color'] = df.apply (lambda x: get_color (x['ID'], df_knn), axis = 1)
    df['knn_type'] = df.apply (lambda x: get_knn_type (x['ID'], df_knn), axis = 1)



    return (df)


def check_rnd_knn (df):
    check = True
    #df = pd.read_csv ('../../data/rnd_5_app_drugs_drugbank_knn_5.tab', sep = '\t')

    query_ids = {}
    target_ids = {}


    for id in list(df['knn_query_id']):
        query_ids[id] = True

    query_ids = list(query_ids.keys())

    if len(query_ids) != 5:
        check = False

    for id in query_ids:
        if id in list(df['knn_target_id']):
            check = False


    for id in list(df['knn_target_id']):
        if id not in list(target_ids.keys()):
            target_ids[id] = True
        else:
            check = False

    return (check)

"""
def save_mol_depiction (smiles, id, dataset_name, path):
    DrawingOptions.atomLabelFontSize = 55
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 3.0
    mol = Chem.MolFromSmiles(smiles)
    fname = path + dataset_name + '_' + id + '.png'
    Chem.Draw.MolToFile(mol, fname, size=(1000, 1000))
"""

def get_letter (id):
    return(data_labels[id])

def assign_data_labels (df):
    df['data_label'] = df.apply(lambda x: get_letter(x['knn_target_id']), axis = 1)
    
    return (df)
    

# It's sufficient to perform this step once, from any of the datasets created at the space-embedding step, as the structural information is constant, and the structures with nonsense/empty 
# Bemis-Murcko scaffodls were already eliminated.

df = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')

df_valid = knn.filter_valid_mols (df, 'Structure')





df_all_knns = knn.morgan_knn (df_valid, df_valid, str1_col = 'Structure', str2_col = 'Structure', id1_col = 'ID', id2_col = 'ID', radius = 3, fplength = 2048, k = 5)
df_all_knns.to_csv ('../../data/app_drugs_drugbank_all_knn_5.tab', sep = '\t', index = False)


df_rnd = df_valid.sample (5, random_state = 55555)

df_knns = knn.morgan_knn (df_rnd, df_valid, str1_col = 'Structure', str2_col = 'Structure', id1_col = 'ID', id2_col = 'ID', radius = 3, fplength = 2048, k = 5)

print (df_knns.head)

df_knns = assign_data_labels(df_knns)
df_knns.to_csv ('../../data/rnd_5_app_drugs_drugbank_knn_5.tab', sep = '\t', index = False)

# Test 1: Making sure that using the 55555 seed and K = 5 results in KNNs with no overlap (i.e. a each molecule in the union of all KNN sets is only member of one KNN set.
# Test 2: Making sure that the randomly selected molecules are not part of any KNN set (to exclude non-sense coloring)):
res = check_rnd_knn (df_knns)
if not res:
    print ('[ERROR]: Randomly selected molecules for KNN analysis is not an ideal set. Check that no query molecule appears as target molecule and that no target molecule is replicated in the target set.')
    print ('Terminating...')
    sys.exit(-1)
else:
    print ('[CHECK] Result: OK. Randomly selected KNN set is all good.')


df_2 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2.tab', sep = '\t')
df_6 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2.tab', sep = '\t')
df_7 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2.tab', sep = '\t')
df_8 = pd.read_csv ('../../data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')


df_knns = separate_query_from_target_mols (df_knns)

df_2_color = color_knns(df_2, df_knns)
df_3_color = color_knns(df_3, df_knns)
df_4_color = color_knns(df_4, df_knns)
df_5_color = color_knns(df_5, df_knns)
df_6_color = color_knns(df_6, df_knns)
df_7_color = color_knns(df_7, df_knns)
df_8_color = color_knns(df_8, df_knns)

plot_multi (df_2_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 2, ChEMBL 24_1 BMSs')
plot_multi (df_3_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 3, ChEMBL 24_1 BMSs')
plot_multi (df_4_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 4, ChEMBL 24_1 BMSs')
plot_multi (df_5_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 5, ChEMBL 24_1 BMSs')
plot_multi (df_6_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 6, ChEMBL 24_1 BMSs')
plot_multi (df_7_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 7, ChEMBL 24_1 BMSs')
plot_multi (df_8_color, df_knns, '../../plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.png', 'KNN, k=5, Approved Drugs DrugBank\nPHC order 8, ChEMBL 24_1 BMSs')


#print (df_knns[['knn_query_id', 'knn_target_id', 'knn_sim', 'knn_color']])



# Save KNN query and target molecule images
depiction_path = '../../plots/knn/molecules/'
dataset_name = 'app_drugbank'

df_2_color['phc_order'] = 2
df_3_color['phc_order'] = 3
df_4_color['phc_order'] = 4
df_5_color['phc_order'] = 5
df_6_color['phc_order'] = 6
df_7_color['phc_order'] = 7
df_8_color['phc_order'] = 8


df_coords = df_2_color.append(df_3_color, ignore_index = True)
df_coords = df_coords.append(df_4_color, ignore_index = True)
df_coords = df_coords.append(df_5_color, ignore_index = True)
df_coords = df_coords.append(df_6_color, ignore_index = True)
df_coords = df_coords.append(df_7_color, ignore_index = True)
df_coords = df_coords.append(df_8_color, ignore_index = True)

df_coords.to_csv ('../../data/knn_coords_app_drugs_drugbank_chembl_24_1_bms.tab', sep = '\t', index = False)






