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

#import knn
import pandas as pd


def plot_multi (df, outputfile, title):
    print(df.columns)

    #print(df.head)
    

    # Comment this out if you want image titles.
    #plt.title (title)


    tab10_palette = ['#1179B0', '#F58C30', '#74BB5A', '#BC412C', '#795B9A', '#764A0C', '#D37DB5', '#7A7A7A', '#B8C449']

    df_1 = df[df['color'] == 1]
    df_2 = df[df['color'] == 2]
    df_3 = df[df['color'] == 3]
    df_4 = df[df['color'] == 4]
    df_5 = df[df['color'] == 5]
    df_6 = df[df['color'] == 6]
    df_7 = df[df['color'] == 7]
    df_8 = df[df['color'] == 8]
    df_9 = df[df['color'] == 9]


    color_1 = list(df_1['color'])[0] - 1
    color_2 = list(df_2['color'])[0] - 1
    color_3 = list(df_3['color'])[0] - 1
    color_4 = list(df_4['color'])[0] - 1
    color_5 = list(df_5['color'])[0] - 1
    color_6 = list(df_6['color'])[0] - 1
    color_7 = list(df_7['color'])[0] - 1
    color_8 = list(df_8['color'])[0] - 1
    color_9 = list(df_9['color'])[0] - 1
 
    #print("color")
    #print(color_1)




    plt.figure()
 

    plt.scatter(df_1['Dim_1'].values, df_1['Dim_2'].values, c = tab10_palette[color_1], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_2['Dim_1'].values, df_2['Dim_2'].values, c = tab10_palette[color_2], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_3['Dim_1'].values, df_3['Dim_2'].values, c = tab10_palette[color_3], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_4['Dim_1'].values, df_4['Dim_2'].values, c = tab10_palette[color_4], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_5['Dim_1'].values, df_5['Dim_2'].values, c = tab10_palette[color_5], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_6['Dim_1'].values, df_6['Dim_2'].values, c = tab10_palette[color_6], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_7['Dim_1'].values, df_7['Dim_2'].values, c = tab10_palette[color_7], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_8['Dim_1'].values, df_8['Dim_2'].values, c = tab10_palette[color_8], alpha = 0.3, marker = 'o', s = 100)
    plt.scatter(df_9['Dim_1'].values, df_9['Dim_2'].values, c = tab10_palette[color_9], alpha = 0.3, marker = 'o', s = 100)









    plt.savefig (outputfile, dpi=300)


"""
#    for i, txt in enumerate(list(df_1_nn['data_label'])):
#        plt.annotate('  ' + txt, ((df_1_nn['Dim_1'].values[i], df_1_nn['Dim_2'].values[i])))

    plt.scatter(df_2_parent['Dim_1'].values, df_2_parent['Dim_2'].values, c = 'green', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_2_nn['Dim_1'].values, df_2_nn['Dim_2'].values, c = 'green', alpha = 0.7, marker = '+', s = 100)

#    for i, txt in enumerate(list(df_2_nn['data_label'])):
#        plt.annotate('  ' + txt, ((df_2_nn['Dim_1'].values[i], df_2_nn['Dim_2'].values[i])))



    plt.scatter(df_3_parent['Dim_1'].values, df_3_parent['Dim_2'].values, c = 'purple', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_3_nn['Dim_1'].values, df_3_nn['Dim_2'].values, c = 'purple', alpha = 0.7, marker = '+', s = 100)

#    for i, txt in enumerate(list(df_3_nn['data_label'])):
#        plt.annotate('  ' + txt, ((df_3_nn['Dim_1'].values[i], df_3_nn['Dim_2'].values[i])))



    plt.scatter(df_4_parent['Dim_1'].values, df_4_parent['Dim_2'].values, c = 'orangered', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_4_nn['Dim_1'].values, df_4_nn['Dim_2'].values, c = 'orangered', alpha = 0.7, marker = '+', s = 100)


#    for i, txt in enumerate(list(df_4_nn['data_label'])):
#        plt.annotate('  ' + txt, ((df_4_nn['Dim_1'].values[i], df_4_nn['Dim_2'].values[i])))



    plt.scatter(df_5_parent['Dim_1'].values, df_5_parent['Dim_2'].values, c = 'aqua', alpha = 0.7, marker = 'x', s = 200)
    plt.scatter(df_5_nn['Dim_1'].values, df_5_nn['Dim_2'].values, c = 'aqua', alpha = 0.7, marker = '+', s = 100)

#    for i, txt in enumerate(list(df_5_nn['data_label'])):
#        plt.annotate('  ' + txt, ((df_5_nn['Dim_1'].values[i], df_5_nn['Dim_2'].values[i])))

"""




    
def color_scaffolds (df_coords, df_color):
    df = df_coords.merge(df_color, on = 'scaffold_id', how = 'inner')
    
    return (df)


"""
def save_mol_depiction (smiles, id, dataset_name, path):
    DrawingOptions.atomLabelFontSize = 55
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 3.0
    mol = Chem.MolFromSmiles(smiles)
    fname = path + dataset_name + '_' + id + '.png'
    Chem.Draw.MolToFile(mol, fname, size=(1000, 1000))
"""

def get_letter (idx, labels):
    return(labels[idx])

def assign_data_labels (df, labels):
    df = df.reset_index(drop = True)
    df['label_index'] = df.index
    df['data_label'] = df.apply(lambda x: get_letter(x['label_index'], labels), axis = 1)
    df = df.drop(columns = ['label_index'])
    
    return (df)
    







# Logic:

# Join color of selected scaffolds with embedding data.

# Create plots



# It's sufficient to perform this step once, from any of the datasets created at the space-embedding step, as the structural information is constant, and the structures with nonsense/empty 
# Bemis-Murcko scaffodls were already eliminated.


df_cp = pd.read_csv ('../../data/cherrypicked_scaffolds.tab', sep = '\t')

df_2 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_2_dim_2.tab', sep = '\t')
df_3 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_3_dim_2.tab', sep = '\t')
df_4 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_4_dim_2.tab', sep = '\t')
df_5 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_5_dim_2.tab', sep = '\t')
df_6 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_6_dim_2.tab', sep = '\t')
df_7 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_7_dim_2.tab', sep = '\t')
df_8 = pd.read_csv ('../../data/cp_scaffolds_chembl_24_1_bms_ord_8_dim_2.tab', sep = '\t')


df_2_color = color_scaffolds(df_2, df_cp)
df_3_color = color_scaffolds(df_3, df_cp)
df_4_color = color_scaffolds(df_4, df_cp)
df_5_color = color_scaffolds(df_5, df_cp)
df_6_color = color_scaffolds(df_6, df_cp)
df_7_color = color_scaffolds(df_7, df_cp)
df_8_color = color_scaffolds(df_8, df_cp)

plot_multi (df_2_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_2_dim_2.png', 'Cherry-Picked BMSs\nPHC order 2, ChEMBL 24_1 BMSs')
plot_multi (df_3_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_3_dim_2.png', 'Cherry-Picked BMSs\nPHC order 3, ChEMBL 24_1 BMSs')
plot_multi (df_4_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_4_dim_2.png', 'Cherry-Picked BMSs\nPHC order 4, ChEMBL 24_1 BMSs')
plot_multi (df_5_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_5_dim_2.png', 'Cherry-Picked BMSs\nPHC order 5, ChEMBL 24_1 BMSs')
plot_multi (df_6_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_6_dim_2.png', 'Cherry-Picked BMSs\nPHC order 6, ChEMBL 24_1 BMSs')
plot_multi (df_7_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_7_dim_2.png', 'Cherry-Picked BMSs\nPHC order 7, ChEMBL 24_1 BMSs')
plot_multi (df_8_color, '../../plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_8_dim_2.png', 'Cherry-Picked BMSs\nPHC order 8, ChEMBL 24_1 BMSs')













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

df_coords.to_csv ('../../data/cp_scaffolds_all_embedding_coords_chembl_24_1_bms.tab', sep = '\t', index = False)




