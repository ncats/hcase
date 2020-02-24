# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov

# Reference scaffold set: ChEMBL

# 1. Import DrugBank KNN set with SKs.
# 2. Import Refernce scaffold set with SK.
# 3. Import tSNE embedding of all reference scaffolds.
# 4. Determine closest reference scaffold based on SKs, compounds assume coordinates of closest reference BMS.
# 5. Generate plots for compounds at each parameter configuration.

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import scaffold_keys as sk
import rdkit
from rdkit import Chem

def read_data(fname):
    df = pd.read_csv(fname, sep ='\t')

    return (df)


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

def generate_bms (df):
    df['bms'] = df.apply (lambda x: sk.smiles2bmscaffold (x['smiles']), axis = 1)
    
    return (df)

def generate_scaffold_keys (df):
    df['sk'] = df.apply (lambda x: sk.smiles2scaffoldkey (x['bms'], trailing_inchikey=True), axis = 1)
    
    return (df)

def closest_scaffold (sk_struct, df_space):
    df = df_space
    df['sk_struct'] = sk_struct
    df['sk_distance'] = df.apply (lambda x: sk.sk_distance (x['sk_struct'], x['scaffold_key']), axis = 1)
    df = df.sort_values (['sk_distance'])
    closest_scaffold_order = df['scaffold_id'].values[0]
    
    return (closest_scaffold_order)




def plot_multi (df_knn, df_embedded, outputfile, title):
    #df = df [df[df['hc_order'] == hc_order]
    #df_other = df[df['knn_type'] == 'other']
    
    df = df_knn.merge (df_embedded, left_on = 'closest_scaffold_id', right_on = 'scaffold_id', how = 'inner')
    print (df)
 
    df_1 = df[df['knn_color'] == 1]
    df_1_parent = df_1[df_1['knn_type'] == 'parent']
    df_1_nn = df_1[df_1['knn_type'] == 'nn']

    #print (df.head)
    #print (df_knn.head)
    #df_1_nn = df_1_nn.drop(columns = ['knn_color'])
    #df_1_nn = df_1_nn.merge(df_knn, left_on = id_col, right_on = 'knn_target_id', how = 'inner')



    df_2 = df[df['knn_color'] == 2]
    df_2_parent = df_2[df_2['knn_type'] == 'parent']
    df_2_nn = df_2[df_2['knn_type'] == 'nn']

    #df_2_nn = df_2_nn.drop(columns = ['knn_color'])
    #df_2_nn = df_2_nn.merge(df_knn, left_on = id_col, right_on = 'knn_target_id', how = 'inner')



    df_3 = df[df['knn_color'] == 3]
    df_3_parent = df_3[df_3['knn_type'] == 'parent']
    df_3_nn = df_3[df_3['knn_type'] == 'nn']

    #df_3_nn = df_3_nn.drop(columns = ['knn_color'])
    #df_3_nn = df_3_nn.merge(df_knn, left_on = id_col, right_on = 'knn_target_id', how = 'inner')


    df_4 = df[df['knn_color'] == 4]
    df_4_parent = df_4[df_4['knn_type'] == 'parent']
    df_4_nn = df_4[df_4['knn_type'] == 'nn']

    #df_4_nn = df_4_nn.drop(columns = ['knn_color'])
    #df_4_nn = df_4_nn.merge(df_knn, left_on = id_col, right_on = 'knn_target_id', how = 'inner')


    df_5 = df[df['knn_color'] == 5]
    df_5_parent = df_5[df_5['knn_type'] == 'parent']
    df_5_nn = df_5[df_5['knn_type'] == 'nn']

    #df_5_nn = df_5_nn.drop(columns = ['knn_color'])
    #df_5_nn = df_5_nn.merge(df_knn, left_on = id_col, right_on = 'knn_target_id', how = 'inner')


    plt.figure()

    # Comment this out if you want image titles.
    #plt.title (title)

    #plt.scatter(df_other['Dim_1'].values, df_other['Dim_2'].values, c = 'gray', alpha = 0.3, marker = '.')

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


def color_knns (df_knn, df_embedded):
    df = df_knn.merge (df_embedded, on = 'scaffold_id', how = 'inner')
    
    return (df)

fname_knn = '../../data/rnd_5_app_drugs_drugbank_knn_5.tab'

fname_ref_bms = '../../data/hc_space.tab'

tsne_5_fname = '../../data/tsne_per_10.tab'
tsne_10_fname = '../../data/tsne_per_20.tab'
tsne_20_fname = '../../data/tsne_per_30.tab'
tsne_30_fname = '../../data/tsne_per_40.tab'
tsne_40_fname = '../../data/tsne_per_50.tab'
tsne_50_fname = '../../data/tsne_per_5.tab'



df_knn = read_data (fname_knn) # knn_query_structure , knn_target_structure, knn_query_id, knn_target_id, knn_color

df_knn = separate_query_from_target_mols (df_knn)

# Structure parsing errors are ignored throughout the script as all molecules must have been processed correctly in previous steps, so that they can appear in the input files.
df_knn = generate_bms(df_knn)
df_knn = generate_scaffold_keys (df_knn)
df_ref_bms = read_data(fname_ref_bms) # structure       order   scaffold_id     scaffold_key
#print(df_ref_bms.head)

df_knn['closest_scaffold_id'] = df_knn.apply(lambda x: closest_scaffold(x['sk'], df_ref_bms), axis = 1)
print (df_knn)



df_5 = read_data(tsne_5_fname)
df_10 = read_data(tsne_10_fname)
df_20 = read_data(tsne_20_fname)
df_30 = read_data(tsne_30_fname)
df_40 = read_data(tsne_40_fname)
df_50 = read_data(tsne_50_fname)

plot_multi (df_knn, df_5, '../../plots/tsne/tsne_knn_perplexity_5.png', 'tSNE KNN (k=5), Perplexity = 5')
plot_multi (df_knn, df_10, '../../plots/tsne/tsne_knn_perplexity_10.png', 'tSNE KNN (k=5), Perplexity = 10')
plot_multi (df_knn, df_20, '../../plots/tsne/tsne_knn_perplexity_20.png', 'tSNE KNN (k=5), Perplexity = 20')
plot_multi (df_knn, df_30, '../../plots/tsne/tsne_knn_perplexity_30.png', 'tSNE KNN (k=5), Perplexity = 30')
plot_multi (df_knn, df_40, '../../plots/tsne/tsne_knn_perplexity_40.png', 'tSNE KNN (k=5), Perplexity = 40')
plot_multi (df_knn, df_50, '../../plots/tsne/tsne_knn_perplexity_50.png', 'tSNE KNN (k=5), Perplexity = 50')







