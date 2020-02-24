# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
#
# Ref: https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jaccard.html#scipy.spatial.distance.jaccard
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/24426410/
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

import sys
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from scaffold_keys import smiles2bmscaffold, smiles2scaffoldkey, sk_distance

#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
#import seaborn as sns

import numpy as np
from sklearn.manifold import TSNE

import math

from knn import get_mol, get_fingerprint, is_valid

def read_data(fname):
    df = pd.read_csv(fname, sep ='\t')

    return (df)


def tr_expand_coords (df, source_col, id_col, delimiter):
    df_orig = df
    df = df[source_col].str.split(delimiter, expand = True)
    nr_cols = len (df.columns)
    columns = []
    for i in range (nr_cols):
        columns.append('Dim_' + str(i + 1))
    
    df.columns = columns
    df = df.astype('int32')
    #df[id_col] = df_orig[id_col]
    
    df = pd.concat([df_orig, df], axis = 1)

    return (df)

def get_coordinates (hc, bucket_id):
#    print (bucket_id)
    coordinates = []
    coordinates = hc.coordinates_from_distance(bucket_id - 1)
    coordinate_str = ''
    nr_dim = len (coordinates)

    for i in range(nr_dim):
        coordinate_str += str(coordinates[i]) + ';'

    coordinate_str = coordinate_str[:-1]

    return (coordinate_str)


def to_bitstring (fp):
    fpstr = 'NA'
    try:
        fpstr = fp.ToBitString()
    except:
        fpstr = 'NA'
    
    return (fpstr)

def fp_gen_with_errohandling (mol):
    fp = None
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = 3, nBits = 2048)
    except:
        fp = None
    
    return (fp)
 

def generate_fp (df):
    df['is_valid'] = df.apply (lambda x: is_valid(x['Structure']), axis = 1)
    df = df[df['is_valid'] == True].copy()
    df['mol'] = df.apply(lambda x: get_mol (x['Structure']), axis = 1)
    print (df.dtypes)
    print(df.columns)
    print (df.head)
    df['fp'] = df.apply (lambda x: fp_gen_with_errohandling(x['mol']), axis = 1)
    print (df.columns)
    df = df[df['fp'] != None].copy()
    df['fp_str'] = df.apply (lambda x: to_bitstring(x['fp']), axis = 1)
    df = df[df['fp_str'] != 'NA'].copy()


    return (df)
    

def get_fp_np_array (fps):
    first = True
    fpl = 0
    all_fp = []
    fp_array = []
    for i in range(len(fps)):
        fp_array = []
        fp = fps[i]
        if first:
            fpl = len(fp)
            first = False
        else:
            if len(fp) != fpl:
               print ('[ERROR] Fingerprint length mismatch. Terminating ...')
               sys.exit (-1)
        
        for j in range(len(fp)):
            fp_array.append(int(fp[j]))
        
        all_fp.append(fp_array)

    all_fp = np.array(all_fp)

    return (all_fp)


def embed_compounds (df, tsne_model):
    df = generate_fp (df)

    print (df.head)

    X = list(df['fp_str'])
    X = get_fp_np_array (X)

    #print (X)

    X_embedded = tsne_model.fit_transform(X)
    print (X_embedded)

    ids = list(df['ID'])

    df_embedded = pd.DataFrame ({'ID': ids, 'Dim_1': X_embedded[:,0], 'Dim_2': X_embedded[:,1]})
    df_embedded = df_embedded.merge (df, on = 'ID', how = 'inner')

    return (df_embedded)



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
    
    
def plot_multi (df_knn, df_embedded, data_outputfile, outputfile, title):
    #df = df [df[df['hc_order'] == hc_order]
    #df_other = df[df['knn_type'] == 'other']
    embedding_fname = data_outputfile.split ('/')[3]
    embedding_fname = '../../data/' + 'tsne_embedding_coords_all_cmpd_' + embedding_fname  

    df_embedded.to_csv (embedding_fname, sep = '\t', index = False)

    
    df = df_knn.merge (df_embedded, on = 'ID', how = 'inner')
    df.to_csv (data_outputfile, sep = '\t', index = False)
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



def do_analysis (df_drugs, df_knn, perplexity_val):
    tsne_model = TSNE(metric = 'jaccard', random_state = 55555, perplexity = perplexity_val, learning_rate = 200.0, n_iter = 1000)


    df_embedded = embed_compounds (df_drugs, tsne_model)

    print (df_embedded.head)
    print (df_knn.head)
    print (df_embedded.head)

    plot_fname = '../../plots/tsne/orig_tsne_perplexity_' + str(perplexity_val) + '.png'
    data_fname = '../../data/orig_tsne_perplexity_' + str(perplexity_val) + '.tab'


    title_str = 'Originak tSNE KNN (k=5)\n Perplexity = ' +  str(perplexity_val)
    plot_multi (df_knn, df_embedded, data_fname, plot_fname, title_str)

    #log_fname = '../../log/tsne_log_perplexity_' + str(perplexity_val) + '.txt'

    #fplog = open (log_fname, 'w+')
    #fplog.write('[Done.]\n')
    #fplog.close()

def do_left_out_analysis (df_drugs, df_knn, perplexity_val):
    tsne_model = TSNE(metric = 'jaccard', random_state = 55555, perplexity = perplexity_val, learning_rate = 200.0, n_iter = 1000)

    df_drugs_orig = df_drugs
    nr = df_drugs.shape[0]

    nr = math.ceil(nr * 0.9)
    df_drugs = df_drugs.sample (n = nr, random_state = 55555)

    df_check = df_knn[~df_knn['ID'].isin(list(df_drugs['ID']))]
    df_add = df_drugs_orig[df_drugs_orig['ID'].isin(list(df_check['ID']))]
    
    df_drugs = df_drugs.append(df_add, ignore_index = True)

    # Extra check
    df_check = df_knn[~df_knn['ID'].isin(list(df_drugs['ID']))]
 
    if not df_check.empty:
        print (df_check)
        print ('[ERROR] Some of the selected KNN molecules were not included in the set where 10 percent was left out.')
        sys.exit(-2)

    df_embedded = embed_compounds (df_drugs, tsne_model)

    print (df_embedded.head)
    print (df_knn.head)
    print (df_embedded.head)

    plot_fname = '../../plots/tsne/lo_orig_tsne_perplexity_' + str(perplexity_val) + '.png'
    data_fname = '../../data/lo_orig_tsne_perplexity_' + str(perplexity_val) + '.tab'


    title_str = 'Left-out Original tSNE KNN (k=5)\n Perplexity = ' +  str(perplexity_val)
    plot_multi (df_knn, df_embedded, data_fname, plot_fname, title_str)



########################

print ('[SYNTAX] analysis_tsne.py')


perplexity_values = [5, 10, 20, 30, 40, 50]


drugs = '../../data/STD_drugbank_approved_structures_v5.txt'
fname_knn = '../../data/rnd_5_app_drugs_drugbank_knn_5.tab'


#output_data_fname = sys.argv[2]
#output_plot_fname = sys.argv[3]



df_drugs = pd.read_csv (drugs, sep = '\t')
df_drugs = df_drugs[['ID', 'Structure']].copy()
#df_drugs = df_drugs.sample(n = 5, random_state = 55555)
#print (df_drugs)


df_knn = read_data (fname_knn) # knn_query_structure , knn_target_structure, knn_query_id, knn_target_id, knn_color
df_knn = separate_query_from_target_mols (df_knn)

for i in range(len(perplexity_values)):
    perplexity_val = perplexity_values[i]
    do_analysis (df_drugs, df_knn, perplexity_val)
    do_left_out_analysis (df_drugs, df_knn, perplexity_val)


print ('[Done.]')

