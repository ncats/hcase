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
from scaffold_keys import smiles2bmscaffold, smiles2scaffoldkey, sk_distance

#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
#import seaborn as sns

import numpy as np
from sklearn.manifold import TSNE

import math

from knn import get_mol, get_fingerprint

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

def closest_scaffold (sk_struct, df_space, idx, nr_structures):
    print ('[*] Processing structure %d out of %d .' % (idx, nr_structures))
    df = df_space
    df['sk_struct'] = sk_struct
    df['sk_distance'] = df.apply (lambda x: sk_distance (x['sk_struct'], x['scaffold_key']), axis = 1)
    df = df.sort_values (['sk_distance'])
    closest_scaffold_order = df['order'].values[0]
    
    return (closest_scaffold_order)
    

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


def altered_tsne (df_structures, str_colname, id_colname, df_space, hc_order, n_dim):
    hilbert_curve = HilbertCurve(hc_order, n_dim)


    df_structures = df_structures[[id_colname, str_colname]].copy()
    df_structures['bms'] = df_structures.apply (lambda x: smiles2bmscaffold (x[str_colname]), axis = 1)

    # filter out invalid of nonsense/empty scaffolds:
    df_structures = df_structures[df_structures['bms'] != 'NA']

    #df_space['jk'] = 1

    df_space = df_space.rename (columns = {
        'structure': 'ref_scaffold_smiles'
    })

    nr_scaffolds = df_space.shape[0]
    bucket_nr = math.pow(math.pow(2, hc_order), n_dim)
    bucket_size = float(nr_scaffolds / (bucket_nr - 1))
    
    
    #df_structures['jk'] = 1
    df_structures['sk_struct'] = df_structures.apply (lambda x: smiles2scaffoldkey(x['bms'], trailing_inchikey = False), axis = 1)
    print ('nr_scaffolds: %d, bucket_nr: %d,  bucket_size %f' % (nr_scaffolds, bucket_nr, bucket_size))


    print ('[*] Number of input structures: %d' % (df_structures.shape[0]))
    df_structures = df_structures[df_structures['sk_struct'] != 'NA']
    print ('[*] Number of structures for which scaffold_key was generated: %d' % (df_structures.shape[0]))
    #df3 = df_space.merge (df_structures, on = 'jk', how = 'inner')

    #df = df3.copy()
    #df = df.reset_index (drop = True)

    df_structures = df_structures.reset_index()
    df_structures['idx'] = df_structures.index + 1


    nr_structures = df_structures.shape[0]
    df_structures['closest_order'] = df_structures.apply (lambda x: closest_scaffold (x['sk_struct'], df_space, x['idx'], nr_structures), axis = 1)
    df_structures['bucket_id'] = df_structures.apply (lambda x: get_bucket_id(x['closest_order'], bucket_size), axis = 1)  
    df_structures['embedded_hs_coordinates'] = df_structures.apply (lambda x: get_hilbert_coordinates(hilbert_curve, x['bucket_id']), axis = 1)  

    df = df_structures
 
    #df = df_structures.merge (df_space, left_on = 'closest_order', right_on = 'order', how = 'inner')
    

    # ignore ->
    """
    df = df.sort_values(['sk_distance'])
    df = df.groupby([str_colname], as_index = False).agg ({
        id_colname: 'first',
        'hs_coordinates': 'first',
        'scaffold_id': 'first',
        'ref_scaffold_smiles': 'first'
    })
    """
    # <- ignore
    
    df = tr_expand_coords (df, 'embedded_hs_coordinates', id_colname, delimiter = ';')
    
    return (df)   


def to_bitstring (fp):
    return (fp.ToBitString())

def generate_fp (df):
    df['bms_mol'] = df.apply(lambda x: get_mol (x['structure']), axis = 1)
    df['bms_fp'] = df.apply (lambda x: get_fingerprint(x['bms_mol'], radius = 3, fplength = 2048), axis = 1)
    df['bms_fp_str'] = df.apply (lambda x: to_bitstring(x['bms_fp']), axis = 1)

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


def embed_ref_scaffolds (df_ref_bms, tsne_model):
    df_ref_bms = generate_fp (df_ref_bms)

    print (df_ref_bms.head)

    X = list(df_ref_bms['bms_fp_str'])
    X = get_fp_np_array (X)

    #print (X)

    X_embedded = tsne_model.fit_transform(X)
    print (X_embedded)

    ids = list(df_ref_bms['scaffold_id'])

    df_embedded = pd.DataFrame ({'scaffold_id': ids, 'Dim_1': X_embedded[:,0], 'Dim_2': X_embedded[:,1]})
    df_embedded = df_embedded.merge (df_ref_bms, on = 'scaffold_id', how = 'inner')

    return (df_embedded)

print ('[SYNTAX] python altered_tsne_analysis.py <perplexity> <output data filename> <output plot filename>')


ref_scaffolds = '../../data/hc_space.tab'
cherry_picked_scaffolds = '../../data/cherrypicked_scaffolds.tab'
drugs = '../../data/STD_drugbank_approved_structures_v5.txt'

perplexity_val = float(sys.argv[1])
output_data_fname = sys.argv[2]
output_plot_fname = sys.argv[3]


# Workflow:

# 1. Generate RDKit Morgan FPs for Scaffolds, radius = 3, length = 2048
# 2. Generate t-SNE for reference scaffolds (ChEBML):
#	a) at default parameters, except using Jaccard distance metric
#       b) at iteration steps = 5000, epsilon = 10, perplexity = [5, 10, 20, 30, 40, 50], Jaccard distance
# 3. Generate plots for Cherry-picked scaffolds on each t-SNE embedding.

df_ref_bms = pd.read_csv (ref_scaffolds, sep = '\t')
df_cp = pd.read_csv (cherry_picked_scaffolds, sep = '\t')
df_drugs = pd.read_csv (drugs, sep = '\t')


# Comment this out for full scale experiment - START
#df_ref_bms = df_ref_bms.sample (n=1000, random_state = 55555)
# Comment this out for full scale experiment - END

#tsne_model = TSNE(metric = 'jaccard', random_state = 55555)

# TSNE(n_components=2, perplexity=30.0, early_exaggeration=12.0, learning_rate=200.0, n_iter=1000, n_iter_without_progress=300, min_grad_norm=1e-07, metric='euclidean', init='random', verbose=0, random_state=None, method='barnes_hut', angle=0.5, n_jobs=None


tsne_model = TSNE(metric = 'jaccard', random_state = 55555, perplexity = perplexity_val, learning_rate = 200.0, n_iter = 1000)


df_embedded = embed_ref_scaffolds (df_ref_bms, tsne_model)


df_embedded_cp = df_embedded[df_embedded['scaffold_id'].isin(df_cp['scaffold_id'])]
df_embedded_not_cp = df_embedded[~df_embedded['scaffold_id'].isin(df_cp['scaffold_id'])]

df_embedded_cp = df_embedded_cp.merge(df_cp, on = 'scaffold_id', how = 'inner')
df_embedded_cp = df_embedded_cp[['scaffold_id', 'Dim_1', 'Dim_2', 'color']].copy()


df_embedded_not_cp = df_embedded_not_cp[['scaffold_id', 'Dim_1', 'Dim_2']].copy()
df_embedded_not_cp['color'] = -1

df_embedded = df_embedded_cp.append (df_embedded_not_cp, ignore_index = True)

df_embedded.to_csv (output_data_fname, sep ='\t', index = False)

print (df_embedded.head)


tab10_palette = ['#1179B0', '#F58C30', '#74BB5A', '#BC412C', '#795B9A', '#764A0C', '#D37DB5', '#7A7A7A', '#B8C449']



df_1 = df_embedded[df_embedded['color'] == 1]
df_2 = df_embedded[df_embedded['color'] == 2]
df_3 = df_embedded[df_embedded['color'] == 3]
df_4 = df_embedded[df_embedded['color'] == 4]
df_5 = df_embedded[df_embedded['color'] == 5]
df_6 = df_embedded[df_embedded['color'] == 6]
df_7 = df_embedded[df_embedded['color'] == 7]
df_8 = df_embedded[df_embedded['color'] == 8]
df_9 = df_embedded[df_embedded['color'] == 9]

title = 't-SNE Analysis of ChEMBL Scaffolds at Perplexity = ' + str(perplexity_val)

plt.figure()

#plt.title (title)

if not df_1.empty:
    color_1 = list(df_1['color'])[0] - 1
    plt.scatter(df_1['Dim_1'].values, df_1['Dim_2'].values, c = tab10_palette[color_1], alpha = 0.3, marker = 'o', s = 100)

if not df_2.empty:
    color_2 = list(df_2['color'])[0] - 1
    plt.scatter(df_2['Dim_1'].values, df_2['Dim_2'].values, c = tab10_palette[color_2], alpha = 0.3, marker = 'o', s = 100)

if not df_3.empty:
    color_3 = list(df_3['color'])[0] - 1
    plt.scatter(df_3['Dim_1'].values, df_3['Dim_2'].values, c = tab10_palette[color_3], alpha = 0.3, marker = 'o', s = 100)

if not df_4.empty:
    color_4 = list(df_4['color'])[0] - 1
    plt.scatter(df_4['Dim_1'].values, df_4['Dim_2'].values, c = tab10_palette[color_4], alpha = 0.3, marker = 'o', s = 100)

if not df_5.empty:
    color_5 = list(df_5['color'])[0] - 1
    plt.scatter(df_5['Dim_1'].values, df_5['Dim_2'].values, c = tab10_palette[color_5], alpha = 0.3, marker = 'o', s = 100)

if not df_6.empty:
    color_6 = list(df_6['color'])[0] - 1
    plt.scatter(df_6['Dim_1'].values, df_6['Dim_2'].values, c = tab10_palette[color_6], alpha = 0.3, marker = 'o', s = 100)

if not df_7.empty:
    color_7 = list(df_7['color'])[0] - 1
    plt.scatter(df_7['Dim_1'].values, df_7['Dim_2'].values, c = tab10_palette[color_7], alpha = 0.3, marker = 'o', s = 100)

if not df_8.empty:
    color_8 = list(df_8['color'])[0] - 1
    plt.scatter(df_8['Dim_1'].values, df_8['Dim_2'].values, c = tab10_palette[color_8], alpha = 0.3, marker = 'o', s = 100)

if not df_9.empty:
    color_9 = list(df_9['color'])[0] - 1
    plt.scatter(df_9['Dim_1'].values, df_9['Dim_2'].values, c = tab10_palette[color_9], alpha = 0.3, marker = 'o', s = 100)



plt.savefig (output_plot_fname, dpi=300)

log_fname = '../../log/tsne_log_perplexity_' + str(perplexity_val) + '.txt'

fplog = open (log_fname, 'w+')
fplog.write('[Done.]\n')
fplog.close()

print ('[Done.]')

