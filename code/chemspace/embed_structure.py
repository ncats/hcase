# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
#
# Ref: https://pypi.org/project/hilbertcurve/
# Ref: https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.append.html
# Ref: https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
# Ref: https://www.dataquest.io/blog/settingwithcopywarning/
# Ref: https://maxpowerwastaken.github.io/blog/pandas_view_vs_copy/
# Ref: https://engineering.hexacta.com/pandas-by-example-columns-547696ff78dd
# Ref: https://realpython.com/python-rounding/
# Ref: https://github.com/matplotlib/matplotlib/issues/3466/
#
#
#
# New trick: to create NxM matrix with dataframes:
# 1. Create a column of an arbitrary but constant value in both dataframes.
# 2. Inner-join these two columns on that column.
# That's it!
# GZK

import sys
import pandas as pd
import rdkit
from rdkit import Chem
from hilbertcurve.hilbertcurve import HilbertCurve
from scaffold_keys import smiles2bmscaffold, smiles2scaffoldkey, sk_distance

#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
#import seaborn as sns

import math

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
    df = df.sort_values (['sk_distance', 'sk_struct'])
    closest_scaffold_order = df['order'].values[0]
    
    return (closest_scaffold_order)
    
def get_bucket_id (closest_order, bucket_size):
    bucket_id = int(round (closest_order / bucket_size)) + 1
    print ('Closest order: %d, bucket id: %d' % (closest_order, bucket_id))

    return (bucket_id)


def get_hilbert_coordinates (hc, bucket_id):
#    print (bucket_id)
    coordinates = []
    coordinates = hc.coordinates_from_distance(bucket_id - 1)
    coordinate_str = ''
    nr_dim = len (coordinates)

    for i in range(nr_dim):
        coordinate_str += str(coordinates[i]) + ';'

    coordinate_str = coordinate_str[:-1]

    return (coordinate_str)


def embed_structures_hcase (df_structures, str_colname, id_colname, df_space, hc_order, n_dim):
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
    

print ('[SYNTAX] python embed_structure.py <tab-separated file of structures> <Column name in structures file containing SMILES> <Column name of structure IDs> <Hilbert-Curve Mapped Scaffold Embeding> <Order of Hilbert-Curve> <Number of Dimensions> <Name stem of output data frame and plot>')

fname_structures = sys.argv[1]
str_colname = sys.argv[2]
id_colname = sys.argv[3]
fname_space = sys.argv[4]
hc_order = int(sys.argv[5])
n_dim = int(sys.argv[6])
outfile_stem = sys.argv[7]

df_structures = pd.read_csv (fname_structures, sep = '\t')
df_space = pd.read_csv (fname_space, sep = '\t')

print (df_space.head())
print (df_space.shape)




df = embed_structures_hcase (df_structures, str_colname, id_colname, df_space, hc_order, n_dim)

#output_df_fname = '../../data/' + outfile_stem + '.tab'
output_df_fname = outfile_stem + '.tab'

df = df.drop (columns = ['index'])

df.to_csv (output_df_fname, index = False, sep = '\t')


"""
#output_plot_fname = '../../plots/' + outfile_stem + '.png'
output_plot_fname = outfile_stem + '.png'



plt.scatter(df['Dim_1'].values, df['Dim_2'].values)
plt.savefig (output_plot_fname)
"""
print ('[Done.]')

