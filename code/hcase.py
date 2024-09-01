# Author: Gergely Zahoranszky-Kohalmi, PhD
#
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
#
#
#
# References
#
#
# Ref: https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
# Ref: https://engineering.hexacta.com/pandas-by-example-columns-547696ff78dd
# Ref: https://github.com/matplotlib/matplotlib/issues/3466/
# Ref: https://maxpowerwastaken.github.io/blog/pandas_view_vs_copy/
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.append.html
# Ref: https://pubs.acs.org/doi/10.1021/ci5001983
# Ref: https://pypi.org/project/hilbertcurve/
# Ref: https://realpython.com/python-rounding/
# Ref: https://stackoverflow.com/questions/38862293/how-to-add-incremental-numbers-to-a-new-column-using-pandas/38862389
# Ref: https://towardsdatascience.com/dockerizing-jupyter-projects-39aad547484a
# Ref: https://www.dataquest.io/blog/settingwithcopywarning/
# Ref: https://www.geeksforgeeks.org/log-functions-python/
# Ref: https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
# Ref: https://www.w3schools.com/python/ref_math_ceil.asp
# ChatGPT 4.0 Palantir Instance
# ChatGPT 4o www.openai.com


"""
By ChatGPT 4.0 Palantir Instance
"""

import sys
import math

import pandas as pd

import rdkit
from rdkit import Chem

import warnings
from rdkit import RDLogger


from hilbertcurve.hilbertcurve import HilbertCurve

from scaffold_keys import smiles2bmscaffold, smiles2scaffoldkey, sk_distance, onestring

from multiprocessing import Pool
from tqdm import tqdm

import numpy as np

import pdb



lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)



def apply_sk_on_split(data_split):
    return (data_split.apply(sk))



def parallel_apply_sk(df, func, n_cores=4):

        
    data_split = np.array_split(df['structure'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.map(func, data_split), total=n_cores))
    pool.close()
    pool.join()
    return (data)


def apply_sk_one_on_split(data_split):
    return (data_split.apply(sk_one))


def parallel_apply_sk_one(df, func, n_cores=4):

        
    data_split = np.array_split(df['sk'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.map(func, data_split), total=n_cores))
    pool.close()
    pool.join()
    return (data)


def apply_smiles2bmscaffold_on_split(data_split):
    return (data_split.apply(smiles2bmscaffold))


def parallel_apply_smiles2bmscaffold(df, func, n_cores=4):

        
    data_split = np.array_split(df['structure'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.map(func, data_split), total=n_cores))
    pool.close()
    pool.join()
    return (data)



def apply_smiles2scaffoldkey_on_split(data_split, trailing_inchikey = False):
    return (data_split.apply(smiles2scaffoldkey, trailing_inchikey=trailing_inchikey))

def parallel_apply_smiles2scaffoldkey (df, func, trailing_inchikey = False, n_cores = 4):
    data_split = np.array_split(df['bms'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, trailing_inchikey) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return (data)



def apply_closest_scaffold_on_split(data_split, df_space, nr_structures):
    return data_split.apply(lambda x: closest_scaffold(x['sk_struct'], df_space, x['idx'], nr_structures), axis=1)

def parallel_apply_closest_scaffold(df, df_space, nr_structures, func, n_cores=4):
    data_split = np.array_split(df[['sk_struct', 'idx']], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, df_space, nr_structures) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return data

def apply_get_bucket_id_on_split(data_split, bucket_size):
    return data_split.apply(lambda x: get_bucket_id(x, bucket_size))

def parallel_apply_get_bucket_id(df, bucket_size, func, n_cores=4):
    data_split = np.array_split(df['closest_order'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, bucket_size) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return data


def apply_get_hilbert_coordinates_on_split(data_split, hilbert_curve):
    return data_split.apply(lambda x: get_hilbert_coordinates(hilbert_curve, x))

def parallel_apply_get_hilbert_coordinates(df, hilbert_curve, func, n_cores=4):
    data_split = np.array_split(df['bucket_id'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, hilbert_curve) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return data



def order_scaffolds(df, n_cores):
    """
        df: Pandas data frame structure (columns):

            pattern_id    structure    ptype    hash



            pattern_id: includes the pattern type and a unique ID separated by a "dot"
                    example: scaffold.4

            structure: SMILES

            ptype:  type of scaffold
                    example: scaffold

            hash: inchikey
                   example: NOWKCMXCCJGMRR-UHFFFAOYSA-N

    """
    print ('[*] Ordering reference scaffolds ..')

    df = df[df['ptype'] == 'scaffold']
    nr_orig_scaffolds = df.shape[0]

    # df = df.sample (100, random_state = 55555)

    # Sequential
    #df['sk'] = df.apply(lambda x: sk(x['structure'], trailing_inchikey=True), axis=1)

    # Parallel
    df['sk'] = parallel_apply_sk(df, apply_sk_on_split, n_cores)




    
    df = df[df['sk'] != 'NA']


    
    # Sequential
    #df['sk_one'] = df.apply(lambda x: sk_one(x['sk'], has_inchikey=True), axis=1)

    # Parallel
    df['sk_one'] = parallel_apply_sk_one(df, apply_sk_one_on_split, n_cores)


    
    

    
    df = df.sort_values(['sk_one'])

    df = df.groupby(['sk_one'], as_index=False).agg({
        'sk': 'first',
        'pattern_id': 'first',
        'structure': 'first'
    })

    df = df.drop(columns=['sk_one'])

    df = df.rename(columns={
        'sk': 'scaffold_key',
        'pattern_id': 'scaffold_id'
    })

    df = df.reset_index()
    df['order'] = df.index + 1

    df = df[['structure', 'order', 'scaffold_id', 'scaffold_key']].copy()

    print('[*] Number of scaffolds in input:')
    print(nr_orig_scaffolds)

    print('[*] Number of unique reference scaffolds:')
    print(df.shape[0])

    print ('[*] Done.')
    
    return (df)


def define_reference_scaffolds(df):
    """
        df: Pandas data frame structure (columns):

            pattern_id    structure    ptype    hash



            pattern_id: includes the pattern type and a unique ID separated by a "dot"
                    example: scaffold.4

            structure: SMILES

            ptype:  type of scaffold
                    example: scaffold

            hash:   inchikey
                      example: NOWKCMXCCJGMRR-UHFFFAOYSA-N

    """

    df_ref = df[['pattern_id', 'structure', 'ptype', 'hash']].copy()

    return (df_ref)


def sk(smiles, trailing_inchikey=True):
    sk = smiles2scaffoldkey(smiles, trailing_inchikey)

    return (sk)


def sk_one(sk, has_inchikey=True):
    sk = onestring(sk, has_inchikey)

    return (sk)


def tr_expand_coords(df, source_col, id_col, delimiter):
    df_orig = df
    df = df[source_col].str.split(delimiter, expand=True)
    nr_cols = len(df.columns)
    columns = []
    for i in range(nr_cols):
        columns.append('Dim_' + str(i + 1))

    df.columns = columns
    df = df.astype('int32')
    # df[id_col] = df_orig[id_col]

    df = pd.concat([df_orig, df], axis=1)

    return (df)


def closest_scaffold(sk_struct, df_space, idx, nr_structures):
    # print ('[*] Finding closest reference scaffold for structure %d out of %d .' % (idx, nr_structures))
    df = df_space
    df['sk_struct'] = sk_struct
    df['sk_distance'] = df.apply(lambda x: sk_distance(x['sk_struct'], x['scaffold_key']), axis=1)
    df = df.sort_values(['sk_distance', 'sk_struct'])
    closest_scaffold_order = df['order'].values[0]

    return (int(closest_scaffold_order))


def get_bucket_id(closest_order, bucket_size):
    bucket_id = int(round(closest_order / bucket_size)) + 1
    # print ('Closest order: %d, bucket id: %d' % (closest_order, bucket_id))

    return (bucket_id)


def get_hilbert_coordinates(hc, bucket_id):
    #    print (bucket_id)
    coordinates = []
    coordinates = hc.point_from_distance(bucket_id - 1)
    coordinate_str = ''
    nr_dim = len(coordinates)

    for i in range(nr_dim):
        coordinate_str += str(coordinates[i]) + ';'

    coordinate_str = coordinate_str[:-1]

    return (coordinate_str)


def train(df, n_cores):
    """
        df: Pandas data frame structure (columns):

            pattern_id    structure    ptype    hash



            pattern_id: includes the pattern type and a unique ID separated by a "dot"
                    example: scaffold.4

            structure: SMILES

            ptype:  type of scaffold
                    example: scaffold

            hash: inchikey
                   example: NOWKCMXCCJGMRR-UHFFFAOYSA-N

    """

    # order scaffolds

    df = define_reference_scaffolds(df)

    # Extract reference scaffolds set (unique scaffolds) and order them by their Scaffold Keys (SKs)

    

    df = order_scaffolds(df, n_cores)

    

    df_space = df


    
    return (df_space)


def compute_max_phc_order(df_space):
    """
        Eq: ceil [log_4 (M)]  , derived from Eq 1-3 of manuscript.

            M is the number of reference scaffolds.
    """

    log_base = 4

    max_z = -1

    M = df_space.shape[0]

    max_z = math.ceil(math.log(M, log_base))

    return (int(max_z))






    

def embed(df_space, df_structures, n_dim, n_cores = 4):
    """
        df_space: Pandas data frame, generated by hcase.train() method

        df_structures: Pandas data frame, structure (columns):

            structure id

            structure: SMILES

            id: unique identifier of compound
    """

    # Compute max order of Pseudo Hilbert Curves (PHCs)

    max_z = compute_max_phc_order(df_space)
    n_dim = 2
    str_colname = 'structure'
    id_colname = 'id'

    # Invariant part, i.e. independent of the order of PHC (parameter "z")

    df_structures = df_structures[[id_colname, str_colname]].copy()

    print (df_structures.columns)


    print ('[*] Generating Bemis-Murcko scaffolds for compounds ..')

    # Sequential
    #df_structures['bms'] = df_structures.apply(lambda x: smiles2bmscaffold(x[str_colname]), axis=1)


    # Parallel
    df_structures['bms'] = parallel_apply_smiles2bmscaffold(df_structures, apply_smiles2bmscaffold_on_split, n_cores)
    

    # filter out invalid of nonsense/empty scaffolds:
    df_structures = df_structures[df_structures['bms'] != 'NA']



    print ('[*] .. done')

    
    # df_space['jk'] = 1

    df_space = df_space.rename(columns={
        'structure': 'ref_scaffold_smiles'
    })

    nr_scaffolds = df_space.shape[0]

    # df_structures['jk'] = 1


    print ('[*] Generating Scaffold-Keys for the Bemis-Murcko scaffolds of compounds ..')
    
    # Sequential
    #df_structures['sk_struct'] = df_structures.apply
    #    lambda x: smiles2scaffoldkey(x['bms'], trailing_inchikey=False), axis=1)

    # Parallel
    df_structures['sk_struct'] = parallel_apply_smiles2scaffoldkey (df_structures,
                                                                  apply_smiles2scaffoldkey_on_split,
                                                                  trailing_inchikey = False,
                                                                  n_cores = n_cores)

    print ('[*] .. done')




    print('[*] Number of input structures: %d' % (df_structures.shape[0]))
    df_structures = df_structures[df_structures['sk_struct'] != 'NA']
    print('[*] Number of structures for which scaffold_key was generated: %d' % (df_structures.shape[0]))
    # df3 = df_space.merge (df_structures, on = 'jk', how = 'inner')

    # df = df3.copy()
    # df = df.reset_index (drop = True)

    df_structures = df_structures.reset_index()
    df_structures['idx'] = df_structures.index + 1

    nr_structures = df_structures.shape[0]
    
    

    print ('[*] Identifying the closest reference scaffolds of compounds ..')

    # Sequential:
    #df_structures['closest_order'] = df_structures.apply(
    #    lambda x: closest_scaffold(x['sk_struct'], df_space, x['idx'], nr_structures), axis=1)


    # Parallel:
    df_structures['closest_order'] = parallel_apply_closest_scaffold(df_structures, df_space, nr_structures, apply_closest_scaffold_on_split, n_cores)
 

    print ('[*] .. done') 


    
    df_res = pd.DataFrame()
    first = True

    # Iterating over the parameter z (from 2 to z, z included)

    for hc_order in range(2, max_z + 1):

        bucket_nr = math.pow(math.pow(2, hc_order), n_dim)

        bucket_size = float(nr_scaffolds / (bucket_nr - 1))

        print('Generating HCASE embedding at parameter z (PHC order): %d, nr_scaffolds: %d, bucket_nr: %d,  bucket_size %f ..' % (
            hc_order, nr_scaffolds, bucket_nr, bucket_size))

        hilbert_curve = HilbertCurve(hc_order, n_dim)

        df_hilbert = df_structures.copy()


        
        print (f'[*] Mapping compounds to pseudo-Hilbert-Curve of z={hc_order} ..')
        
        # Sequential
        #df_hilbert['bucket_id'] = df_hilbert.apply(lambda x: get_bucket_id(x['closest_order'], bucket_size), axis=1)


        # Parallel
        df_hilbert['bucket_id'] = parallel_apply_get_bucket_id(df_hilbert, bucket_size, apply_get_bucket_id_on_split)
        
        
        print ('[*] .. done')    



        
        print ('[*] Determining the 2D coordinate of compounds in the HCASE map the bucket id of compounds ..')        

        # Sequential
        #df_hilbert['embedded_hs_coordinates'] = df_hilbert.apply(
        #    lambda x: get_hilbert_coordinates(hilbert_curve, x['bucket_id']), axis=1)


        # Parallel
        df_hilbert['embedded_hs_coordinates'] = parallel_apply_get_hilbert_coordinates(df_hilbert, hilbert_curve, apply_get_hilbert_coordinates_on_split)


        print ('[*] .. done')   
        
        
        df = df_hilbert

        # df = df_structures.merge (df_space, left_on = 'closest_order', right_on = 'order', how = 'inner')

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

        df = tr_expand_coords(df, 'embedded_hs_coordinates', id_colname, delimiter=';')

        df = df.drop(columns=['index'])

        df['hc_order'] = hc_order

        if first:

            df_res = df
            first = False

        else:

            # this might be depricated
            df_res = pd.concat([df_res, df], ignore_index=True)

        print('.. done.')

    return (df_res)
