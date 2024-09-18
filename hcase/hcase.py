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

try:
    from rdkit import Chem
except ImportError:
    raise ImportError(
        "RDKit is required for this functionality. Please install it via conda: 'conda install -c conda-forge rdkit'")

import math
import pandas as pd
from typing import Callable, Any, Optional
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

from rdkit import RDLogger

from hilbertcurve.hilbertcurve import HilbertCurve
from hcase.scaffold_keys import smiles2bmscaffold, smiles2scaffoldkey, sk_distance, onestring

from logging import getLogger


lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

logger = getLogger()

# def process_chunk(chunk: pd.DataFrame, func: Callable[[Any], Any], kwargs: dict, column: Optional[str] = None) -> pd.Series:
#     """
#     Processes a chunk of rows or columns by applying the given function.
    
#     Args:
#         chunk: A chunk of the DataFrame to process.
#         func: The function to apply to each row or column element.
#         kwargs: Additional arguments to pass to the function.
#         column: If provided, applies the function to a specific column of the chunk; otherwise applies to the whole row.
    
#     Returns:
#         A Series with the results of applying the function to each row or column element.
#     """
#     if column:
#         # Apply function to a specific column in the chunk
#         return pd.Series([func(x, **kwargs) for x in chunk[column]]).reset_index(drop=True)
#     else:
#         # Apply function to entire rows
#         return pd.Series([func(row, **kwargs) for _, row in chunk.iterrows()]).reset_index(drop=True)

# def apply_function(
#     df: pd.DataFrame, 
#     func: Callable[[Any], Any], 
#     n_cores: int = 1, 
#     column: Optional[str] = None,  # Allow specifying a column to apply the function on
#     **kwargs
# ) -> pd.DataFrame:
#     """
#     Generalized function to apply `func` either to an entire row or a specific column of a DataFrame.
    
#     Args:
#         df: DataFrame on which to apply the function.
#         func: Function to apply on each row or element of the DataFrame.
#         n_cores: Number of cores to use for parallel processing. If 1, runs sequentially.
#         column: The specific column to apply the function on (optional).
#         kwargs: Additional arguments to pass to `func`.
    
#     Returns:
#         A DataFrame or Series with the applied function.
#     """
#     # If a specific column is provided, apply function to that column
#     if column:
#         if n_cores == 1:
#             # Sequential processing: apply function to the specific column
#             return df[column].map(lambda x: func(x, **kwargs)).reset_index(drop=True)
#         else:
#             data_split = np.array_split(df, n_cores)  # Split the full DataFrame
#             pool = Pool(n_cores)

#             # Process each chunk (with a specific column) and reset index
#             result = pd.concat(
#                 tqdm(pool.starmap(process_chunk, [(chunk, func, kwargs, column) for chunk in data_split]), total=n_cores)
#             ).reset_index(drop=True)

#             pool.close()
#             pool.join()
#             return result
#     else:
#         # If no column is specified, apply function to entire rows
#         if n_cores == 1:
#             return df.apply(lambda row: func(row, **kwargs), axis=1)
#         else:
#             data_split = np.array_split(df, n_cores)
#             pool = Pool(n_cores)

#             # Process each chunk (rows) and reset index
#             result = pd.concat(
#                 tqdm(pool.starmap(process_chunk, [(chunk, func, kwargs, None) for chunk in data_split]), total=n_cores)
#             ).reset_index(drop=True)

#             pool.close()
#             pool.join()
#             return result

# def apply_sk(df: pd.DataFrame, n_cores: int = 1) -> pd.DataFrame:
#     return apply_function(df, sk, n_cores, column='structure')

# def apply_sk_one(df: pd.DataFrame, n_cores: int = 1) -> pd.DataFrame:
#     return apply_function(df, sk_one, n_cores, column='sk')

# # def apply_smiles2bmscaffold(df: pd.DataFrame, n_cores: int = 1) -> pd.DataFrame:
# #     return apply_function(df, smiles2bmscaffold, n_cores, column='bms')

# # def apply_smiles2scaffoldkey(df: pd.DataFrame, trailing_inchikey: bool = False, n_cores: int = 1) -> pd.DataFrame:
# #     return apply_function(df, smiles2scaffoldkey, n_cores=n_cores, column='structure', trailing_inchikey=trailing_inchikey)

# # def apply_closest_scaffold(df: pd.DataFrame, df_space: pd.DataFrame, nr_structures: int, n_cores: int = 1) -> pd.DataFrame:
# #     return apply_function(df, closest_scaffold_func, n_cores=n_cores, df_space=df_space, nr_structures=nr_structures)

# def apply_get_bucket_id(df: pd.DataFrame, bucket_size: int, n_cores: int = 1) -> pd.DataFrame:
#     """
#     Applies the get_bucket_id function to the 'closest_order' column of the DataFrame.

#     Args:
#         df: DataFrame containing the 'closest_order' column.
#         bucket_size: The size of each bucket.
#         n_cores: Number of cores to use for parallel processing. If 1, runs sequentially.

#     Returns:
#         A DataFrame with the bucket IDs assigned.
#     """
#     return apply_function(df, get_bucket_id_func, 'closest_order', n_cores, bucket_size=bucket_size)

# def apply_get_hilbert_coordinates(df: pd.DataFrame, hilbert_curve: HilbertCurve, n_cores: int = 1) -> pd.DataFrame:
#     """
#     Applies the get_hilbert_coordinates function to the 'bucket_id' column of the DataFrame.

#     Args:
#         df: DataFrame containing the 'bucket_id' column.
#         hilbert_curve: The HilbertCurve object to use for mapping.
#         n_cores: Number of cores to use for parallel processing. If 1, runs sequentially.

#     Returns:
#         A DataFrame with the Hilbert space coordinates assigned.
#     """
#     return apply_function(df, get_hilbert_coordinates_func, 'bucket_id', n_cores, hilbert_curve=hilbert_curve)



# def apply_closest_scaffold_on_split(data_split, df_space, nr_structures):
#     return data_split.apply(lambda x: closest_scaffold(x['sk_struct'], df_space, x['idx'], nr_structures), axis=1)

# def parallel_apply_closest_scaffold(df, df_space, nr_structures, func, n_cores=1):
#     data_split = np.array_split(df[['sk_struct', 'idx']], n_cores)
#     pool = Pool(n_cores)
#     data = pd.concat(tqdm(pool.starmap(func, [(split, df_space, nr_structures) for split in data_split]), total=n_cores))
#     pool.close()
#     pool.join()
#     return data

# def apply_smiles2bmscaffold_on_split(data_split):
#     return (data_split.apply(smiles2bmscaffold))


# def parallel_apply_smiles2bmscaffold(df, func, n_cores=1):

        
#     data_split = np.array_split(df['structure'], n_cores)
#     pool = Pool(n_cores)
#     data = pd.concat(tqdm(pool.map(func, data_split), total=n_cores))
#     pool.close()
#     pool.join()
#     return (data)

# def apply_smiles2scaffoldkey_on_split(data_split, trailing_inchikey = False):
#     return (data_split.apply(smiles2scaffoldkey, trailing_inchikey=trailing_inchikey))

# def parallel_apply_smiles2scaffoldkey (df, func, trailing_inchikey = False, n_cores = 4):
#     data_split = np.array_split(df['bms'], n_cores)
#     pool = Pool(n_cores)
#     data = pd.concat(tqdm(pool.starmap(func, [(split, trailing_inchikey) for split in data_split]), total=n_cores))
#     pool.close()
#     pool.join()
#     return (data)

def apply_sk_on_split(data_split):
    return (data_split.apply(sk))



def parallel_apply_sk(df, func, n_cores=1):

        
    data_split = np.array_split(df['structure'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.map(func, data_split), total=n_cores))
    pool.close()
    pool.join()
    return (data)


def apply_sk_one_on_split(data_split):
    return (data_split.apply(sk_one))


def parallel_apply_sk_one(df, func, n_cores=1):

        
    data_split = np.array_split(df['sk'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.map(func, data_split), total=n_cores))
    pool.close()
    pool.join()
    return (data)


def apply_smiles2bmscaffold_on_split(data_split):
    return (data_split.apply(smiles2bmscaffold))


def parallel_apply_smiles2bmscaffold(df, func, n_cores=1):

        
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

def parallel_apply_closest_scaffold(df, df_space, nr_structures, func, n_cores=1):
    data_split = np.array_split(df[['sk_struct', 'idx']], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, df_space, nr_structures) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return data

def apply_get_bucket_id_on_split(data_split, bucket_size):
    return data_split.apply(lambda x: get_bucket_id(x, bucket_size))

def parallel_apply_get_bucket_id(df, bucket_size, func, n_cores=1):
    data_split = np.array_split(df['closest_order'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, bucket_size) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return data


def apply_get_hilbert_coordinates_on_split(data_split, hilbert_curve):
    return data_split.apply(lambda x: get_hilbert_coordinates(hilbert_curve, x))

def parallel_apply_get_hilbert_coordinates(df, hilbert_curve, func, n_cores=1):
    data_split = np.array_split(df['bucket_id'], n_cores)
    pool = Pool(n_cores)
    data = pd.concat(tqdm(pool.starmap(func, [(split, hilbert_curve) for split in data_split]), total=n_cores))
    pool.close()
    pool.join()
    return data


# def closest_scaffold_func(row: pd.Series, df_space: pd.DataFrame, nr_structures: int) -> int:
#     """
#     Finds the closest scaffold for a given structure based on the 'sk_struct' and 'idx' columns of the row.

#     Args:
#         row: A row of the DataFrame containing 'sk_struct' and 'idx'.
#         df_space: DataFrame containing the reference scaffolds.
#         nr_structures: Total number of structures.

#     Returns:
#         The closest scaffold as an integer (order).
#     """
#     return closest_scaffold(row['sk_struct'], df_space, row['idx'], nr_structures)


def order_scaffolds(df: pd.DataFrame, n_cores: int = 1) -> pd.DataFrame:
    """
    Orders scaffolds based on structure and scaffold key.

    Args:
        df: Pandas DataFrame with the following columns:
            - pattern_id: includes the pattern type and a unique ID separated by a "dot"
              example: scaffold.4
            - structure: SMILES structure
            - ptype: type of scaffold, example: scaffold
            - hash: InChIKey, example: NOWKCMXCCJGMRR-UHFFFAOYSA-N
        n_cores: Number of cores to use for parallel processing. If 1, runs sequentially.

    Returns:
        A DataFrame ordered by scaffold keys with structure, order, scaffold_id, and scaffold_key.
    """
    logger.info('[*] Ordering reference scaffolds ..')

    # Filter the dataframe to include only rows where ptype is 'scaffold'
    df = df[df['ptype'] == 'scaffold']
    nr_orig_scaffolds = df.shape[0]

    # Apply sk function (structure key) in parallel or sequentially
    df['sk'] = parallel_apply_sk(df, apply_sk_on_split, n_cores)

    # Filter out rows where the scaffold key is 'NA'
    df = df[df['sk'] != 'NA']

    # Apply sk_one function (structure key variant) in parallel or sequentially
    df['sk_one'] = parallel_apply_sk_one(df, apply_sk_one_on_split, n_cores)

    # Sort by 'sk_one'
    df = df.sort_values(['sk_one'])

    # Group by 'sk_one' and take the first occurrence of 'sk', 'pattern_id', and 'structure'
    df = df.groupby(['sk_one'], as_index=False).agg({
        'sk': 'first',
        'pattern_id': 'first',
        'structure': 'first'
    })

    # Drop the 'sk_one' column and rename columns for clarity
    df = df.drop(columns=['sk_one']).rename(columns={
        'sk': 'scaffold_key',
        'pattern_id': 'scaffold_id'
    })

    # Reset the index and create a new 'order' column
    df = df.reset_index(drop=True)
    df['order'] = df.index + 1

    # Reorganize the DataFrame to include only the necessary columns
    df = df[['structure', 'order', 'scaffold_id', 'scaffold_key']].copy()

    # Print information about the number of scaffolds
    logger.info('[*] Number of scaffolds in input:')
    logger.info(nr_orig_scaffolds)

    logger.info('[*] Number of unique reference scaffolds:')
    logger.info(df.shape[0])

    logger.info('[*] Done.')

    return df


def define_reference_scaffolds(df: pd.DataFrame) -> pd.DataFrame:
    """
    Defines the reference scaffolds by extracting relevant columns.

    Args:
        df: Pandas DataFrame with the following columns:
            - pattern_id: includes the pattern type and a unique ID separated by a "dot"
              example: scaffold.4
            - structure: SMILES structure
            - ptype: type of scaffold, example: scaffold
            - hash: InChIKey, example: NOWKCMXCCJGMRR-UHFFFAOYSA-N

    Returns:
        A new DataFrame containing the 'pattern_id', 'structure', 'ptype', and 'hash' columns.
    """
    # Create a reference DataFrame by selecting the relevant columns
    df_ref = df[['pattern_id', 'structure', 'ptype', 'hash']].copy()

    return df_ref


def sk(smiles: str, trailing_inchikey: bool = True) -> Optional[str]:
    """
    Converts a SMILES string into a scaffold key.

    Args:
        smiles: The SMILES string representing a chemical structure.
        trailing_inchikey: If True, includes trailing InChIKey in the scaffold key.

    Returns:
        The scaffold key as a string, or None if the conversion fails.
    """
    sk = smiles2scaffoldkey(smiles, trailing_inchikey)
    return sk


def sk_one(sk: str, has_inchikey: bool = True) -> Optional[str]:
    """
    Converts a scaffold key into a "one string" representation.

    Args:
        sk: The scaffold key.
        has_inchikey: If True, indicates that the scaffold key includes an InChIKey.

    Returns:
        The "one string" representation of the scaffold key, or None if the conversion fails.
    """
    sk_one_str = onestring(sk, has_inchikey)
    return sk_one_str


def tr_expand_coords(df: pd.DataFrame, source_col: str, id_col: str, delimiter: str) -> pd.DataFrame:
    """
    Expands coordinates from a delimited string column into multiple dimensions.

    Args:
        df: Pandas DataFrame containing the data.
        source_col: The column containing delimited coordinate strings.
        id_col: The column containing unique identifiers (currently not used in the function).
        delimiter: The delimiter used to split the coordinate string.

    Returns:
        A DataFrame with expanded coordinates in separate columns.
    """
    df_orig = df.copy()
    df_expanded = df[source_col].str.split(delimiter, expand=True)

    # Rename expanded columns
    nr_cols = len(df_expanded.columns)
    columns = [f'Dim_{i+1}' for i in range(nr_cols)]
    df_expanded.columns = columns
    df_expanded = df_expanded.astype('int32')

    # Concatenate original dataframe with expanded coordinates
    df_result = pd.concat([df_orig, df_expanded], axis=1)

    return df_result


def closest_scaffold(sk_struct: Any, df_space: pd.DataFrame, idx: int, nr_structures: int) -> int:
    """
    Finds the closest reference scaffold for a given structure.

    Args:
        sk_struct: The scaffold structure.
        df_space: DataFrame containing scaffold keys and related information.
        idx: The index of the structure.
        nr_structures: Total number of structures.

    Returns:
        The order of the closest scaffold as an integer.
    """
    df = df_space.copy()
    df['sk_struct'] = sk_struct

    df['sk_distance'] = df.apply(lambda x: sk_distance(x['sk_struct'], x['scaffold_key']), axis=1)
    df = df.sort_values(['sk_distance', 'sk_struct'])

    closest_scaffold_order = df['order'].values[0]
    return int(closest_scaffold_order)


def get_bucket_id(closest_order: int, bucket_size: int) -> int:
    """
    Calculates the bucket ID based on the closest scaffold order and bucket size.

    Args:
        closest_order: The order of the closest scaffold.
        bucket_size: The size of each bucket.

    Returns:
        The bucket ID as an integer.
    """
    bucket_id = int(round(closest_order / bucket_size)) + 1
    return bucket_id


def get_hilbert_coordinates(hc: HilbertCurve, bucket_id: int) -> str:
    """
    Retrieves Hilbert curve coordinates based on the bucket ID.

    Args:
        hc: The HilbertCurve object.
        bucket_id: The ID of the bucket.

    Returns:
        A string representing the Hilbert coordinates, separated by semicolons.
    """
    coordinates = hc.point_from_distance(bucket_id - 1)
    coordinate_str = ';'.join(str(coord) for coord in coordinates)

    return coordinate_str


def train(df: pd.DataFrame, n_cores: int = 1) -> pd.DataFrame:
    """
    Trains the model by ordering scaffolds and preparing the reference scaffold set.

    Args:
        df: Pandas DataFrame with the following columns:
            - pattern_id: includes the pattern type and a unique ID separated by a "dot"
              example: scaffold.4
            - structure: SMILES structure
            - ptype: type of scaffold, example: scaffold
            - hash: InChIKey, example: NOWKCMXCCJGMRR-UHFFFAOYSA-N
        n_cores: Number of cores to use for parallel processing. If 1, runs sequentially.

    Returns:
        A DataFrame representing the ordered reference scaffolds set (df_space).
    """
    # Define reference scaffolds
    df = define_reference_scaffolds(df)

    # Order the reference scaffolds by their Scaffold Keys (SKs)
    df = order_scaffolds(df, n_cores)

    df_space = df
    return df_space


def compute_max_phc_order(df_space: pd.DataFrame) -> int:
    """
    Computes the maximum PHC order based on the number of reference scaffolds.

    Args:
        df_space: DataFrame containing the reference scaffolds.

    Returns:
        The maximum PHC order (int).
    """
    log_base = 4

    # Number of reference scaffolds
    M = df_space.shape[0]

    # Compute the maximum PHC order using log base 4
    max_z = math.ceil(math.log(M, log_base))

    return int(max_z)


def embed(df_space: pd.DataFrame, df_structures: pd.DataFrame, n_dim: int, n_cores: int = 1) -> pd.DataFrame:
    """
    Embeds structures based on reference scaffolds using Pseudo Hilbert Curves (PHCs).

    Args:
        df_space: Pandas DataFrame generated by the hcase.train() method.
        df_structures: Pandas DataFrame containing the structures with the following columns:
            - structure: SMILES representation of the structure.
            - id: unique identifier of the compound.
        n_dim: Number of dimensions for the Hilbert curve.
        n_cores: Number of cores to use for parallel processing, 1 will be sequential (default: 1).

    Returns:
        A DataFrame containing the embedded Hilbert space coordinates for the structures.
    """
    # Compute the maximum order of the Pseudo Hilbert Curves (PHCs)
    max_z = compute_max_phc_order(df_space)
    str_colname = 'structure'
    id_colname = 'id'

    # Work with a copy of df_structures and limit it to the necessary columns
    df_structures = df_structures[[id_colname, str_colname]].copy()

    logger.info(df_structures.columns)
    logger.info('[*] Generating Bemis-Murcko scaffolds for compounds ..')

    # Generate Bemis-Murcko scaffolds for the compounds in parallel
    df_structures['bms'] = parallel_apply_smiles2bmscaffold(df_structures, apply_smiles2bmscaffold_on_split, n_cores)

    # Filter out invalid scaffolds
    df_structures = df_structures[df_structures['bms'] != 'NA']
    logger.info('[*] .. done')

    # Rename structure column in df_space for clarity
    df_space = df_space.rename(columns={'structure': 'ref_scaffold_smiles'})
    nr_scaffolds = df_space.shape[0]

    logger.info('[*] Generating Scaffold-Keys for the Bemis-Murcko scaffolds of compounds ..')

    # Generate scaffold keys for the Bemis-Murcko scaffolds in parallel
    df_structures['sk_struct'] = parallel_apply_smiles2scaffoldkey (df_structures,
                                                                  apply_smiles2scaffoldkey_on_split,
                                                                  trailing_inchikey = False,
                                                                  n_cores = n_cores)

    df_structures = df_structures[df_structures['sk_struct'] != 'NA']
    df_structures = df_structures.reset_index(drop=True)
    df_structures['idx'] = df_structures.index + 1

    nr_structures = df_structures.shape[0]
    logger.info('[*] Identifying the closest reference scaffolds of compounds ..')

    # Identify the closest reference scaffolds in parallel
    # df_structures['closest_order'] = apply_closest_scaffold(df_structures, df_space, nr_structures, n_cores=n_cores)
    df_structures['closest_order'] = parallel_apply_closest_scaffold(df_structures, df_space, nr_structures, apply_closest_scaffold_on_split, n_cores)
 
    logger.info('[*] .. done')

    # Prepare for embedding by iterating over the parameter z (PHC order)
    df_res = pd.DataFrame()
    first = True

    for hc_order in range(2, max_z + 1):
        bucket_nr = math.pow(math.pow(2, hc_order), n_dim)
        bucket_size = int(round(nr_scaffolds / (bucket_nr - 1)))

        logger.info(f'Generating HCASE embedding at PHC order: {hc_order}, '
                    f'nr_scaffolds: {nr_scaffolds}, bucket_nr: {int(bucket_nr)}, bucket_size {bucket_size:.4f} ..')

        # Initialize Hilbert curve
        hilbert_curve = HilbertCurve(hc_order, n_dim)

        df_hilbert = df_structures.copy()

        logger.info(f'[*] Mapping compounds to pseudo-Hilbert-Curve of z={hc_order} ..')

        # Assign bucket IDs based on closest scaffold order in parallel
        df_hilbert['bucket_id'] = parallel_apply_get_bucket_id(df_hilbert, bucket_size, apply_get_bucket_id_on_split, n_cores=n_cores)
        logger.info('[*] .. done')

        logger.info('[*] Determining the 2D/3D coordinates of compounds in the HCASE map ..')

        # Determine Hilbert space coordinates in parallel
        df_hilbert['embedded_hs_coordinates'] = parallel_apply_get_hilbert_coordinates(df_hilbert, hilbert_curve, apply_get_hilbert_coordinates_on_split, n_cores=n_cores)
        logger.info('[*] .. done')

        # Expand Hilbert space coordinates into separate columns
        df_hilbert = tr_expand_coords(df_hilbert, 'embedded_hs_coordinates', id_colname, delimiter=';')

        df_hilbert['hc_order'] = hc_order

        if first:
            df_res = df_hilbert
            first = False
        else:
            df_res = pd.concat([df_res, df_hilbert], ignore_index=True)

        logger.info('.. done.')

    return df_res
