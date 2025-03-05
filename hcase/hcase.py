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

import math
import numpy as np
import cupy as cp

from hilbertcurve.hilbertcurve import HilbertCurve
from hcase.fingerprints import smiles2scaffoldkey, smiles2bmscaffold



def order_scaffolds_np(structures: np.ndarray, pattern_ids: np.ndarray) -> tuple:

    scaffold_keys = np.array([smiles2scaffoldkey(str(s)) for s in structures])

    # Continue with the rest of your code...
    valid_mask = np.array(scaffold_keys) != 'NA'
    structures, scaffold_keys, pattern_ids = (
        structures[valid_mask], scaffold_keys[valid_mask], pattern_ids[valid_mask]
    )

    sort_indices = np.argsort(scaffold_keys)
    structures, scaffold_keys, pattern_ids, scaffold_keys_one = (
        structures[sort_indices], scaffold_keys[sort_indices], pattern_ids[sort_indices], scaffold_keys[sort_indices]
    )

    _, unique_indices = np.unique(scaffold_keys_one, return_index=True)
    structures, scaffold_keys, pattern_ids = (
        structures[unique_indices], scaffold_keys[unique_indices], pattern_ids[unique_indices]
    )

    order = np.arange(1, len(structures) + 1)

    return structures, order, pattern_ids, scaffold_keys


def get_bucket_id(closest_order: int, bucket_size: float) -> int:
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



def compute_max_phc_order(ref_scaffold_smiles: np.ndarray) -> int:
    log_base = 4
    M = ref_scaffold_smiles.shape[0]
    max_z = math.ceil(math.log(M, log_base))
    return int(max_z)


def preprocess_scaffolds(ref_scaffold_smiles: np.ndarray) -> np.ndarray:
    """
    Preprocesses scaffold keys (from smiles) for efficient distance computation.
    
    Args:
        ref_scaffold_smiles: Array containing scaffold keys (SMILES or similar representations).
        
    Returns:
        numpy array: Preprocessed scaffold keys for efficient distance computation.
    """
    # Split scaffold keys into individual elements (assuming space-separated strings)
    scaffold_keys_split = np.array([s.split(' ') for s in ref_scaffold_smiles], dtype=object)
    return scaffold_keys_split

def compute_distances_vectorized(sk_struct_arr: np.ndarray, ref_scaffolds: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Vectorized distance computation."""
    diff = np.abs(ref_scaffolds[:, None, :] - sk_struct_arr[None, :, :]) ** 1.5  # Shape (num_scaffolds, num_structures, num_features)
    distances = np.einsum('ijk,k->ij', diff, weights)  # Shape (num_scaffolds, num_structures)
    
    return distances

def batch_process_find_closest_scaffolds(sk_struct_arr, ref_scaffolds, weights, batch_size, use_cupy=False):
    """Batch processing for finding closest scaffolds."""
    num_structures = sk_struct_arr.shape[0]
    all_distances = []  # We'll store all the distance computations here

    # Iterate through the data in batches
    for i in range(0, num_structures, batch_size):
        sk_struct_batch = sk_struct_arr[i:i + batch_size]
        
        if use_cupy:  # GPU-based computation with CuPy
            sk_struct_batch_cp = cp.asarray(sk_struct_batch, dtype=cp.float32)
            distances_batch = compute_distances_vectorized(
                sk_struct_batch_cp, cp.asarray(ref_scaffolds), cp.asarray(weights)
            )
            all_distances.append(distances_batch.get())  # Transfer result back to CPU
        else:  # CPU-based computation with NumPy
            distances_batch = compute_distances_vectorized(
                sk_struct_batch, ref_scaffolds, weights
            )
            all_distances.append(distances_batch)

    # Concatenate all batches' distances
    all_distances = np.concatenate(all_distances, axis=1)  # Shape (num_scaffolds, num_structures)

    # Find the closest scaffold for each structure after concatenation
    closest_indices = np.argmin(all_distances, axis=0)  # Get closest scaffold indices for each structure

    return closest_indices

def embed(ref_fingerprints: np.ndarray, structures: np.ndarray, ids: np.ndarray, n_dim: int, df_space_order: np.ndarray, use_cupy=False, batch_size=500) -> np.ndarray:
    max_z = compute_max_phc_order(ref_fingerprints)

    # Preprocess scaffold fingerprints
    ref_scaffolds = preprocess_scaffolds(ref_fingerprints)
    ref_scaffolds = ref_scaffolds.astype(float)

    # Precompute weights
    num_features = ref_scaffolds.shape[1]
    weights = 1 / np.arange(1, num_features + 1)

    # Compute `bms` for the scaffolds using `smiles2bmscaffold`
    bms = np.array([smiles2bmscaffold(smiles) for smiles in structures])  # Apply the function to the array of SMILES strings

    # Filter out rows where 'bms' is 'NA'
    valid_indices = bms != 'NA'
    structures = structures[valid_indices]
    ids = ids[valid_indices]
    bms = bms[valid_indices]

    # Compute `sk_struct` using `smiles2scaffoldkey` on `bms`
    sk_struct = np.array([smiles2scaffoldkey(bm) for bm in bms])

    # Filter out rows where 'sk_struct' is 'NA'
    valid_indices = sk_struct != 'NA'
    structures = structures[valid_indices]
    ids = ids[valid_indices]
    sk_struct = sk_struct[valid_indices]
    bms = bms[valid_indices]
    
    # Preprocess sk_struct into NumPy arrays
    sk_struct = np.vstack([np.fromstring(x, sep=' ', dtype=np.float32) for x in sk_struct])

    # Prepare output list for results
    results = []

    # Use batching to process the scaffolds
    # Precompute the closest indices for scaffolds
    closest_indices = batch_process_find_closest_scaffolds(sk_struct, ref_scaffolds, weights, batch_size, use_cupy)

    # Now, compute closest_order using NumPy indexing:
    # closest_order is obtained by indexing into df_space_order using closest_indices
    closest_order = df_space_order[closest_indices]  # Directly index into df_space_order (NumPy array)

    # Iterate over the desired Hilbert curve order
    nr_scaffolds = ref_fingerprints.shape[0]
    for hc_order in range(2, max_z + 1):
        bucket_nr = math.pow(math.pow(2, hc_order), n_dim)
        bucket_size = float(nr_scaffolds / (bucket_nr - 1))

        # Assign a unique bucket ID based on the closest_order
        bucket_ids = np.round(closest_order / bucket_size) + 1

        # Compute embedding coordinates
        hilbert_curve = HilbertCurve(hc_order, n_dim)
        embedded_hs_coordinates = np.array([get_hilbert_coordinates(hilbert_curve, b) for b in bucket_ids])
        # Expand coordinates into individual dimensions
        expanded_coords = np.char.split(embedded_hs_coordinates.astype(str), ';')

        # Convert list of strings to a 2D numpy array of integers
        expanded_coords = np.array([np.array(coords, dtype=int) for coords in expanded_coords])
        # Combine all necessary information for the result
        results.append(np.column_stack([
            ids,                         # ID column
            structures,                  # Structure column
            sk_struct,                   # sk_struct column
            closest_order,               # closest_order column
            bucket_ids,                  # bucket_id column
            expanded_coords              # Dim_1, Dim_2, etc.
        ]))

    # Concatenate all results for different Hilbert curve orders
    df_res = np.vstack(results)

    return df_res



