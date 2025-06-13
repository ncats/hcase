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
from tqdm import tqdm
import os
from typing import List, Dict, Tuple

from hilbertcurve.hilbertcurve import HilbertCurve

def _load_index_npz(cache_dir: str, index_filename: str = "index.npz"):
    index_path = os.path.join(cache_dir, index_filename)
    if os.path.exists(index_path):
        data = np.load(index_path, allow_pickle=True)
        identifiers = data["identifiers"]
        batch_files = data["batch_files"]
        row_in_batch = data["row_in_batch"]
        return identifiers, batch_files, row_in_batch
    else:
        return np.array([], dtype=str), np.array([], dtype=str), np.array([], dtype=int)

def _save_index_npz(identifiers, batch_files, row_in_batch, cache_dir: str, index_filename: str = "index.npz"):
    index_path = os.path.join(cache_dir, index_filename)
    np.savez(index_path, identifiers=identifiers, batch_files=batch_files, row_in_batch=row_in_batch)

def _get_next_batch_number_npz(cache_dir: str) -> int:
    existing = [f for f in os.listdir(cache_dir) if f.startswith("batch_") and f.endswith(".npz")]
    if not existing:
        return 1
    nums = [int(f.split("_")[1].split(".")[0]) for f in existing]
    return max(nums) + 1

def _add_closest_indices_batch_npz(identifiers: np.ndarray, closest_indices: np.ndarray, cache_dir: str):
    assert len(identifiers) == closest_indices.shape[0]
    batch_num = _get_next_batch_number_npz(cache_dir)
    batch_file = f"batch_{batch_num:05d}.npz"
    batch_path = os.path.join(cache_dir, batch_file)
    np.savez(batch_path, identifiers=identifiers, closest_indices=closest_indices)
    return batch_file

def _get_closest_indices_from_cache_npz(identifiers: np.ndarray, cache_dir: str, idx_identifiers: np.ndarray, idx_batch_files: np.ndarray, idx_row_in_batch: np.ndarray) -> Dict[str, np.ndarray]:
    # Map: batch_file -> list of (identifier, row_in_batch)
    result = {}
    if len(idx_identifiers) == 0:
        return result
    # Find indices in index for requested identifiers
    id_to_idx = {id_: i for i, id_ in enumerate(idx_identifiers)}
    found = [id_ for id_ in identifiers if id_ in id_to_idx]
    if not found:
        return result
    # Group by batch_file
    batch_map = {}
    for id_ in found:
        idx = id_to_idx[id_]
        batch_file = idx_batch_files[idx]
        row = idx_row_in_batch[idx]
        batch_map.setdefault(batch_file, []).append((id_, row))
    for batch_file, id_row_list in batch_map.items():
        batch_path = os.path.join(cache_dir, batch_file)
        with np.load(batch_path, allow_pickle=True) as data:
            batch_indices = data["closest_indices"]
            for identifier, row_idx in id_row_list:
                result[identifier] = batch_indices[row_idx]
    return result
from hcase.distance_metrics import weighted_minkowski_1_5_distance, euclidean_distance, manhattan_distance, cosine_distance


distance_fn_mapping = {
    'euclidean': euclidean_distance,
    'manhattan': manhattan_distance,
    'cosine': cosine_distance,
    'minkowski_1_5': weighted_minkowski_1_5_distance
}


def compute_max_phc_order(ref_scaffold_smiles: np.ndarray) -> int:
    """
    Compute the maximum PHC order for a given set of reference scaffolds.

    The maximum PHC order is the highest power of 4 that is less than or equal to the number of reference scaffolds.

    The equation is Ceil(log(M, 4)), where M is the number of reference scaffolds.

    Parameters
    ----------
    ref_scaffold_smiles : ndarray
        A 1D numpy array of SMILES strings, where each element is a reference scaffold.

    Returns
    -------
    int
        The maximum PHC order.
    """
    log_base = 4
    M = ref_scaffold_smiles.shape[0]
    max_z = math.ceil(math.log(M, log_base))
    return int(max_z)


def compute_hilbert_embeddings(bucket_ids: np.ndarray, hc_order: int, n_dim: int):
    """
    Compute the Hilbert space embeddings for given bucket IDs.

    This function calculates the coordinates of points in a multidimensional
    Hilbert space using the specified Hilbert Curve order and dimensionality.

    Parameters
    ----------
    bucket_ids : np.ndarray
        A 1D numpy array of integers representing the bucket IDs. Each bucket ID
        corresponds to a point on the Hilbert curve.
    hc_order : int
        The order of the Hilbert Curve, which determines its resolution.
    n_dim : int
        The number of dimensions of the Hilbert space.

    Returns
    -------
    np.ndarray
        A 2D numpy array where each row contains the coordinates of a point in
        the Hilbert space corresponding to a bucket ID.
    """

    hilbert_curve = HilbertCurve(hc_order, n_dim)
    embedded_hs_coordinates = np.array([np.array(hilbert_curve.point_from_distance(b - 1)) for b in bucket_ids])
    return embedded_hs_coordinates


def process_buckets(bucket_size: int, closest_order: np.ndarray, hc_order: int, n_dim: int):

    """
    Processes bucket information and computes Hilbert space embeddings.

    This function assigns bucket IDs to the given closest order values based 
    on the specified bucket size. It then computes the corresponding Hilbert 
    space coordinates for these bucket IDs using the provided Hilbert Curve 
    order and dimensionality.

    Parameters
    ----------
    bucket_size : int
        The size of each bucket used for assigning bucket IDs.
    closest_order : np.ndarray
        A 1D numpy array of order values for which bucket IDs need to be 
        assigned.
    hc_order : int
        The order of the Hilbert Curve, which determines its resolution.
    n_dim : int
        The number of dimensions of the Hilbert space.

    Returns
    -------
    tuple
        A tuple containing:
        - bucket_ids: np.ndarray
            A 1D numpy array of integers representing the assigned bucket IDs.
        - embedded_hs_coordinates: np.ndarray
            A 2D numpy array where each row contains the coordinates of a point 
            in the Hilbert space corresponding to a bucket ID.
    """

    bucket_ids = (np.round(closest_order / bucket_size) + 1).astype(int)  # Assign bucket ID
    embedded_hs_coordinates = compute_hilbert_embeddings(bucket_ids, hc_order, n_dim)

    return bucket_ids, embedded_hs_coordinates


def embed(target_fingerprints: np.ndarray, ref_fingerprints: np.ndarray = None, target_ids: List[str] = None,
          ref_ids: List[str] = None,
          n_dim: int = 2, 
          use_cupy=False, target_batch_size=1000, ref_batch_size=50000, hc_order=None, 
          distance_fn='minkowski_1_5', max_hc_order_only=True, cache_dir: str = "embedding_cache"):

    if distance_fn not in distance_fn_mapping:
        raise ValueError(f"Unknown distance function: {distance_fn}. Choose from {list(distance_fn_mapping.keys())}")
    distance_fn = distance_fn_mapping[distance_fn]

    # 1. Sort target fingerprints & target_ids
    sort_idx_target = np.lexsort(target_fingerprints.T[::-1])
    target_fingerprints = target_fingerprints[sort_idx_target]
    if target_ids is not None:
        target_ids = np.array(target_ids)[sort_idx_target]
    else:
        target_ids = np.array([str(i) for i in range(target_fingerprints.shape[0])])

    # 2. Create ref_fingerprints and unique indices from sorted target_fingerprints
    if ref_fingerprints is None:
        ref_fingerprints, unique_indices = np.unique(target_fingerprints, axis=0, return_index=True)
        if ref_ids is None and target_ids is not None:
            ref_ids = target_ids[unique_indices]

    if use_cupy:
        import cupy as cp
        target_fingerprints = cp.asarray(target_fingerprints, dtype=cp.float32)

    num_targets = target_fingerprints.shape[0]
    num_refs = ref_fingerprints.shape[0]
    if num_targets < target_batch_size:
        target_batch_size = num_targets
    if num_refs < ref_batch_size:
        ref_batch_size = num_refs
    max_z = compute_max_phc_order(ref_fingerprints)

    hc_order_vals = [max_z] if max_hc_order_only else (range(2, max_z + 1) if hc_order is None else [hc_order])

    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        idx_identifiers, idx_batch_files, idx_row_in_batch = _load_index_npz(cache_dir)
    else:
        idx_identifiers = np.array([], dtype=str)
        idx_batch_files = np.array([], dtype=str)
        idx_row_in_batch = np.array([], dtype=int)

    # Determine which target_ids are already cached
    id_to_idx = {id_: i for i, id_ in enumerate(idx_identifiers)}
    cached_mask = np.array([id_ in id_to_idx for id_ in target_ids])
    uncached_mask = ~cached_mask

    num_cached = np.sum(cached_mask)
    num_uncached = np.sum(uncached_mask)
    print(f"[HCASE] {num_cached} targets found in cache, {num_uncached} to compute.")

    # Prepare arrays for all results
    closest_indices_all = np.full(num_targets, -1, dtype=int)

    # Fill in cached results
    if np.any(cached_mask):
        cached_ids = target_ids[cached_mask]
        print(f"[HCASE] Loading cached closest_indices for {len(cached_ids)} targets.")
        cached_closest_indices = _get_closest_indices_from_cache_npz(cached_ids, cache_dir, idx_identifiers, idx_batch_files, idx_row_in_batch)
        for i, id_ in enumerate(target_ids):
            if cached_mask[i] and id_ in cached_closest_indices:
                closest_indices_all[i] = cached_closest_indices[id_]

    # Compute for uncached targets
    uncached_indices = np.where(uncached_mask)[0]
    if len(uncached_indices) > 0:
        print(f"[HCASE] Computing closest_indices for {len(uncached_indices)} uncached targets...")
        for batch_start in tqdm(range(0, len(uncached_indices), target_batch_size), desc="Processing uncached batches", unit="batch"):
            batch_idx = uncached_indices[batch_start:batch_start + target_batch_size]
            target_fingerprint_batch = target_fingerprints[batch_idx]
            batch_size = target_fingerprint_batch.shape[0]
            min_distances = np.full(batch_size, np.inf, dtype=np.float32)
            closest_indices = np.full(batch_size, -1, dtype=int)
            for j in tqdm(range(0, num_refs, ref_batch_size), desc="Processing reference batches"):
                ref_fingerprint_batch = ref_fingerprints[j:j + ref_batch_size]
                if use_cupy:
                    target_fingerprint_batch_cp = cp.asarray(target_fingerprint_batch, dtype=cp.float32)
                    ref_fingerprint_batch_cp = cp.asarray(ref_fingerprint_batch, dtype=cp.float32)
                    distances_batch = distance_fn(target_fingerprint_batch_cp, ref_fingerprint_batch_cp, use_cupy=True).get()
                else:
                    distances_batch = distance_fn(target_fingerprint_batch, ref_fingerprint_batch)
                batch_closest_indices = np.argmin(distances_batch, axis=0) + j
                batch_min_distances = np.min(distances_batch, axis=0)
                update_mask = batch_min_distances < min_distances
                min_distances[update_mask] = batch_min_distances[update_mask]
                closest_indices[update_mask] = batch_closest_indices[update_mask]
            # Save this batch to cache
            batch_target_ids = target_ids[batch_idx]
            batch_file = _add_closest_indices_batch_npz(batch_target_ids, closest_indices, cache_dir)
            # Update index arrays
            new_identifiers = np.concatenate([idx_identifiers, batch_target_ids])
            new_batch_files = np.concatenate([idx_batch_files, np.array([batch_file] * len(batch_target_ids), dtype=str)])
            new_row_in_batch = np.concatenate([idx_row_in_batch, np.arange(len(batch_target_ids), dtype=int)])
            _save_index_npz(new_identifiers, new_batch_files, new_row_in_batch, cache_dir)
            # Update in-memory index for this run
            idx_identifiers = new_identifiers
            idx_batch_files = new_batch_files
            idx_row_in_batch = new_row_in_batch
            # Fill in results
            closest_indices_all[batch_idx] = closest_indices

    order = np.arange(1, num_refs + 1)
    closest_order = order[closest_indices_all]
    all_embeddings = []
    ref_only_embeddings = []
    all_hc_orders = []

    for hc_order_val in hc_order_vals:
        bucket_nr = math.pow(2, hc_order_val) ** n_dim
        bucket_size = num_refs / (bucket_nr - 1)
        bucket_ids, embedded_hs_coordinates = process_buckets(bucket_size, closest_order, hc_order_val, n_dim)
        bucket_ids, ref_embedded_hs_coordinates = process_buckets(bucket_size, order, hc_order_val, n_dim)
        all_embeddings.append(embedded_hs_coordinates)
        ref_only_embeddings.append(ref_embedded_hs_coordinates)
        all_hc_orders.append(hc_order_val)

    return all_embeddings, ref_only_embeddings, all_hc_orders, ref_ids, closest_indices_all
