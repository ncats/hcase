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
import pandas as pd

from hilbertcurve.hilbertcurve import HilbertCurve
from hcase.distance_metrics import weighted_minkowski_1_5_distance, euclidean_distance, manhattan_distance, cosine_distance


distance_fn_mapping = {
    'euclidean': euclidean_distance,
    'manhattan': manhattan_distance,
    'cosine': cosine_distance,
    'minkowski_1_5': weighted_minkowski_1_5_distance
}


def compute_max_phc_order(ref_scaffold_smiles: np.ndarray) -> int:
    log_base = 4
    M = ref_scaffold_smiles.shape[0]
    max_z = math.ceil(math.log(M, log_base))
    return int(max_z)


def compute_hilbert_embeddings(bucket_ids: np.ndarray, hc_order: int, n_dim: int):
    hilbert_curve = HilbertCurve(hc_order, n_dim)
    embedded_hs_coordinates = np.array([np.array(hilbert_curve.point_from_distance(b - 1)) for b in bucket_ids])
    return embedded_hs_coordinates


def process_buckets(bucket_size: int, closest_order: np.ndarray, hc_order: int, n_dim: int):

    bucket_ids = (np.round(closest_order / bucket_size) + 1).astype(int)  # Assign bucket ID
    embedded_hs_coordinates = compute_hilbert_embeddings(bucket_ids, hc_order, n_dim)
    
    return bucket_ids, embedded_hs_coordinates



def embed(ref_fingerprints: np.ndarray, n_dim: int = 2, 
          use_cupy=False, target_batch_size=1500, ref_batch_size=50000, hc_order=None, 
          distance_fn='minkowski_1_5', target_fingerprints: np.ndarray = None, 
          max_hc_order_only=True, output_dir="embedding_results"):
    """
    Embeds new structures based on reference fingerprints.
    """
        
    if distance_fn not in distance_fn_mapping:
        raise ValueError(f"Unknown distance function: {distance_fn}. Choose from {list(distance_fn_mapping.keys())}")
    distance_fn = distance_fn_mapping[distance_fn]

    if use_cupy:
        import cupy as cp
        ref_fingerprints = cp.asarray(ref_fingerprints, dtype=cp.float32)

    if target_fingerprints is None:
        target_fingerprints = ref_fingerprints

    num_targets = target_fingerprints.shape[0]
    num_refs = ref_fingerprints.shape[0]
    max_z = compute_max_phc_order(ref_fingerprints)

    hc_order_vals = [max_z] if max_hc_order_only else (range(2, max_z + 1) if hc_order is None else [hc_order])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i in tqdm(range(0, num_targets, target_batch_size), desc="Processing target batches", unit="batch"):
        all_exist = True
        for hc_order_val in hc_order_vals:
            if not os.path.exists(f"{output_dir}/embedding_batch_{i}_{hc_order_val}.csv"):
                all_exist = False
        if all_exist:
            continue 

        target_fingerprint_batch = target_fingerprints[i:i + target_batch_size]

        # Initialize storage for the closest reference index and min distances
        closest_indices = np.full(target_fingerprint_batch.shape[0], -1, dtype=int)
        min_distances = np.full(target_fingerprint_batch.shape[0], np.inf, dtype=np.float32)

        for j in tqdm(range(0, num_refs, ref_batch_size), desc="Processing reference batches"):
            ref_fingerprint_batch = ref_fingerprints[j:j + ref_batch_size]

            if use_cupy:
                target_fingerprint_batch_cp = cp.asarray(target_fingerprint_batch, dtype=cp.float32)
                ref_fingerprint_batch_cp = cp.asarray(ref_fingerprint_batch, dtype=cp.float32)
                distances_batch = distance_fn(target_fingerprint_batch_cp, ref_fingerprint_batch_cp, use_cupy=True).get()
            else:
                distances_batch = distance_fn(target_fingerprint_batch, ref_fingerprint_batch)

            # Find the closest reference for each target in this batch
            batch_closest_indices = np.argmin(distances_batch, axis=0)
            batch_min_distances = np.min(distances_batch, axis=0)

            batch_closest_indices += j  # Adjust indices to global reference index space

            # Update only where we found a closer reference
            update_mask = batch_min_distances < min_distances
            min_distances[update_mask] = batch_min_distances[update_mask]
            closest_indices[update_mask] = batch_closest_indices[update_mask]

        # Convert closest indices to hierarchical order
        order = np.arange(1, num_refs + 1)
        closest_order = order[closest_indices]


        # Process hierarchical ordering
        for hc_order_val in hc_order_vals:
            bucket_nr = math.pow(math.pow(2, hc_order_val), n_dim)
            bucket_size = float(num_refs / (bucket_nr - 1))
            bucket_ids, embedded_hs_coordinates = process_buckets(bucket_size, closest_order, hc_order_val, n_dim)

            # Save batch results to disk
            df = pd.DataFrame({
                'closest_order': closest_order,
                'bucket_ids': bucket_ids,
                'hc_order_val': [hc_order_val] * len(closest_order)
            })

            for dim in range(embedded_hs_coordinates.shape[1]):
                df[f'coord_{dim}'] = embedded_hs_coordinates[:, dim]

            df.to_csv(f"{output_dir}/embedding_batch_{i}_{hc_order_val}.csv", index=False)

   # Combine final results by hc_order
    combined_results_by_hc_order = {}

    for filename in os.listdir(output_dir):
        if filename.endswith(".csv"):
            df = pd.read_csv(os.path.join(output_dir, filename))
            
            # Iterate over each unique hc_order_val in the batch
            for hc_order_val in np.unique(df['hc_order_val'].to_numpy()):
                if hc_order_val not in combined_results_by_hc_order:
                    # Initialize a new entry for this hc_order_val
                    combined_results_by_hc_order[hc_order_val] = {
                        'closest_order': [],
                        'bucket_ids': [],
                        'embedded_hs_coordinates': [],
                        'hc_order_val': []
                    }

                # Append results for this hc_order_val
                mask = df['hc_order_val'] == hc_order_val
                combined_results_by_hc_order[hc_order_val]['closest_order'].append(df[mask]['closest_order'].to_numpy())
                combined_results_by_hc_order[hc_order_val]['bucket_ids'].append(df[mask]['bucket_ids'].to_numpy())

                # hc_order_val is the same across the batch, so we can just repeat it
                hc_order_val_array = np.full(df[mask]['hc_order_val'].shape, hc_order_val)
                combined_results_by_hc_order[hc_order_val]['hc_order_val'].append(hc_order_val_array)

                # Concatenate coordinates (all columns that start with 'coord_')
                embedded_hs_coordinates = df[mask][[col for col in df.columns if col.startswith('coord_')]].to_numpy()
                combined_results_by_hc_order[hc_order_val]['embedded_hs_coordinates'].append(embedded_hs_coordinates)

    # Now, concatenate arrays within each hc_order
    final_combined_results = {
        'closest_order': [],
        'bucket_ids': [],
        'embedded_hs_coordinates': [],
        'hc_order_val': []
    }

    for hc_order_val, result in combined_results_by_hc_order.items():
        final_combined_results['closest_order'].append(np.concatenate(result['closest_order']))
        final_combined_results['bucket_ids'].append(np.concatenate(result['bucket_ids']))
        final_combined_results['embedded_hs_coordinates'].append(np.concatenate(result['embedded_hs_coordinates']))
        final_combined_results['hc_order_val'].append(np.concatenate(result['hc_order_val']))

    return (final_combined_results['closest_order'], final_combined_results['bucket_ids'],
            final_combined_results['embedded_hs_coordinates'], final_combined_results['hc_order_val'])