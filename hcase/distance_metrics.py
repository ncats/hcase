import numpy as np
import cupy as cp

def weighted_minkowski_1_5_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    array_lib = cp if use_cupy else np
    
    # Generate the weights using either cupy or numpy
    weights = 1 / array_lib.arange(1, ref_fingerprint_arr.shape[1] + 1)
    diff = np.abs(ref_fingerprint_arr[:, None, :] - new_fingerprint_arr[None, :, :]) ** 1.5
    distances = np.einsum('ijk,k->ij', diff, weights)
    
    return distances


def euclidean_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    diff = ref_fingerprint_arr[:, None, :] - new_fingerprint_arr[None, :, :]
    distances = np.linalg.norm(diff, axis=2)
    return distances


def manhattan_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    diff = np.abs(ref_fingerprint_arr[:, None, :] - new_fingerprint_arr[None, :, :])
    distances = np.sum(diff, axis=2)
    return distances

def cosine_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    dot_product = np.dot(ref_fingerprint_arr, new_fingerprint_arr.T)
    norm_ref = np.linalg.norm(ref_fingerprint_arr, axis=1)
    norm_new = np.linalg.norm(new_fingerprint_arr, axis=1)
    similarity = dot_product / (norm_ref[:, None] * norm_new[None, :])
    distances = 1 - similarity
    return distances
