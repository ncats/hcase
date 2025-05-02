import numpy as np

try:
    import cupy as cp
except ImportError:
    cp = None

def weighted_minkowski_1_5_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    if use_cupy and cp is not None:
        array_lib = cp
    else:
        array_lib = np
    
    # Generate the weights using either cupy or numpy
    weights = 1 / array_lib.arange(1, ref_fingerprint_arr.shape[1] + 1)
    diff = np.abs(ref_fingerprint_arr[:, None, :] - new_fingerprint_arr[None, :, :]) ** 1.5
    distances = np.einsum('ijk,k->ij', diff, weights)
    
    return distances


def euclidean_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    array_lib = cp if use_cupy and cp is not None else np
    diff = array_lib.subtract(ref_fingerprint_arr[:, None, :], new_fingerprint_arr[None, :, :])
    distances = array_lib.linalg.norm(diff, axis=2)
    return distances


def manhattan_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    array_lib = cp if use_cupy and cp is not None else np
    diff = array_lib.abs(array_lib.subtract(ref_fingerprint_arr[:, None, :], new_fingerprint_arr[None, :, :]))
    distances = array_lib.sum(diff, axis=2)
    return distances

def cosine_distance(new_fingerprint_arr: np.ndarray, ref_fingerprint_arr: np.ndarray, use_cupy = False) -> np.ndarray:
    array_lib = cp if use_cupy and cp is not None else np
    dot_product = array_lib.dot(ref_fingerprint_arr, new_fingerprint_arr.T)
    norm_ref = array_lib.linalg.norm(ref_fingerprint_arr, axis=1)
    norm_new = array_lib.linalg.norm(new_fingerprint_arr, axis=1)
    similarity = dot_product / (norm_ref[:, None] * norm_new[None, :])
    distances = 1 - similarity
    return distances
