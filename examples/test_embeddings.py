import os
import numpy as np
import pandas as pd

# Define file paths
reference_file = os.path.join(os.path.dirname(__file__), '../data/drugs_emb_hcase_chembl_original.tab')
new_file = os.path.join(os.path.dirname(__file__), '../data/drugs_emb_hcase_chembl.tab')

def load_embeddings(file_path):
    """
    Helper function to load the embeddings from the given file path.
    Assumes the file is a tab-delimited file with embeddings as numerical values.
    """
    # Load the data using pandas to handle non-numeric columns and numeric columns separately
    df = pd.read_csv(file_path, sep='\t')
    
    # Extract numeric columns (the embeddings columns should be numeric, ignoring non-numeric columns like 'id', 'structure', etc.)
    numeric_columns = df.select_dtypes(include=[np.number]).columns
    embeddings = df[numeric_columns].values  # Convert the numeric columns to a numpy array
    
    return embeddings

def test_embeddings():
    # Check if the files exist
    assert os.path.exists(reference_file), f"Reference file {reference_file} does not exist."
    assert os.path.exists(new_file), f"New embeddings file {new_file} does not exist."
    
    # Load the embeddings
    reference_embeddings = load_embeddings(reference_file)
    new_embeddings = load_embeddings(new_file)

    # Ensure the embeddings have the same shape
    assert reference_embeddings.shape == new_embeddings.shape, (
        f"Shape mismatch: reference has {reference_embeddings.shape}, "
        f"new has {new_embeddings.shape}")

    # Compare the embeddings (You can adjust this threshold based on your needs)
    difference = np.abs(reference_embeddings - new_embeddings)

    # Here we check if any values deviate significantly. Adjust the threshold as necessary.
    max_difference = np.max(difference)
    threshold = 1e-4  # Allowable threshold for deviation

    assert max_difference < threshold, f"Embeddings differ by more than {threshold}. Max difference: {max_difference}"

    print("Embeddings match within the allowable threshold.")
