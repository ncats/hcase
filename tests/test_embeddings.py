import os
import numpy as np
import pandas as pd

# Define file paths
reference_file = os.path.join(os.path.dirname(__file__), '../data/drugs_emb_hcase_chembl_original.tab')
new_file = os.path.join(os.path.dirname(__file__), '../data/drugs_emb_hcase_chembl.tab')

def load_embeddings(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]  # Drop 'Unnamed' columns
    return df  

def test_embeddings():
    assert os.path.exists(reference_file), f"Reference file {reference_file} does not exist."
    assert os.path.exists(new_file), f"New embeddings file {new_file} does not exist."

    reference_df = load_embeddings(reference_file)
    new_df = load_embeddings(new_file)

    common_ids = np.intersect1d(reference_df['id'].values, new_df['id'].values)

    reference_subset = reference_df[reference_df['id'].isin(common_ids)][['id', 'hc_order', 'Dim_1', 'Dim_2']]
    new_subset = new_df[new_df['id'].isin(common_ids)][['id', 'hc_order', 'Dim_1', 'Dim_2']]

    reference_subset_sorted = reference_subset.sort_values(by=['id', 'hc_order'])
    new_subset_sorted = new_subset.sort_values(by=['id', 'hc_order'])

    assert reference_subset_sorted.shape == new_subset_sorted.shape, (
        f"Shape mismatch: reference has {reference_subset_sorted.shape}, "
        f"new has {new_subset_sorted.shape}"
    )

    reference_embeddings = reference_subset_sorted[['Dim_1', 'Dim_2']].values
    new_embeddings = new_subset_sorted[['Dim_1', 'Dim_2']].values

    distance = np.linalg.norm(reference_embeddings - new_embeddings, axis=1)  

    threshold = 1e-4
    exceeding_indices = np.where(distance >= threshold)[0]  # Find all indices exceeding threshold

    if len(exceeding_indices) > 0:
        print(f"\n{len(exceeding_indices)} embeddings exceed the threshold of {threshold}.\n")

        for idx in exceeding_indices:
            problematic_id = reference_subset_sorted['id'].iloc[idx]
            problematic_hc_order = reference_subset_sorted['hc_order'].iloc[idx]
            original_value = reference_embeddings[idx]
            new_value = new_embeddings[idx]
            max_distance = distance[idx]

            print(f"ID: {problematic_id}, hc_order: {problematic_hc_order}")
            print(f"Max Euclidean distance: {max_distance}")
            print(f"Original embedding: {original_value}")
            print(f"New embedding: {new_value}")
            
            # Retrieve full rows before filtering
            original_full_row = reference_df[
                (reference_df['id'] == problematic_id) & 
                (reference_df['hc_order'] == problematic_hc_order)
            ].iloc[0]

            new_full_row = new_df[
                (new_df['id'] == problematic_id) & 
                (new_df['hc_order'] == problematic_hc_order)
            ].iloc[0]

            print("\nOriginal full row (before filtering):")
            for col, value in original_full_row.items():
                print(f"{col}: {value}")

            print("\nNew full row (before filtering):")
            for col, value in new_full_row.items():
                print(f"{col}: {value}")
            
            print("-" * 50)  # Separator for readability

        raise AssertionError(f"{len(exceeding_indices)} embeddings differ by more than {threshold}.")
    
    print("Embeddings match within the allowable threshold.")

