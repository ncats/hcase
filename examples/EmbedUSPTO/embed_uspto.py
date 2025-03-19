import pandas as pd
import re
from tqdm.auto import tqdm
tqdm.pandas()
import cupy as cp
from hcase.fingerprints import compute_substructure_batched, compute_fingerprints_batched, undo_stringify_fingerprint
from hcase import embed

def extract_reaction_smiles_from_csv(file_path):
    """
    Extracts all ReactionSmiles from a CSV file.
    
    Args:
        file_path (str): The path to the CSV file.
        
    Returns:
        list: A list of extracted ReactionSmiles strings.
    """
    df = pd.read_csv(file_path, delimiter='\t')
    return df['ReactionSmiles'].tolist()


def extract_product_smiles(reaction_smiles):
    """
    Extracts the product SMILES from a reaction SMILES string by returning
    everything after the last '>' symbol, and removes any unwanted annotations 
    like '|f:0.1,2.3,4.5,8.9|'.
    
    Args:
        reaction_smiles (str): The reaction SMILES string.
        
    Returns:
        str: The cleaned product SMILES string.
    """
    if '>' in reaction_smiles:
        # Extract everything after the last '>'
        product_smiles = reaction_smiles.split('>')[-1]
    else:
        # If no '>' is found, return the original string (could be a single molecule)
        product_smiles = reaction_smiles

    # Remove any unwanted parts such as '|f:0.1,2.3,4.5,8.9|'
    product_smiles = re.sub(r'\|f:[\d.,]+\|', '', product_smiles)
    
    return product_smiles.strip()


file_name = '2001_Sep2016_USPTOapplications_smiles.rsmi'
substructure_file = 'substructures.csv'
fingerprint_file = 'fingerprints.csv'
skip_fingerprint_generation = True


if not skip_fingerprint_generation:

    ref_smiles = extract_reaction_smiles_from_csv(file_name)
    product_smiles = [extract_product_smiles(reaction_smiles) for reaction_smiles in ref_smiles]

    # at some point computing BMS was getting "Killed" signal so turned off parrallelization
    compute_substructure_batched(product_smiles, csv_file=substructure_file, parallel=False)
    compute_fingerprints_batched(substructure_file=substructure_file, csv_file=fingerprint_file, parallel=False)



df_fingerprints = pd.read_csv(fingerprint_file, usecols=['fingerprint'])
ref_fingerprints = df_fingerprints['fingerprint'].values
ref_fingerprints = undo_stringify_fingerprint(ref_fingerprints)
use_cupy = cp.cuda.is_available()
embedded_data = embed(
    ref_fingerprints,
    n_dim=2,
    use_cupy=use_cupy,
    target_batch_size=1000, 
    ref_batch_size=50000,  
    max_hc_order_only=True
)
