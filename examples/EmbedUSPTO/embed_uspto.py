import pandas as pd
import re
from tqdm.auto import tqdm
tqdm.pandas()
try:
    import cupy as cp
except:
    pass
from rdkit import Chem
from rdkit.Chem import inchi
import os
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)
import numpy as np
import time
import argparse

from hcase.fingerprints import (
    compute_substructure_batched,
    compute_fingerprints_batched,
    undo_stringify_fingerprint,
)
from hcase import embed
from hcase.plot import plot_interactive_hs

def get_unique_fingerprint_indices_per_file(all_fps_dfs):
    """
    For each file, return the indices of fingerprints unique to that file
    relative to all previous files (first file: all unique, next: unique vs previous, etc).
    Returns a list of lists of indices (one list per file).
    """
    seen = set()
    unique_indices_per_file = []
    for df in all_fps_dfs:
        fps = df['fingerprint'].apply(lambda x: tuple(map(int, x.strip('[]').replace('\n', '').split())))
        indices = []
        for idx, fp in enumerate(fps):
            if fp not in seen:
                indices.append(idx)
                seen.add(fp)
        unique_indices_per_file.append(indices)
    return unique_indices_per_file

def report_fingerprint_intersections(all_fps_dfs, filenames):
    """
    For each file, compute unique fingerprints and intersection across all files.
    Report statistics and associated substructures.
    """
    # Convert fingerprints to hashable type (tuple of ints)
    fps_sets = []
    fps_to_substructs = []
    for df in all_fps_dfs:
        # Convert stringified fingerprint to tuple
        fps = df['fingerprint'].apply(lambda x: tuple(map(int, x.strip('[]').replace('\n', '').split())))
        fps_sets.append(set(fps))
        # Map fingerprint tuple to substructure(s)
        fp2sub = {}
        for fp, sub in zip(fps, df['substructure']):
            fp2sub.setdefault(fp, set()).add(sub)
        fps_to_substructs.append(fp2sub)

    # Intersection across all files
    intersection = set.intersection(*fps_sets)
    print(f"\n=== Fingerprint Intersection Statistics ===")
    print(f"Number of files: {len(filenames)}")
    print(f"Fingerprints common to all files (intersection): {len(intersection)}")
    if intersection:
        # Collect all substructures associated with intersection fingerprints
        substructs = set()
        for fp in intersection:
            for fp2sub in fps_to_substructs:
                substructs.update(fp2sub.get(fp, []))
        substructs_list = list(substructs)
        print("Associated substructures (intersection):")
        for sub in substructs_list[:20]:
            print(f"  - {sub}")
        if len(substructs_list) > 20:
            print(f"  ... ({len(substructs_list)-20} more not shown)")

    # Unique, intersection, and "shared with some but not all" fingerprints per file
    unique_indices_per_file = get_unique_fingerprint_indices_per_file(all_fps_dfs)
    for i, (df, indices, fname) in enumerate(zip(all_fps_dfs, unique_indices_per_file, filenames)):
        # Print unique IDs, structures, substructures, and unique fingerprints
        unique_ids = set(df['ID'])
        unique_structures = set(df['structure'])
        unique_substructures = set(df['substructure'])
        unique_fingerprints = set(df['fingerprint'].apply(lambda x: tuple(map(int, x.strip('[]').replace('\n', '').split()))))
        print(f"\nFile: {fname}")
        print(f"  Number of unique IDs: {len(unique_ids)}")
        print(f"  Number of unique structures: {len(unique_structures)}")
        print(f"  Number of unique substructures: {len(unique_substructures)}")
        print(f"  Number of unique fingerprints: {len(unique_fingerprints)}")

        fps = df['fingerprint'].apply(lambda x: tuple(map(int, x.strip('[]').replace('\n', '').split())))
        fps_set = set(fps)
        # Unique: only in this file
        others = set.union(*(set(df2['fingerprint'].apply(lambda x: tuple(map(int, x.strip('[]').replace('\n', '').split())))) for j, df2 in enumerate(all_fps_dfs) if j != i))
        unique_fps = fps_set - others
        # Intersection: in all files
        intersection_fps = set.intersection(*[set(df2['fingerprint'].apply(lambda x: tuple(map(int, x.strip('[]').replace('\n', '').split())))) for df2 in all_fps_dfs])
        # Shared with some but not all: in this file, not unique, not in intersection
        shared_some_fps = fps_set - unique_fps - intersection_fps

        print(f"\nFile: {fname}")
        print(f"  Unique fingerprints: {len(unique_fps)}")
        print(f"  In intersection (common to all): {len(fps_set & intersection_fps)}")
        print(f"  Shared with some but not all: {len(shared_some_fps)}")

        if unique_fps:
            substructs = set(df.loc[fps.isin(unique_fps), 'substructure'])
            substructs_list = list(substructs)
            print("  Associated substructures (unique):")
            for sub in substructs_list[:20]:
                print(f"    - {sub}")
            if len(substructs_list) > 20:
                print(f"    ... ({len(substructs_list)-20} more not shown)")
        if shared_some_fps:
            substructs = set(df.loc[fps.isin(shared_some_fps), 'substructure'])
            substructs_list = list(substructs)
            print("  Associated substructures (shared with some but not all):")
            for sub in substructs_list[:20]:
                print(f"    - {sub}")
            if len(substructs_list) > 20:
                print(f"    ... ({len(substructs_list)-20} more not shown)")


def canonicalize_rxn_smiles(rxn_smiles: str) -> str:
    def canon_and_sort(smiles_part):
        mols = [Chem.MolFromSmiles(smi) for smi in smiles_part.split('.') if smi]
        canons = [Chem.MolToSmiles(mol, canonical=True) for mol in mols if mol is not None]
        return '.'.join(sorted(canons))

    parts = rxn_smiles.split('>')
    if len(parts) == 3:
        reactants, agents, products = parts
    elif len(parts) == 2:
        reactants, products = parts
        agents = ''
    else:
        return None

    reactants = canon_and_sort(reactants)
    agents = canon_and_sort(agents)
    products = canon_and_sort(products)

    return f"{reactants}>{agents}>{products}"


def extract_reaction_smiles_from_csv(file_path, id_prefix="USPTO"):
    # Auto-detect delimiter
    import csv
    with open(file_path, 'r', encoding='utf-8') as f:
        sample = f.read(2048)
        sniffer = csv.Sniffer()
        delimiter = sniffer.sniff(sample).delimiter
    df = pd.read_csv(file_path, delimiter=delimiter, usecols=['ReactionSmiles'], dtype={'ReactionSmiles': str})
    df = df.dropna().reset_index(drop=True)
    df['ID'] = [f"{id_prefix}{i+1:07d}" for i in range(len(df))]
    return df[['ID', 'ReactionSmiles']]


def extract_product_smiles(reaction_smiles):
    if '>' in reaction_smiles:
        product_smiles = reaction_smiles.split('>')[-1]
    else:
        product_smiles = reaction_smiles
    product_smiles = re.sub(r'\|f:[\d.,]+\|', '', product_smiles)
    return product_smiles.strip()


def smiles_to_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except:
        return None


def process(file_path, substructure_file='substructures.csv',
        fingerprint_file='fingerprints.csv',
        preprocessed_file='preprocessed_rxns.csv',
        fps_inchikey_file='df_fps_with_inchikey.csv', id_prefix='USPTO'):

    # Get absolute directory path of the current file
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # Update file paths to absolute paths relative to this script
    file_path = os.path.join(base_dir, file_path)
    substructure_file = os.path.join(base_dir, substructure_file)
    fingerprint_file = os.path.join(base_dir, fingerprint_file)
    preprocessed_file = os.path.join(base_dir, preprocessed_file)
    fps_inchikey_file = os.path.join(base_dir, fps_inchikey_file)

    # --- Preprocessing ---
    if os.path.exists(preprocessed_file):
        df_rxns = pd.read_csv(preprocessed_file)
        original_count = None  # Can't know original if only loading preprocessed
    else:
        df_rxns_raw = extract_reaction_smiles_from_csv(file_path, id_prefix)
        original_count = len(df_rxns_raw)
        print(f"  [Preprocessing] Rows after loading CSV: {len(df_rxns_raw)}")
        df_rxns = df_rxns_raw.copy()
        df_rxns['canonical_rxn'] = df_rxns['ReactionSmiles'].parallel_apply(canonicalize_rxn_smiles)
        print(f"  [Preprocessing] Rows after canonicalization (dropna): {df_rxns['canonical_rxn'].notna().sum()}")
        df_rxns = df_rxns.dropna(subset=['canonical_rxn'])

        df_rxns['product_smiles'] = df_rxns['canonical_rxn'].apply(extract_product_smiles)
        print(f"  [Preprocessing] Rows after extracting product_smiles: {len(df_rxns)}")
        df_rxns = df_rxns[df_rxns['product_smiles'].str.contains(r'^[^\.]+$')]  # single-product only
        print(f"  [Preprocessing] Rows after filtering for single-product: {len(df_rxns)}")
        df_rxns = df_rxns.drop_duplicates(subset=['ID']).reset_index(drop=True)
        print(f"  [Preprocessing] Rows after drop_duplicates on ID: {len(df_rxns)}")

        print(f"Saving preprocessed data to {preprocessed_file}")
        df_rxns.to_csv(preprocessed_file, index=False)
    

    # --- Fingerprints ---
    if os.path.exists(fps_inchikey_file):
        print(f"Loading cached fingerprint+InChIKey data from {fps_inchikey_file}")
        df_fps = pd.read_csv(fps_inchikey_file)
    else:
        if not os.path.exists(substructure_file) or not os.path.exists(fingerprint_file):
            compute_substructure_batched(df_rxns['product_smiles'].tolist(), csv_file=substructure_file, parallel=False)
            compute_fingerprints_batched(substructure_file=substructure_file, csv_file=fingerprint_file, parallel=False)

        df_fps = pd.read_csv(fingerprint_file, usecols=['structure', 'substructure', 'fingerprint'])

        structure_to_id = dict(zip(df_rxns['product_smiles'], df_rxns['ID']))
        df_fps['ID'] = df_fps['structure'].map(structure_to_id)

        print(f"  [Fingerprinting] Rows in df_fps: {len(df_fps)}")
        print(f"  [Fingerprinting] Unique IDs in df_fps: {df_fps['ID'].nunique()}")
        print(f"  [Fingerprinting] Unique structures in df_fps: {df_fps['structure'].nunique()}")

        print("Computing InChIKeys...")
        df_fps['InChIKey'] = df_fps['substructure'].parallel_apply(smiles_to_inchikey)

        print(f"Saving fingerprint+InChIKey data to {fps_inchikey_file}")
        df_fps.to_csv(fps_inchikey_file, index=False)

    # --- Embedding ---
    # (Embedding and sorting are now performed only after concatenation of all files in the main block.)
    postprocessed_count = len(df_rxns) if 'df_rxns' in locals() else None
    return original_count, postprocessed_count, df_fps



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process and embed reaction SMILES.')
    parser.add_argument('--file_paths', type=str, nargs='+', required=True, help='List of input SMILES CSV files')
    parser.add_argument('--special_file', type=str, default=None, help='File path of the special file whose InChIKeys will be highlighted')
    args = parser.parse_args()

    all_fps_dfs = []
    substructure_sets = []
    substructure_filenames = []
    special_ref_inchikeys = None
    timings = {}
    start = time.time()
    start_total = start
    timings['Preprocessing'] = time.time() - start
    for file_path in args.file_paths:
        prefix = os.path.splitext(os.path.basename(file_path))[0]
        substructure_file = f"{prefix}_substructures.csv"
        fingerprint_file = f"{prefix}_fingerprints.csv"
        preprocessed_file = f"{prefix}_preprocessed_rxns.csv"
        fps_inchikey_file = f"{prefix}_df_fps_with_inchikey.csv"

        # Run up to fingerprint+InChIKey generation
        original_count, postprocessed_count, df_fps = process(
            file_path=file_path,
            substructure_file=substructure_file,
            fingerprint_file=fingerprint_file,
            preprocessed_file=preprocessed_file,
            fps_inchikey_file=fps_inchikey_file,
            id_prefix=prefix
        )

        print(f"\n=== File: {file_path} ===")
        if original_count is not None:
            print(f"  Original input count: {original_count}")
        else:
            print(f"  Original input count: (unknown, loaded from preprocessed)")
        print(f"  Post-processed count: {postprocessed_count}")

        all_fps_dfs.append(df_fps)
        substructure_sets.append(set(df_fps['substructure']))
        substructure_filenames.append(file_path)
        # If this is the special file, extract its InChIKeys
        if args.special_file is not None and os.path.abspath(file_path) == os.path.abspath(args.special_file):
            special_ref_inchikeys = df_fps['InChIKey'].tolist()

    timings['Fingerprint + InChIKey'] = time.time() - start

    # Report fingerprint intersection and uniqueness statistics
    report_fingerprint_intersections(all_fps_dfs, substructure_filenames)

    # Concatenate only unique fingerprints (relative to previous files)
    unique_indices_per_file = get_unique_fingerprint_indices_per_file(all_fps_dfs)
    filtered_dfs = [df.iloc[indices] for df, indices in zip(all_fps_dfs, unique_indices_per_file)]
    combined_fps = pd.concat(filtered_dfs, ignore_index=True)

    # Compute special_unique_ref_inchikeys if special file is provided
    special_unique_ref_inchikeys = None
    if args.special_file is not None:
        # Find index of special file in substructure_filenames
        special_idx = None
        for i, fname in enumerate(substructure_filenames):
            if os.path.abspath(fname) == os.path.abspath(args.special_file):
                special_idx = i
                break
        if special_idx is not None:
            df_special = all_fps_dfs[special_idx]
            unique_idxs = unique_indices_per_file[special_idx]
            special_unique_ref_inchikeys = df_special.iloc[unique_idxs]['InChIKey'].tolist()

    # Embedding and plotting on the combined data
    target_fingerprints = undo_stringify_fingerprint(combined_fps['fingerprint'].values)
    target_ids = np.array(combined_fps['substructure'].values)
    parent_target_ids = np.array(combined_fps['structure'].values)
    parent_reaction_ids = np.array(combined_fps['ID'].values)
    target_inchikeys = np.array(combined_fps['InChIKey'].values)

    sort_idx = np.lexsort(target_fingerprints.T[::-1])
    target_fingerprints = target_fingerprints[sort_idx]
    target_ids = target_ids[sort_idx]
    parent_target_ids = parent_target_ids[sort_idx]
    parent_reaction_ids = parent_reaction_ids[sort_idx]
    target_inchikeys = target_inchikeys[sort_idx]

    ref_fingerprints, unique_indices = np.unique(target_fingerprints, axis=0, return_index=True)
    ref_ids = target_ids[unique_indices]
    parent_ref_ids = parent_target_ids[unique_indices]
    ref_inchikeys = target_inchikeys[unique_indices]
    try:
        use_cupy = cp.cuda.is_available()
    except:
        use_cupy = False




    embedded, ref_only_embedded, hc_orders, ref_ids, closest_indices_all = embed(
        target_fingerprints,
        ref_fingerprints=ref_fingerprints,
        ref_ids=ref_ids,
        target_ids=target_ids,
        n_dim=2,
        use_cupy=use_cupy,
        target_batch_size=1000,
        ref_batch_size=40000,
        max_hc_order_only=True
    )

    timings['Embedding'] = time.time() - start

    # --- Total ---
    timings['Total'] = time.time() - start_total

    print("\n=== Timing Summary ===")
    for step, secs in timings.items():
        print(f"{step:<25}: {secs:.2f} seconds")

    # Print all unique and common special InChIKeys before plotting
    if special_ref_inchikeys is not None:
        unique_set = set(special_unique_ref_inchikeys) if special_unique_ref_inchikeys is not None else set()
        all_set = set(special_ref_inchikeys)
        common_set = all_set - unique_set
        print("\nSpecial Unique InChIKeys:")
        for ik in sorted(unique_set):
            print(f"  {ik}")
        print("\nSpecial Common InChIKeys:")
        for ik in sorted(common_set):
            print(f"  {ik}")

    plot_interactive_hs(
        embedded,
        target_ids,
        hc_orders,
        ref_ids,
        parent_target_ids,
        parent_ref_ids,
        closest_indices_all,
        parent_reaction_ids,
        target_inchikeys,
        ref_inchikeys,
        special_ref_inchikeys=special_ref_inchikeys,
        special_unique_ref_inchikeys=special_unique_ref_inchikeys
    )
