from rdkit import Chem
from rdkit import RDLogger
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from tqdm import tqdm
import os


from hcase.substructures import smiles2bmscaffold

# Suppress all RDKit messages (including warnings and errors)
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)


heteroatoms_nr = '#3,#4,#5,#7,#8,#9,#12,#13,#14,#15,#16,#17,#25,#26,#27,#28,#29,#30,#31,#32,#33,#34,#35,#45,#46,#47,#48,#50,#52,#53,#78,#80'
heteroatoms_symbol = 'Li,Be,B,N,O,F,Mg,Al,Si,P,S,Cl,Zn,As,Se,Br,Te,I,Pt,Hg,Mn,Fe,Co,Ni,Cu,Ga,Ge,Rh,Pd,Ag,Cd,Sn'
heteroatoms_reduced_nr = '#3,#4,#5,#9,#12,#13,#14,#15,#17,#25,#26,#27,#28,#29,#30,#31,#32,#33,#34,#35,#45,#46,#47,#48,#50,#52,#53,#78,#80'

def initialize_pandarallel(num_workers=None):
    if num_workers is None:
        # Detect available CPUs
        num_cores = os.cpu_count()
        
        # Leave 10% of the cores free
        num_workers = max(1, int(num_cores * 0.9))  # 90% of available cores

    # Initialize pandarallel with the calculated number of workers
    pandarallel.initialize(nb_workers=num_workers, progress_bar=True)

initialize_pandarallel()

def compute_substructure(structures: np.ndarray, method=smiles2bmscaffold):
    """Computes substructures (default: BMS scaffold) and filters invalid ones."""
    df = pd.DataFrame({'structure': structures})
    df['substructure'] = df['structure'].parallel_apply(method)
    
    # Filter out invalid substructures
    df_valid = df[df['substructure'] != 'NA'].copy()
    return df_valid

def compute_substructure_batched(
    structures,
    method=smiles2bmscaffold,
    batch_size=100000,
    csv_file='substructures.csv',
    parallel=True,
    deduplicate=True
):
    """Computes substructures in batches, appending to a CSV file and resuming from previous progress."""

    # Check existing results
    processed_structures = set()
    if os.path.exists(csv_file):
        df_existing = pd.read_csv(csv_file, usecols=['structure'])
        processed_structures = set(df_existing['structure'])

    # Filter only new structures
    structures = [s for s in structures if s not in processed_structures]

    # If everything is processed, exit early
    if not structures:
        return

    # Process in batches
    first_batch = not os.path.exists(csv_file)
    batches = [structures[i:i + batch_size] for i in range(0, len(structures), batch_size)]

    for batch in tqdm(batches, desc="Processing substructure batches", unit="batch"):
        df = pd.DataFrame({'structure': batch})
        if parallel:
            df['substructure'] = df['structure'].parallel_apply(method)
        else:
            df['substructure'] = df['structure'].apply(method)
        df = df[df['substructure'] != 'NA']

        # Append results, preserving existing data
        df.to_csv(csv_file, mode='w' if first_batch else 'a', header=first_batch, index=False)
        first_batch = False

    # Deduplicate structure→substructure, keeping only the first instance
    if deduplicate:
        df_all = pd.read_csv(csv_file)
        df_dedup = df_all.drop_duplicates(subset='substructure', keep='first')
        df_dedup.to_csv(csv_file, index=False)


def compute_fingerprints(df: pd.DataFrame, method=lambda x: smiles2scaffoldkey(x, return_as_string=True)):
    """Computes fingerprints from substructures and filters invalid ones."""
    df['fingerprint'] = df['substructure'].parallel_apply(method)
    
    # Filter out invalid fingerprints
    df_valid = df[df['fingerprint'] != 'NA'].copy()
    return df_valid

def compute_fingerprints_batched(substructure_file='substructures.csv',
                                 method=lambda x: smiles2scaffoldkey(x, return_as_string=True),
                                 batch_size=10000,
                                 csv_file='fingerprints.csv',
                                 parallel=True):
    # Load processed substructures iteratively to save memory
    processed_substructures = set()
    if os.path.exists(csv_file):
        for chunk in pd.read_csv(csv_file, usecols=['substructure'], dtype=str, chunksize=10000):
            processed_substructures.update(chunk['substructure'].values)

    # Get total line count (excluding header)
    with open(substructure_file) as f:
        file_size = sum(1 for _ in f) - 1  

    # Adjust total_chunks based on batch_size
    total_chunks = (file_size // batch_size) + (1 if file_size % batch_size else 0)

    reader = pd.read_csv(substructure_file, dtype=str, chunksize=batch_size)
    
    first_batch = not os.path.exists(csv_file)

    with tqdm(total=total_chunks, desc="Processing substructure chunks", unit="chunk") as pbar:
        for df_substructures in reader:
            # Filter only new substructures
            new_substructures = df_substructures[~df_substructures['substructure'].isin(processed_substructures)]
            new_substructures = new_substructures.copy()
            if new_substructures.empty:
                pbar.update(1)
                continue  # Skip empty batches

            # Compute fingerprints
            if parallel:
                new_substructures['fingerprint'] = new_substructures['substructure'].parallel_apply(method)
            else:
                new_substructures['fingerprint'] = new_substructures['substructure'].apply(method)

            # Filter out 'NA' values
            new_substructures = new_substructures[new_substructures['fingerprint'] != 'NA']

            # Append results
            new_substructures.to_csv(csv_file, mode='w' if first_batch else 'a', header=first_batch, index=False)
            first_batch = False
            
            # Update progress bar
            pbar.update(1)
    

def prepare_fingerprints(structures: np.ndarray):
    """Prepares scaffold-based fingerprints for embedding, filtering out invalid entries."""
    df_valid_substructures = compute_substructure(structures)
    df_valid_fingerprints = compute_fingerprints(df_valid_substructures)

    # Extract and return as NumPy arrays
    return (
        df_valid_fingerprints['fingerprint'].values,
        df_valid_fingerprints['substructure'].values,
        df_valid_fingerprints['structure'].values,
        df_valid_fingerprints.index  # Keep index for filtering `ids`
    )


def order_scaffolds(df: pd.DataFrame) -> pd.DataFrame:

    # Apply the function in parallel
    df['sk'] = df['structure'].parallel_apply(lambda x: smiles2scaffoldkey(x, return_as_string=True))

    # Filter out rows where the scaffold key is 'NA'
    df = df[df['sk'] != 'NA']

    # Apply sk_one function (structure key variant) in parallel or sequentially
    df['sk_one'] = df['sk'].apply(stringify_fingerprint)

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


    return df

def undo_stringify_fingerprint(fingerprints, remove_trailing_inchikey=False):
    if remove_trailing_inchikey:
        fingerprints = [fingerprint.rsplit(' ', 1)[0] for fingerprint in fingerprints]

    # Split each string by spaces and convert to integers
    split_fingerprints = np.array([list(map(int, fingerprint.split())) for fingerprint in fingerprints], dtype=int)
    return split_fingerprints

def stringify_fingerprint(scaffold_key):
    return ' '.join(map(str, scaffold_key))

def stringify_fingerprint_all(scaffold_keys):
    # Convert the input list of lists into a NumPy array
    scaffold_array = np.array(scaffold_keys, dtype=int)
    
    # Convert each number to a string
    string_array = scaffold_array.astype(str)
    
    # Join each row with spaces to form the fingerprint string
    scaffold_strings = np.array([' '.join(row) for row in string_array])
    
    return scaffold_strings


def smiles2scaffoldkey(smiles, return_as_string=False):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 'NA'  # Return early if invalid molecule

    scaffold_keys = generate_scaffold_key(mol)

    if scaffold_keys == 'NA':
        return 'NA'
    if return_as_string:
        scaffold_keys = stringify_fingerprint(scaffold_keys)

    return scaffold_keys

def generate_scaffold_key(mol):
    sc_keys = []

    smiles = Chem.MolToSmiles(mol).replace('=', '-')
    bms = smiles2bmscaffold(smiles)
    bms = Chem.MolFromSmiles(bms)
    if bms == None:
        return "NA"

    sc_keys.append(sk_num_ring_and_linker_atoms(bms))
    sc_keys.append(sk_num_linker_atoms(bms))
    sc_keys.append(sk_num_linker_bonds(bms))
    sc_keys.append(sk_num_rings(mol))
    sc_keys.append(sk_num_spiro_atoms(mol))
    sc_keys.append(sk_size_largest_ring(mol))

    res = sk_conjugation(mol)
    sc_keys.append(res[0])
    sc_keys.append(res[1])

    sc_keys.append(sk_num_heteroatoms_in_rings(mol))
    sc_keys.append(sk_num_heteroatoms_in_rings_no_N_O_S(mol))
    sc_keys.append(sk_num_S_in_rings(mol))
    sc_keys.append(sk_num_O_in_rings(mol))
    sc_keys.append(sk_num_N_in_rings(mol))
    sc_keys.append(sk_num_heteroatoms(mol))
    sc_keys.append(sk_num_heteroatoms_no_N_S_O(mol))
    sc_keys.append(sk_num_S_atoms(mol))
    sc_keys.append(sk_num_O_atoms(mol))
    sc_keys.append(sk_num_N_atoms(mol))
    sc_keys.append(sk_num_multiple_linker_bonds(bms))
    sc_keys.append(sk_num_2_adjacent_heteroatoms(mol))
    sc_keys.append(sk_num_3_adjacent_heteroatoms(mol))
    sc_keys.append(sk_num_2_heteroatoms_separated_by_single_carbon(mol))
    sc_keys.append(sk_num_2_heteroatoms_separated_by_two_carbons(mol))
    sc_keys.append(sk_num_double_bonds_with_at_least_one_heteroatom(mol))
    sc_keys.append(sk_num_heteroatoms_adjacent_to_non_aromatic_double_bond(mol))
    sc_keys.append(sk_num_conjugated_pairs_of_double_bonds_nonaromatic(mol))
    sc_keys.append(sk_num_pairs_of_adjacent_branched_atoms(mol))
    sc_keys.append(sk_num_branched_atoms_separated_by_single_nonbranched_atom(mol))
    sc_keys.append(sk_num_three_of_adjacent_branched_atoms(mol))
    sc_keys.append(sk_num_branched_atoms_separated_by_any_2_atoms(mol))
    sc_keys.append(sk_num_exocyclic_and_exolinker_atoms(mol))
    sc_keys.append(sk_num_heteroatoms_with_more_than_2_bonds(mol))

    return (sc_keys)

def sk_num_ring_and_linker_atoms(bms):

    return (bms.GetNumAtoms())

def sk_num_linker_atoms(bms):

    patt = Chem.MolFromSmarts('[!R]')
    pm = bms.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_linker_bonds(bms):

    patt = Chem.MolFromSmarts('[!R][!#1]')
    ap = bms.GetSubstructMatches(patt)

    return (len(ap))

def sk_num_rings(mol):
    rings = Chem.rdmolops.GetSSSR(mol)
    
    return len(rings)

def sk_num_spiro_atoms(mol):
    nr = Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)

    return (nr)

def sk_size_largest_ring(mol):
    rings = Chem.rdmolops.GetSymmSSSR(mol)
    max_size = 0
    for i in range(len(rings)):
        rs = len(list(rings[i]))
        if rs > max_size:
            max_size = rs

    return (max_size)

def is_ring_fully_conjugated(ring, suppl):
    first_idx_conj = -1
    first = True

    for i in range(len(ring)):
        atom_idx = ring[i]
        atom_conj_grp_idx = Chem.ResonanceMolSupplier.GetAtomConjGrpIdx(suppl, atom_idx)

        if first:
            first_idx_conj = atom_conj_grp_idx

            # If atom is not in conjugated system, the conjugated group index is
            # supposed to be -1, however, I experienced this value to be a large
            # positive integer instead of -1 when doing a negative control (atom
            # not in conjugated system). Therefore both a large value and -1 are
            # checked as indication of atom not being part of conjugated system.

            if first_idx_conj > 99999 or first_idx_conj == -1:
                return (False)

            first = False

        else:
            if atom_conj_grp_idx != first_idx_conj:
                return (False)

    return (True)


def get_num_multiple_bond_in_ring(mol, ring):
    atom_i = -1
    atom_j = -1
    bond = None
    nr_multiple_bonds = 0
    multiple_bond_types = [2, 3, 12]

    for i in range(len(ring)):
        atom_i = ring[i]

        for j in range(0, i):
            atom_j = ring[j]
            bond = mol.GetBondBetweenAtoms(atom_i, atom_j)
            if bond is not None:
                bond_type = Chem.Bond.GetBondType(bond)
                if int(bond_type) in multiple_bond_types:
                    nr_multiple_bonds += 1

    return (nr_multiple_bonds)


def sk_conjugation(mol):
    extracted_rings = []
    nr_conjugated_rings = 0
    nr_multiple_bonds_infcr = 0
    suppl = Chem.ResonanceMolSupplier(mol)
    rings = Chem.GetSymmSSSR(mol)

    for i in range(len(rings)):
        extracted_rings.append(list(rings[i]))

    for ring in extracted_rings:

        if is_ring_fully_conjugated(ring, suppl):
            nr_conjugated_rings += 1
        else:
            nr_multiple_bonds_infcr += get_num_multiple_bond_in_ring(mol, ring)

    return ((nr_conjugated_rings, nr_multiple_bonds_infcr))

def sk_num_heteroatoms_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[' + heteroatoms_nr + ']~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_heteroatoms_in_rings_no_N_O_S(mol):
    patt = Chem.MolFromSmarts('[*]~[' + heteroatoms_reduced_nr + ']~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_S_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[#16]~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_O_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[#8]~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_N_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[#7]~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_heteroatoms(mol):
    patt = Chem.MolFromSmarts('[' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_heteroatoms_no_N_S_O(mol):
    patt = Chem.MolFromSmarts('[' + heteroatoms_reduced_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


def sk_num_S_atoms(mol):
    patt = Chem.MolFromSmarts('[#16]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_O_atoms(mol):
    patt = Chem.MolFromSmarts('[#8]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_N_atoms(mol):
    patt = Chem.MolFromSmarts('[#7]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_multiple_linker_bonds(bms):
    patt = Chem.MolFromSmarts('[D3,D4;!R][!#1]')
    ap = bms.GetSubstructMatches(patt)

    return (len(ap))

def sk_num_2_adjacent_heteroatoms(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_3_adjacent_heteroatoms(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][' + heteroatoms_nr + '][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_2_heteroatoms_separated_by_single_carbon(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][#6][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_2_heteroatoms_separated_by_two_carbons(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][#6][#6][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_double_bonds_with_at_least_one_heteroatom(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + ']=[' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_heteroatoms_adjacent_to_non_aromatic_double_bond(mol):
    patt = Chem.MolFromSmarts(
        '[C,' + heteroatoms_symbol + ']=[C,' + heteroatoms_symbol + '][' + heteroatoms_nr + ']')

    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_conjugated_pairs_of_double_bonds_nonaromatic(mol):
    patt = Chem.MolFromSmarts('[!a]=[!a][!a]=[!a]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_pairs_of_adjacent_branched_atoms(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


def sk_num_branched_atoms_separated_by_single_nonbranched_atom(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][!#1D2][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)

def sk_num_three_of_adjacent_branched_atoms(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][$([!#1]([!#1])([!#1])[!#1])][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


def sk_num_branched_atoms_separated_by_any_2_atoms(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][*][*][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


def sk_num_exocyclic_and_exolinker_atoms(mol):
    smiles = Chem.MolToSmiles (mol)
    bms = smiles2bmscaffold (smiles)
    bms = Chem.MolFromSmiles (bms)

    patt_exocyclic = Chem.MolFromSmarts('[R]=[!#1;!R]')
    patt_exolinker = Chem.MolFromSmarts('[*][!R](=[!#1;!R])[*]')

    pm_exoclyclic = bms.GetSubstructMatches(patt_exocyclic)
    nr_exoclyclic = len(pm_exoclyclic)

    pm_exolinker = bms.GetSubstructMatches(patt_exolinker)
    nr_exolinker = len(pm_exolinker)

    return (nr_exoclyclic + nr_exolinker)


def sk_num_heteroatoms_with_more_than_2_bonds(mol):
    patt = Chem.MolFromSmarts('[' + heteroatoms_nr + ';H0;!X1&!X2]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


