from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

# Suppress all RDKit messages (including warnings and errors)
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)


heteroatoms_nr = '#3,#4,#5,#7,#8,#9,#12,#13,#14,#15,#16,#17,#25,#26,#27,#28,#29,#30,#31,#32,#33,#34,#35,#45,#46,#47,#48,#50,#52,#53,#78,#80'
heteroatoms_symbol = 'Li,Be,B,N,O,F,Mg,Al,Si,P,S,Cl,Zn,As,Se,Br,Te,I,Pt,Hg,Mn,Fe,Co,Ni,Cu,Ga,Ge,Rh,Pd,Ag,Cd,Sn'
heteroatoms_reduced_nr = '#3,#4,#5,#9,#12,#13,#14,#15,#17,#25,#26,#27,#28,#29,#30,#31,#32,#33,#34,#35,#45,#46,#47,#48,#50,#52,#53,#78,#80'

def smiles2scaffoldkey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 'NA'  # Return early if invalid molecule

    scaffold_keys = generate_scaffold_key(mol)
    if scaffold_keys == 'NA':
        return 'NA'
    sk = "".join(map(str, scaffold_keys))

    return sk.strip()

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
    sc_keys.append(sk_num_branched_atoms_separated_by_any_2_atoms(mol))

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


def smiles2bmscaffold(smiles):
    if smiles == 'NA':
        bms = 'NA'
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            bms = MurckoScaffold.GetScaffoldForMol(mol)
            bms = Chem.MolToSmiles(bms)

        except:
            bms = 'NA'

    if len(bms) == 0:
        bms = 'NA'

    return (bms)


def smiles2inchikey(smiles):
    if smiles == 'NA':
        inchi = 'NA'
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            inchi = Chem.MolToInchi(mol)
        except:
            inchi = 'NA'

    if inchi == 'NA':
        inchikey = 'NA'
    else:
        try:
            inchikey = Chem.InchiToInchiKey(inchi)
        except:
            inchikey = 'NA'

    return (inchikey)

