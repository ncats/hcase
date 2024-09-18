# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Scope: Generate Scaffold Keys, own implementation based on paper by Peter Ertl
#
# Reference: https://pubs.acs.org/doi/10.1021/ci5001983
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
#
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html
# Ref: https://pubs.acs.org/doi/suppl/10.1021/ci5006614
# Ref: https://askubuntu.com/questions/742782/how-to-install-cpickle-on-python-3-4
# Ref: https://sourceforge.net/p/rdkit/mailman/message/35745115/
# Ref: https://github.com/rdkit/rdkit/issues/2018
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/24426410/
# Ref: https://www.rdkit.org/docs/Cookbook.html
# Ref: https://github.com/rdkit/rdkit/blob/master/Code/JavaWrappers/gmwrapper/src-test/org/RDKit/FingerprintsTests.java
# Ref: https://www.rdkit.org/docs/cppapi/structRDKit_1_1ReactionFingerprintParams.html
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.inchi.html
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html
# Ref: https://stackabuse.com/append-vs-extend-in-python-lists/
# Ref: https://www.guru99.com/python-regular-expressions-complete-tutorial.html
# Ref: https://docs.python.org/3/library/re.html
# Ref: https://stackoverflow.com/questions/28417293/sample-datasets-in-pandas
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.append.html
# Ref: https://pandas.pydata.org/pandas-docs/version/0.25/reference/api/pandas.DataFrame.from_dict.html
# Ref: https://www.geeksforgeeks.org/create-a-new-column-in-pandas-dataframe-based-on-the-existing-columns/
# Ref: https://stackoverflow.com/questions/40353519/how-to-apply-custom-function-to-pandas-data-frame-for-each-row
# Ref: https://stackoverflow.com/questions/19914937/applying-function-with-multiple-arguments-to-create-a-new-pandas-column
# Ref: https://stackoverflow.com/questions/39441484/pandas-groupby-and-aggregate-without-losing-the-column-which-was-grouped
# Ref: https://stackoverflow.com/questions/23690284/pandas-apply-function-that-returns-multiple-values-to-rows-in-pandas-dataframe
# Ref: https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
# Ref: https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
# Ref: https://stackoverflow.com/questions/23481954/turning-a-two-dimensional-array-into-a-two-column-dataframe-pandas
# Ref: https://datascience.stackexchange.com/questions/26333/convert-a-list-of-lists-into-a-pandas-dataframe
# Ref: https://www.geeksforgeeks.org/python-split-string-into-list-of-characters/
# Ref: https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html
# Ref: https://pbpython.com/pandas_dtypes.html
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.pivot.html
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.fillna.html
# Ref: https://pubs.acs.org/doi/suppl/10.1021/ci5006614
# Ref: https://askubuntu.com/questions/742782/how-to-install-cpickle-on-python-3-4
# Ref: https://sourceforge.net/p/rdkit/mailman/message/35745115/
# Ref: https://github.com/rdkit/rdkit/issues/2018
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/24426410/
# Ref: https://www.rdkit.org/docs/Cookbook.html
# Ref: https://github.com/rdkit/rdkit/blob/master/Code/JavaWrappers/gmwrapper/src-test/org/RDKit/FingerprintsTests.java
# Ref: https://www.rdkit.org/docs/cppapi/structRDKit_1_1ReactionFingerprintParams.html
# Ref: https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.Draw.SimilarityMaps.html
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
# Ref: https://stackoverflow.com/questions/19078325/naming-returned-columns-in-pandas-aggregate-function
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.fillna.html
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.reset_index.html
# Ref: https://stackoverflow.com/questions/12203901/pandas-crashes-on-repeated-dataframe-reset-index/12204428
# Ref: https://stackoverflow.com/questions/42698322/cannot-resolve-column-numeric-column-name-in-spark-dataframe
# Ref: https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
# Ref: https://stackoverflow.com/questions/14940743/selecting-excluding-sets-of-columns-in-pandas
# Ref: https://www.linkedin.com/pulse/dimensionality-reduction-using-tsne-python-deepak-kumar
# Ref: https://seaborn.pydata.org/generated/seaborn.heatmap.html
# Ref: https://medium.com/@vladbezden/how-to-set-seaborn-plot-size-in-jupyter-notebook-63ffb1415431
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.fillna.html
# Ref: https://likegeeks.com/seaborn-heatmap-tutorial/#Adjust-heatmap-font-size
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sort_values.html
# Ref: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.rank.html
# Ref: https://pandas.pydata.org/pandas-docs/version/0.24/reference/api/pandas.DataFrame.rank.html
# Ref: https://stackoverflow.com/questions/41519991/how-to-make-seaborn-heatmap-larger-normal-size
# Ref: https://datascience.stackexchange.com/questions/17540/make-seaborn-heatmap-bigger
# Ref: https://stackoverflow.com/questions/14406214/moving-x-axis-to-the-top-of-a-plot-in-matplotlib
# Ref: https://stackoverflow.com/questions/12444716/how-do-i-set-the-figure-title-and-axes-labels-font-size-in-matplotlib
# Ref: https://stackoverflow.com/questions/10998621/rotate-axis-text-in-python-matplotlib
# Ref: https://stackoverflow.com/questions/6774086/why-is-my-xlabel-cut-off-in-my-matplotlib-plo
# Ref: https://stackoverflow.com/questions/27037241/changing-the-rotation-of-tick-labels-in-seaborn-heatmap
# Ref: https://stackoverflow.com/questions/34706845/change-xticklabels-fontsize-of-seaborn-heatmap
# Ref: https://stackoverflow.com/questions/26540035/rotate-label-text-in-seaborn-factorplot/34722235
# Ref: https://stackoverflow.com/questions/2969867/how-do-i-add-space-between-the-ticklabels-and-the-axes-in-matplotlib
# Ref: https://stackoverflow.com/questions/6406368/matplotlib-move-x-axis-label-downwards-but-not-x-axis-ticks
# Ref: https://stackoverflow.com/questions/37233108/seaborn-change-font-size-of-the-colorbar
# Ref: https://seaborn.pydata.org/generated/seaborn.clustermap.html
# Ref: https://stackoverflow.com/questions/34706845/change-xticklabels-fontsize-of-seaborn-heatmap
# Ref: https://seaborn.pydata.org/examples/structured_heatmap.html
# Ref: https://stackoverflow.com/questions/49254337/how-do-i-add-a-title-to-a-seaborn-clustermap
# Ref: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.suptitle.html
# Ref: https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
# Ref: https://www.c-sharpcorner.com/article/a-complete-python-seaborn-tutorial/
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
# Ref: https://stackoverflow.com/questions/46234158/how-to-remove-x-and-y-axis-labels-in-a-clustermap
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
# Ref: https://www.rdkit.org/docs/RDKit_Book.html
# Ref: https://buildmedia.readthedocs.org/media/pdf/rdkit/release_2017_03_1/rdkit.pdf
# Ref: https://www.britannica.com/science/periodic-table
# Ref: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
# Ref: https://tel.archives-ouvertes.fr/tel-01427636v1/document
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/28157259/
# Ref: https://www.rdkit.org/UGM/2012/Evans_Papadatos_RDKit_UGM.pdf
# Ref: https://git.durrantlab.pitt.edu/jdurrant/gypsum_dl/blob/f610eb34b08868a41b9837b5fcb5a6307d9ce71f/gypsum/molvs/resonance.py
# Ref: http://rdkit.org/docs_temp/source/rdkit.Chem.rdchem.html
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html
# Ref: https://docs.python.org/2/library/math.html
# Ref: http://www.datasciencemadesimple.com/join-merge-data-frames-pandas-python/
# Ref: https://docs.python.org/3/library/functions.html#round
# Ref: https://www.daylight.com/meetings/summerschool01/course/basics/smarts.html


import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdChemReactions import CreateStructuralFingerprintForReaction, ReactionFingerprintParams, FingerprintType
import math


heteroatoms_nr = '#3,#4,#5,#7,#8,#9,#12,#13,#14,#15,#16,#17,#25,#26,#27,#28,#29,#30,#31,#32,#33,#34,#35,#45,#46,#47,#48,#50,#52,#53,#78,#80'
heteroatoms_symbol = 'Li,Be,B,N,O,F,Mg,Al,Si,P,S,Cl,Zn,As,Se,Br,Te,I,Pt,Hg,Mn,Fe,Co,Ni,Cu,Ga,Ge,Rh,Pd,Ag,Cd,Sn'
heteroatoms_reduced_nr = '#3,#4,#5,#9,#12,#13,#14,#15,#17,#25,#26,#27,#28,#29,#30,#31,#32,#33,#34,#35,#45,#46,#47,#48,#50,#52,#53,#78,#80'


def smiles2bmscaffold(smiles):
    if smiles == 'NA':
        bms = 'NA'
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            # print (Chem.MolToSmiles(mol))
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


# 1    number of ring and linker atoms    RL    20.029    7.556
def sk_num_ring_and_linker_atoms(mol):
    smiles = Chem.MolToSmiles(mol).replace('=', '-')
    bms = smiles2bmscaffold(smiles)
    bms = Chem.MolFromSmiles(bms)

    return (bms.GetNumAtoms())


# 2    number of linker atoms    L    2.518    3.481
def sk_num_linker_atoms(mol):
    smiles = Chem.MolToSmiles(mol).replace('=', '-')
    bms = smiles2bmscaffold(smiles)
    bms = Chem.MolFromSmiles(bms)

    patt = Chem.MolFromSmarts('[!R]')
    pm = bms.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 3    number of linker bonds    L    3.993    3.897
def sk_num_linker_bonds(mol):
    smiles = Chem.MolToSmiles(mol).replace('=', '-')
    bms = smiles2bmscaffold(smiles)
    bms = Chem.MolFromSmiles(bms)
    patt = Chem.MolFromSmarts('[!R][!#1]')
    ap = bms.GetSubstructMatches(patt)  # ap: atom pairs

    return (len(ap))

    """
    Alternative to sk_num_linker_bonds(mol):

    smi_mod = Chem.MolToSmiles(mol).replace('=', '-')
    m = Chem.MolFromSmiles(smi_mod)
    bms = smiles2bmscaffold (Chem.MolToSmiles(m))
    bms = Chem.MolFromSmiles (bms)
    nr_all_bonds = bms.GetNumBonds()
    bms_mod = Chem.DeleteSubstructs(bms, Chem.MolFromSmarts('[!R]'))
    nr_ring_bonds = bms_mod.GetNumBonds()


    return (nr_all_bonds - nr_ring_bonds)
    """


# 4    number of rings         3.348    3.156
def sk_num_rings(mol):
    nr = len(Chem.rdmolops.GetSSSR(mol))
    return (nr)


# 5    number of spiro atoms    R    0.031    0.193
def sk_num_spiro_atoms(mol):
    nr = Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)

    return (nr)


# 6    size of the largest ring    R    6.241    1.905
def sk_size_largest_ring(mol):
    rings = Chem.rdmolops.GetSymmSSSR(mol)
    max_size = 0
    for i in range(len(rings)):
        rs = len(list(rings[i]))
        if rs > max_size:
            max_size = rs

    return (max_size)


# 7    number of bonds in fully conjugated rings    R    13.824    6.201
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
    nr_multiple_bonds_infcr = 0    # infcr: in not fully conjugated ring
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


# 8    number of multiple bonds in not fully conjugated rings    R    0.112    0.383
# Implemented as part of Rule 7 above.


# 9    number of heteroatoms in rings    R    2.177    1.640
def sk_num_heteroatoms_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[' + heteroatoms_nr + ']~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 10    number of heteroatoms other than N, S, O in rings    R    0.003    0.061
def sk_num_heteroatoms_in_rings_no_N_O_S(mol):
    patt = Chem.MolFromSmarts('[*]~[' + heteroatoms_reduced_nr + ']~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 11    number of S ring atoms    R    0.143    0.388
def sk_num_S_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[#16]~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 12    number of O ring atoms    R    0.310    0.703
def sk_num_O_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[#8]~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 13    number of N ring atoms    R    1.721    1.507
def sk_num_N_in_rings(mol):
    patt = Chem.MolFromSmarts('[*]~[#7]~[*]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 14    number of heteroatoms    A    4.248    2.921
def sk_num_heteroatoms(mol):
    patt = Chem.MolFromSmarts('[' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 15    number of heteroatoms other than N, S, O    A    0.009    0.131
def sk_num_heteroatoms_no_N_S_O(mol):
    patt = Chem.MolFromSmarts('[' + heteroatoms_reduced_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 16    number of S atoms    A    0.289    0.540
def sk_num_S_atoms(mol):
    patt = Chem.MolFromSmarts('[#16]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 17    number of O atoms    A    1.603    1.695
def sk_num_O_atoms(mol):
    patt = Chem.MolFromSmarts('[#8]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 18    number of N atoms    A    2.347    1.789
def sk_num_N_atoms(mol):
    patt = Chem.MolFromSmarts('[#7]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 19    number of multiple linker bonds    A    0.109    0.351
# Remark (GZK): Definition of 'multiple linker' was not provided. Hence, I quantified this property as the number of bonds associated with the branched linker atom.
def sk_num_multiple_linker_bonds(mol):
    smiles = Chem.MolToSmiles(mol).replace('=', '-')
    bms = smiles2bmscaffold(smiles)
    bms = Chem.MolFromSmiles(bms)
    patt = Chem.MolFromSmarts('[D3,D4;!R][!#1]')
    ap = bms.GetSubstructMatches(patt)  # ap: atom pairs

    return (len(ap))


# 20    count of two adjacent heteroatoms    AO    0.575    1.162
def sk_num_2_adjacent_heteroatoms(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 21    count of 3 adjacent heteroatoms    AO    0.350    1.169
def sk_num_3_adjacent_heteroatoms(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][' + heteroatoms_nr + '][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 22    count of 2 heteroatoms separated by a single carbon    AO    1.804    1.953
def sk_num_2_heteroatoms_separated_by_single_carbon(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][#6][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 23    count of 2 heteroatoms separated by 2 carbons    AO    1.505    2.564
def sk_num_2_heteroatoms_separated_by_two_carbons(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + '][#6][#6][' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 24    number of double bonds with at least one heteroatom    AO    1.235    1.433
def sk_num_double_bonds_with_at_least_one_heteroatom(mol):
    patt = Chem.MolFromSmarts(
        '[' + heteroatoms_nr + ']=[' + heteroatoms_nr + ']')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 25    number of heteroatoms adjacent to a double (nonaromatic) bond    AO    1.439    1.861
def sk_num_heteroatoms_adjacent_to_non_aromatic_double_bond(mol):
    patt = Chem.MolFromSmarts(
        '[C,' + heteroatoms_symbol + ']=[C,' + heteroatoms_symbol + '][' + heteroatoms_nr + ']')
    # Maybe?:
    # patt = Chem.MolFromSmarts('[!a]=[!a][#3,#4,#5,#7,#8,#9,#12,#13,#14,#15,#16,#17,#30,#33,#34,#35,#52,#53,#78,#80]')

    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 26    count of pairs of conjugated double (nonaromatic) bonds    AO    0.094    0.380
def sk_num_conjugated_pairs_of_double_bonds_nonaromatic(mol):
    patt = Chem.MolFromSmarts('[!a]=[!a][!a]=[!a]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 27    count of pairs of adjacent branched atoms    AO    2.860    2.320
def sk_num_pairs_of_adjacent_branched_atoms(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 28    count of branched atoms separated by a single nonbranched atom    AO    1.504    1.467
def sk_num_branched_atoms_separated_by_single_nonbranched_atom(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][!#1D2][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 29    count of 3 adjacent branched atoms    AO    1.734    2.591
def sk_num_three_of_adjacent_branched_atoms(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][$([!#1]([!#1])([!#1])[!#1])][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 30    count of branched atoms separated by any 2 atoms    AO    4.294    4.409
def sk_num_branched_atoms_separated_by_any_2_atoms(mol):
    patt = Chem.MolFromSmarts(
        '[$([!#1]([!#1])([!#1])[!#1])][*][*][$([!#1]([!#1])([!#1])[!#1])]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


# 31    number of exocyclic and exolinker atoms    !R, !L    1.170    1.425
def sk_num_exocyclic_and_exolinker_atoms(mol):
    smiles = Chem.MolToSmiles(mol)
    bms = smiles2bmscaffold(smiles)
    bms = Chem.MolFromSmiles(bms)

    patt_exocyclic = Chem.MolFromSmarts('[R]=[!#1;!R]')
    patt_exolinker = Chem.MolFromSmarts('[*][!R](=[!#1;!R])[*]')

    pm_exoclyclic = bms.GetSubstructMatches(patt_exocyclic)
    nr_exoclyclic = len(pm_exoclyclic)

    pm_exolinker = bms.GetSubstructMatches(patt_exolinker)
    nr_exolinker = len(pm_exolinker)

    return (nr_exoclyclic + nr_exolinker)


# 32    number of heteroatoms with more than 2 bonds
def sk_num_heteroatoms_with_more_than_2_bonds(mol):
    patt = Chem.MolFromSmarts('[' + heteroatoms_nr + ';H0;!X1&!X2]')
    pm = mol.GetSubstructMatches(patt)
    nr = len(pm)

    return (nr)


def generate_scaffold_key(mol):
    sc_keys = []

    sc_keys.append(sk_num_ring_and_linker_atoms(mol))
    sc_keys.append(sk_num_linker_atoms(mol))
    sc_keys.append(sk_num_linker_bonds(mol))
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
    sc_keys.append(sk_num_multiple_linker_bonds(mol))
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


def smiles2scaffoldkey(smiles, trailing_inchikey=False):
    """
        Original rule-set of Scaffold Keys does not include inchikey, but I included it an optional argument, which can help deal with collisions (scaffolds of identical Scaffold Keys).
    """
    sk = ''
    sc_keys = []
    try:
        mol = Chem.MolFromSmiles(smiles)
        bms = smiles2bmscaffold(smiles)
        bms = Chem.MolFromSmiles(bms)
        inchikey = smiles2inchikey(smiles)
        sc_keys = generate_scaffold_key(mol)

        for i in range(len(sc_keys)):
            sk += (str(sc_keys[i]) + ' ')

        if trailing_inchikey:
            sk += inchikey

        sk = sk.strip()

    except:
        # print('[WARNING] SMILES: %s cannot be processed by RDKit.' % (smiles))
        sk = 'NA'

    return (sk)


def onestring(scaffold_key, has_inchikey=False):
    """
        This function can be used to convert scaffold key into a fixed-length string so it can be used to sort scaffold keys simply alphanumerically.
    """
    tmp = scaffold_key.split(' ')
    tmp2 = []
    res = ''

    if has_inchikey:
        tmp2 = tmp[:-1]
    else:
        tmp2 = tmp

    for i in range(len(tmp2)):
        dim_val = tmp2[i]
        if (len(dim_val) == 1):
            dim_val = '000' + dim_val
        elif (len(dim_val) == 2):
            dim_val = '00' + dim_val
        elif (len(dim_val) == 3):
            dim_val = '0' + dim_val
        elif (len(dim_val) > 4):
            print('[ERROR]: Current implementation does not handle values greater than 9999, please contact the developer. Terminating ...')

        res += dim_val
    res = res.strip()

    if has_inchikey:
        res += tmp[(len(tmp) - 1)]

    return (res)


def sk_distance(sk1, sk2):



    distance = 0.0
    tmp1 = sk1.split(' ')
    tmp2 = sk2.split(' ')

    # If inchikey is trailing the scaffold key then remove it:
    if len(tmp1) == 32:
        i = 1
    elif len(tmp1) == 33:
        tmp1 = tmp1[:-1]
    else:
        print('[ERROR:] Length of scaffold key is neither 32 nor 33 (inchikey-appended). Problematic scaffold key: %s. Terminating' % (sk1))

    # If inchikey is trailing the scaffold key then remove it:
    if len(tmp2) == 32:
        i = 1
    elif len(tmp2) == 33:
        tmp2 = tmp2[:-1]
    else:
        print('[ERROR:] Length of scaffold key is neither 32 nor 33 (inchikey-appended). Problematic scaffold key: %s. Terminating' % (sk2))

    sk_length = len(tmp1)

    for i in range(sk_length):
        n_1 = int(tmp1[i])
        n_2 = int(tmp2[i])
        distance += (math.pow(math.fabs(n_1 - n_2), 1.5)) / (i + 1)

    return (distance)

# smiles = 'c1ccccc1OCCCc2ccccc2'
# smiles = 'C=C(CC1CCCCC1)C1=C2NCOC3C(=C)C(=O)SC(N=C1)=C23'
# sk = smiles2scaffoldkey (smiles, trailing_inchikey = True)
# print(sk)

# print (onestring (sk, has_inchikey = True))
