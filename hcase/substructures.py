import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def filter_invalid_bms(new_structures: np.ndarray, ids: np.ndarray, bms: np.ndarray):
    """Filters out invalid BMS entries (where BMS is 'NA')."""
    valid_indices = bms != 'NA'
    return new_structures[valid_indices], ids[valid_indices], bms[valid_indices]


def smiles2bmscaffold(smiles):
    if smiles == 'NA':
        bms = 'NA'
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            bms = MurckoScaffold.GetScaffoldForMol(mol)
            bms = Chem.MolToSmiles(bms, canonical=True)

        except:
            bms = 'NA'

    if len(bms) == 0:
        bms = 'NA'

    return (bms)