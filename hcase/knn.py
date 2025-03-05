# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#


import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def is_valid (smiles):
    is_valid = False
    try:
        mol = Chem.MolFromSmiles(smiles)
        is_valid = True
    except:
        is_valid = False
    
    return (is_valid)

def get_mol (smiles):
    mol = None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        mol = None
    
    return (mol)


def get_fingerprint (mol, radius, fplength):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,radius,nBits=fplength)
    
    return (fp)

def get_Tanimoto (fp1, fp2):
    return (DataStructs.TanimotoSimilarity(fp1, fp2))

def filter_valid_mols (df, str_col):
    df['is_valid'] = df.apply(lambda x: is_valid (x[str_col]), axis = 1)
    df = df[df.is_valid == True]
    
    return (df)


def morgan_knn (df1, df2, str1_col, str2_col, id1_col, id2_col, k, radius, fplength):
    df1['knn_mol'] = df1.apply (lambda x: get_mol (x[str1_col]), axis = 1)
    df2['knn_mol'] = df2.apply (lambda x: get_mol (x[str2_col]), axis = 1)
    df1['knn_fp'] = df1.apply (lambda x: get_fingerprint (x['knn_mol'], radius, fplength), axis = 1)
    df2['knn_fp'] = df2.apply (lambda x: get_fingerprint (x['knn_mol'], radius, fplength), axis = 1)
        
    ids = list(df1[id1_col])
    first = True

    df_res = pd.DataFrame()

    for i in range(len(ids)):
        id = ids[i]
        df_query = df1[df1[id1_col] == id].copy()
        fp_query = list(df_query['knn_fp'])[0]
        df_target = df2.copy()
        df_target['sim'] = df_target.apply (lambda x: get_Tanimoto(x['knn_fp'], fp_query), axis = 1)
        df_target['knn_query_id'] = id
        df_target = df_target[df_target[id2_col] != df_target['knn_query_id']]
        df_target = df_target.sort_values (['sim'], ascending = False)
        df_target = df_target.reset_index(drop = True)
        df_target['knn_rank'] = df_target.index + 1
        df_target = df_target[df_target['knn_rank'] <= k]
        df_target['knn_color'] = (i + 1)
        if first:
            df_res = df_target
            first = False
        else:
            df_res = pd.concat([df_res, df_target], ignore_index=True)
    
    df_res['knn_target_id'] = df_res [id2_col]
    
    df_res = df_res[['knn_query_id', 'knn_target_id', 'sim', 'knn_rank', 'knn_color']].copy()
    df_res = df_res.rename(columns={'sim': 'knn_sim'})
    df_res['knn_fp_type'] = 'Morgan'
    df_res['knn_fp_radius'] = radius
    df_res['fplength'] = fplength
    df_res = df_res.sort_values(['knn_query_id', 'knn_rank'])

    df_struct_query = df1[[id1_col, str1_col]].copy()
    df_struct_query = df_struct_query.rename (columns={
        str1_col: 'knn_query_structure'
    })
    df_res = df_res.merge (df_struct_query, left_on = 'knn_query_id', right_on = id1_col, how = 'inner')



    df_struct_target = df2[[id2_col, str2_col]].copy()
    df_struct_target = df_struct_target.rename (columns={
        str2_col: 'knn_target_structure'
    })
    df_res = df_res.merge (df_struct_target, left_on = 'knn_target_id', right_on = id2_col, how = 'inner')


    return (df_res)


def compute_sim (df, str1_col, str2_col, id1_col, id2_col, radius, fplength):

    df['not_nn_mol_query'] = df.apply (lambda x: get_mol (x[str1_col]), axis = 1)
    df['not_nn_mol_target'] = df.apply (lambda x: get_mol (x[str2_col]), axis = 1)
    df['not_nn_fp_query'] = df.apply (lambda x: get_fingerprint (x['not_nn_mol_query'], radius, fplength), axis = 1)
    df['not_nn_fp_target'] = df.apply (lambda x: get_fingerprint (x['not_nn_mol_target'], radius, fplength), axis = 1)
    df['sim'] = df.apply (lambda x: get_Tanimoto(x['not_nn_fp_query'], x['not_nn_fp_target']), axis = 1)
    
    return (df)

