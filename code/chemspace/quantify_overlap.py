# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
#
# Ref: https://forum.knime.com/t/tanimoto-similarity-using-count-based-fingerprints/12176/3
# Ref: https://pubs.acs.org/doi/full/10.1021/ci9800211
# Ref: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0069-3?optIn=false

import pandas as pd
import math
import sys

# formula
# Two embedding-vectors: X and Y both of N dimensions
#
# numerator: sum 1 through N: Xi * Yi
#
# denominator: sum Xi^2 + sum Yi^2 - numerator


def aggregate_dim (df):
    df_agg = df.groupby(['bucket_id'], as_index = False).agg({
        'ID': 'count'
    })

    df_agg = df_agg.rename (columns = {'ID': 'dim_count'})    

    return (df_agg)


def get_embedding_vector (df, n_dim):
    embedding_vector = []

    for i in range(n_dim):
        df_i = df[df['bucket_id'] == i].copy()

        if (not df_i.empty):
            embedding_vector.append(list(df_i['dim_count'])[0])
        else:
            embedding_vector.append(0)

    return (embedding_vector)
        
def process_vectors (vec1, vec2):
    value_sum = 0

    if len(vec1) != len(vec2):
        print ('[ERROR] Mismatch in vector lengths. Terminating...')
        sys.exit(-1)
    

    n_dim = len (vec1)

    for i in range(n_dim):
        value_sum += vec1[i] * vec2[i]

    return (value_sum)


drug_chembl_beginning = '../../data/app_drugs_drugbank_chembl_24_1_bms_ord_'

drug_natprod_beginning = '../../data/app_drugbank_into_hc_natprod_bms_ord_'

canvass_chembl_beginning = '../../data/canvass_chembl_24_1_bms_ord_'

canvass_natprod_beginning = '../../data/canvass_into_hc_natprod_bms_ord_'

reference_set = 'chembl'
ord_values = []
tanimoto_values = []
df_all = pd.DataFrame ()

for i in range (2, 9):
    p = math.pow(2, i) 
    n_dim = int(math.pow(p,2))
    print (n_dim)

    drug_chembl_fname = drug_chembl_beginning + str(i) + '_dim_2.tab'
    canvass_chembl_fname = canvass_chembl_beginning + str(i) + '_dim_2.tab'

    df_drug = pd.read_csv (drug_chembl_fname, sep = '\t')
    df_canvass = pd.read_csv (canvass_chembl_fname, sep = '\t')

    df_drug_agg = aggregate_dim(df_drug)
    df_canvass_agg = aggregate_dim(df_canvass)

    drug_ev = get_embedding_vector (df_drug_agg, n_dim)
    canvass_ev = get_embedding_vector (df_canvass_agg, n_dim)

    c = process_vectors (drug_ev, canvass_ev)
    a = process_vectors (drug_ev, drug_ev)
    b = process_vectors (canvass_ev, canvass_ev)
    
    T = float(c / (a + b - c))
    
    ord_values.append(int(i))
    tanimoto_values.append(T)

    print (T)

df_all = pd.DataFrame ({'dataset': 'drugbank_canvass', 'reference_set': reference_set, 'phc_order': ord_values, 'tanimoto_sim': tanimoto_values})



reference_set = 'natprod'
ord_values = []
tanimoto_values = []


for i in range (2, 6):
    p = math.pow(2, i) 
    n_dim = int(math.pow(p,2))
    print (n_dim)


    drug_natprod_fname = drug_natprod_beginning + str(i) + '_dim_2.tab'
    canvass_natprod_fname = canvass_natprod_beginning + str(i) + '_dim_2.tab'

    df_drug = pd.read_csv (drug_natprod_fname, sep = '\t')
    df_canvass = pd.read_csv (canvass_natprod_fname, sep = '\t')

    df_drug_agg = aggregate_dim(df_drug)
    df_canvass_agg = aggregate_dim(df_canvass)

    drug_ev = get_embedding_vector (df_drug_agg, n_dim)
    canvass_ev = get_embedding_vector (df_canvass_agg, n_dim)

    c = process_vectors (drug_ev, canvass_ev)
    a = process_vectors (drug_ev, drug_ev)
    b = process_vectors (canvass_ev, canvass_ev)
    
    T = float(c / (a + b - c))
 
    ord_values.append(int(i))
    tanimoto_values.append(T)


    print (T)

df_res = pd.DataFrame ({'dataset': 'drugbank_canvass', 'reference_set': reference_set, 'phc_order': ord_values, 'tanimoto_sim': tanimoto_values})

df_all = df_all.append (df_res, ignore_index = True)
df_all = df_all.reset_index (drop = True)

df_all.to_csv ('../../data/quantified_overlaps.tab', sep = '\t', index = False)


