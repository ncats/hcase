# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
#
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.chebyshev.html
# Ref: https://www.geeksforgeeks.org/python-pandas-dataframe-corr/

import scipy
from scipy import spatial
from scipy.spatial import distance
from scipy.spatial.distance import chebyshev
import math
import numpy as np
import pandas as pd

def chebyshev_distance (x_dim1, x_dim_2, y_dim_1, y_dim_2):
    x = (x_dim1, x_dim_2)
    y = (y_dim_1, y_dim_2)
    d = chebyshev(x, y)
    
    return (d)

def rank_dist (rank1, rank2):
    return (math.fabs(rank1-rank2))


def determine_correlation (df):
    df_x = df.sample(n = 100, random_state = 55555)
    df_x['key'] = 1
    df_x_paired = df_x.merge(df_x, on = 'key', how = 'inner')
    df_x_paired = df_x_paired[df_x_paired['ID_x'] != df_x_paired['ID_y']]
    df_x_paired['chebyshev_dist'] = df_x_paired.apply (lambda x: chebyshev_distance (x['Dim_1_x'], x['Dim_2_x'], x['Dim_1_y'], x['Dim_2_y']), axis = 1)
    df_x_paired['rank_dist'] = df_x_paired.apply (lambda x: rank_dist (x['closest_order_x'], x['closest_order_y']), axis = 1)
    df_x_dist = df_x_paired[['chebyshev_dist', 'rank_dist']].copy()
    df_x_corr = df_x_dist.corr(method ='pearson')
    corr_value = df_x_corr['rank_dist'].values[0]
    
    # Preparation for selecting next random sample
    df = df[~df['ID'].isin(df_x['ID'])].copy()
    
    return (df, corr_value)

df_all = pd.DataFrame()

dataset_name = 'drugbank_chembl_bms'
dataset_fname_begin = '../../data/app_drugs_drugbank_chembl_24_1_bms_ord_'    
corr_mean_values = []
corr_std_values = []
ord_values = []

for i in range(2,9):
    corr_values = []
    ord_values.append(i)

    fname = dataset_fname_begin + str(i) + '_dim_2.tab'
    print ("[*] Processing dataset: %s" %(fname))
    df = pd.read_csv(fname, sep = '\t')

    while (df.shape[0] >= 100):
        (df, cv) = determine_correlation (df)
        corr_values.append(cv)
        #print (cv)

    #print (corr_values)

    corr_values = np.array(corr_values)
    corr_mean_values.append(corr_values.mean())
    corr_std_values.append(corr_values.std())

df_res = pd.DataFrame({'dataset': dataset_name, 'phc_order': ord_values, 'pearson.corr.mean': corr_mean_values, 'pearson.corr.std': corr_std_values})
df_all = df_res

dataset_name = 'drugbank_natprod_bms'
dataset_fname_begin = '../../data/app_drugbank_into_hc_natprod_bms_ord_'    
corr_mean_values = []
corr_std_values = []
ord_values = []

for i in range(2,6):
    corr_values = []
    ord_values.append(i)

    fname = dataset_fname_begin + str(i) + '_dim_2.tab'
    print ("[*] Processing dataset: %s" %(fname))

    df = pd.read_csv(fname, sep = '\t')

    while (df.shape[0] >= 100):
        (df, cv) = determine_correlation (df)
        corr_values.append(cv)
        #print (cv)

    #print (corr_values)

    corr_values = np.array(corr_values)
    corr_mean_values.append(corr_values.mean())
    corr_std_values.append(corr_values.std())

df_res = pd.DataFrame({'dataset': dataset_name, 'phc_order': ord_values, 'pearson.corr.mean': corr_mean_values, 'pearson.corr.std': corr_std_values})

df_all = df_all.append (df_res, ignore_index = True)

df_all = df_all.reset_index (drop = True)
df_all.to_csv ('../../data/chebyshev_stat.tab', sep = '\t', index = False)

print ('[Done.]')
