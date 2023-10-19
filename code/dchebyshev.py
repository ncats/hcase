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


def determine_correlation_glob (df):
    df_x = df
    df_x['key'] = 1
    df_x_paired = df_x.merge(df_x, on = 'key', how = 'inner')
    df_x_paired = df_x_paired[df_x_paired['id_x'] != df_x_paired['id_y']]
    df_x_paired['chebyshev_dist'] = df_x_paired.apply (lambda x: chebyshev_distance (x['Dim_1_x'], x['Dim_2_x'], x['Dim_1_y'], x['Dim_2_y']), axis = 1)
    df_x_paired['rank_dist'] = df_x_paired.apply (lambda x: rank_dist (x['closest_order_x'], x['closest_order_y']), axis = 1)
    df_x_dist = df_x_paired[['chebyshev_dist', 'rank_dist']].copy()
    df_x_corr = df_x_dist.corr(method ='pearson')
    corr_value = df_x_corr['rank_dist'].values[0]
    
    # Preparation for selecting next random sample
    df = df[~df['id'].isin(df_x['id'])].copy()
    
    return (df, corr_value)



def determine_correlation_local (df, sample_size = 100, random_seed = 55555):
    df_x = df.sample(n = sample_size, random_state = random_seed)
    df_x['key'] = 1
    df_x_paired = df_x.merge(df_x, on = 'key', how = 'inner')
    df_x_paired = df_x_paired[df_x_paired['id_x'] != df_x_paired['id_y']]
    df_x_paired['chebyshev_dist'] = df_x_paired.apply (lambda x: chebyshev_distance (x['Dim_1_x'], x['Dim_2_x'], x['Dim_1_y'], x['Dim_2_y']), axis = 1)
    df_x_paired['rank_dist'] = df_x_paired.apply (lambda x: rank_dist (x['closest_order_x'], x['closest_order_y']), axis = 1)
    df_x_dist = df_x_paired[['chebyshev_dist', 'rank_dist']].copy()
    df_x_corr = df_x_dist.corr(method ='pearson')
    corr_value = df_x_corr['rank_dist'].values[0]
    
    # Preparation for selecting next random sample
    df = df[~df['id'].isin(df_x['id'])].copy()
    
    return (df, corr_value)


