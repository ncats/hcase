# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NIH/NCATS)
#
#
# Ref: https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jaccard.html#scipy.spatial.distance.jaccard
# Ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/24426410/
# Ref: https://python-graph-gallery.com/197-available-color-palettes-with-matplotlib/
# Ref: https://stackoverflow.com/questions/57568311/matplotlib-scatter-issue-with-python-3-x
# Ref: https://www.science-emergence.com/Articles/How-to-create-a-scatter-plot-with-several-colors-in-matplotlib-/
# Ref: https://www.pluralsight.com/guides/choosing-color-palettes
# Ref: https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf
# Ref: https://htmlcolorcodes.com/color-picker/
#
#

import sys
import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


import numpy as np

import sklearn
from sklearn.manifold import TSNE

import math

from scaffold_keys import smiles2bmscaffold, smiles2scaffoldkey, sk_distance

from knn import get_mol, get_fingerprint, is_valid




def to_bitstring (fp):
    fpstr = 'NA'
    try:
        fpstr = fp.ToBitString()
    except:
        fpstr = 'NA'
    
    return (fpstr)


    
    

def fp_gen_with_errohandling (mol, fp_radius = 3, fp_length = 2048):
    fp = None
    
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, fp_radius, fp_length)
    
    except:
        fp = None
    
    return (fp)


def generate_fp (df, fp_radius = 3, fp_length = 2048):
    df['is_valid'] = df.apply (lambda x: is_valid(x['structure']), axis = 1)
    
    df = df[df['is_valid'] == True].copy()
    df['mol'] = df.apply(lambda x: get_mol (x['structure']), axis = 1)
    
    #print (df.dtypes)
    #print(df.columns)
    #print (df.head)
    
    df['fp'] = df.apply (lambda x: fp_gen_with_errohandling(x['mol'], fp_radius, fp_length), axis = 1)
    
    #print (df.columns)
    
    df = df[df['fp'] != None].copy()
    df['fp_str'] = df.apply (lambda x: to_bitstring(x['fp']), axis = 1)
    df = df[df['fp_str'] != 'NA'].copy()
    #print (df)
    
    return (df)


def get_fp_np_array (fps):
    first = True
    fpl = 0
    all_fp = []
    fp_array = []
    for i in range(len(fps)):
        fp_array = []
        fp = fps[i]
        if first:
            fpl = len(fp)
            first = False
        else:
            if len(fp) != fpl:
               print ('[ERROR] Fingerprint length mismatch. Terminating ...')
               sys.exit (-1)
        
        for j in range(len(fp)):
            fp_array.append(int(fp[j]))
        
        all_fp.append(fp_array)

    all_fp = np.array(all_fp)

    #print (all_fp.shape)

    return (all_fp)
    
    
def perform_morgan_fp_embedding (df_structures, tsne_model, fp_radius = 3, fp_length = 2048):
    df = generate_fp (df_structures, fp_radius, fp_length)


    #print (df.head)

    X = list(df['fp_str'])
    X = get_fp_np_array (X)
    
    #print (X)

    X_embedded = tsne_model.fit_transform(X)
    #print (X_embedded)

    ids = list(df['id'])

    df_embedded = pd.DataFrame ({'id': ids, 'Dim_1': X_embedded[:,0], 'Dim_2': X_embedded[:,1]})
    df_embedded = df_embedded.merge (df, on = 'id', how = 'inner')

    return (df_embedded)
    
    

def embed (df_structures, perplexity_val, fp_radius = 3, fp_length = 2048, metric = 'jaccard', random_state = 55555,  learning_rate = 200.0, n_iter = 1000, n_jobs = 3):
    """
        df_structures: Pandas dataframe, columns:
            - structure: SMILES
            - id: unique identifier of compounds
    """
    
    print ('[*] Computing t-SNE embedding using Morgan fingerprints (rad: %d, length: %d) applying metric: %s, perplexity %d, learning rate: %f, n_iteration: %d  ... ' % (fp_radius, fp_length, metric, perplexity_val, learning_rate, n_iter))

    n_dim = 2

    tsne_model = TSNE(n_components = n_dim, metric = metric, perplexity = perplexity_val, random_state = random_state, learning_rate = learning_rate, n_iter = n_iter, n_jobs = n_jobs)


    df_embedded = perform_morgan_fp_embedding (df_structures, tsne_model, fp_radius, fp_length)
    
    print (' .. done')
    
    

    return (df_embedded)
   

