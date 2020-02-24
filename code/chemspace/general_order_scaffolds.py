# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Aim: Order Bemis-Murcko scaffolds according to Ertl's intuitive scaffol ordering algorithm: Scaffold-Keys.
#
#
# Reference: https://pubs.acs.org/doi/10.1021/ci5001983
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
#
#
# Ref: https://pypi.org/project/hilbertcurve/
# Ref: https://stackoverflow.com/questions/38862293/how-to-add-incremental-numbers-to-a-new-column-using-pandas/38862389

#from hilbertcurve.hilbertcurve import HilbertCurve
import pandas as pd
from scaffold_keys import smiles2scaffoldkey, onestring, smiles2bmscaffold
import sys



def sk (smiles, trailing_inchikey = True):
    sk = smiles2scaffoldkey(smiles, trailing_inchikey)
    
    return (sk)

def sk_one (sk, has_inchikey = False):
    sk = onestring(sk, has_inchikey)
    
    return(sk)

"""
def embedhs (idx, hilbert_curve):
    coordinates = []
    coordinates = hilbert_curve.coordinates_from_distance(idx)
    coordinate_str = ''
    nr_dim = len (coordinates)

    for i in range(nr_dim):
        coordinate_str += str(coordinates[i]) + ';'
    
    coordinate_str = coordinate_str[:-1]

    return (coordinate_str)
"""

print ('[SYNTAX] python order_scaffolds.py <input> <Compound SMILES column> <ID column> <output>')

#filein = '../../data/scaffolds_chembl_24.tab'
filein = sys.argv[1]
str_col = sys.argv[2]
id_col = sys.argv[3]
fileout = sys.argv[4]



df = pd.read_csv (filein, sep = '\t')

#print (df.head())
df['bms'] = df.apply (lambda x: smiles2bmscaffold(x[str_col]), axis = 1)


nr_orig_scaffolds = df.shape[0]

#df = df.sample (100, random_state = 55555)

df['sk'] = df.apply (lambda x: sk(x['bms'], trailing_inchikey = True), axis = 1)


df = df[df['sk'] != 'NA']

df['sk_one'] = df.apply (lambda x: sk_one(x['sk'], has_inchikey = True), axis = 1)


df = df.sort_values (['sk_one'])


print ('[*] Number of scaffolds in input:')
print (nr_orig_scaffolds)


print ('[*] Number of scaffolds processed by RDKit:')
print (df.shape[0])

df = df.groupby (['sk_one'], as_index = False).agg ({
    'sk': 'first',
    id_col: 'first',
    'bms': 'first'
})


df = df.drop (columns = ['sk_one'])


df = df.rename (columns={
    'sk': 'scaffold_key',
    id_col: 'scaffold_id'
})


df = df.reset_index()
df['order'] = df.index + 1

#df['hs_coordinates'] = df.apply (lambda x: embedhs(x['order'], hilbert_curve), axis = 1)

#df = df[['structure', 'order', 'scaffold_id', 'hs_coordinates', 'scaffold_key']]
df = df[['bms', 'order', 'scaffold_id', 'scaffold_key']].copy()



df.to_csv (fileout, sep = '\t', index = False)

#print ('[*] Dimension of Scaffold-space:')
#print (N)

#print ('[*] Order of Pseudo-Hilbert-Curve:')
#print (p)



print ('[Done.]')



