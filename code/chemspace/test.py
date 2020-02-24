# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Scope: Test script for generting Scaffold Keys
#
# Reference: https://pubs.acs.org/doi/10.1021/ci5001983
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
#


import rdkit
from rdkit import Chem
from scaffold_keys import smiles2scaffoldkey, onestring

smiles = 'C1CCC(CC1)OC1=C2C=CC=CC2=CC=C1'
print (smiles)
#print(smiles2scaffoldkey(smiles))

sk = smiles2scaffoldkey(smiles, trailing_inchikey = True)
print (sk)
print (onestring (sk, has_inchikey = True))

sk = smiles2scaffoldkey(smiles, trailing_inchikey = False)
print (sk)
print (onestring (sk, has_inchikey = False))
