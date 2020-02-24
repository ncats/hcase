# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

source activate my-rdkit-env

python order_scaffolds.py ../../data/scaffolds_chembl_24.tab ../../data/hc_space.tab

source deactivate

echo "[Step 1 done.]"
