# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

source activate my-rdkit-env

python chebyshev_stat.py

python chebyshev_stat_full.py 

source deactivate

echo "[Step 9 done.]"
