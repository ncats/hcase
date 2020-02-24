# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

# Placeholder for creating heatmaps from aggregated embedding data.

source activate my-rdkit-env

python quantify_overlap.py

source deactivate


echo "[Step 8 done.]"
