# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

echo "[*] KNN analysis started."
source activate my-rdkit-env

python analysis_knn.py

echo "[*] KNN Analysis done."


echo "[*] Not NN Analysis started."

python analysis_not_neighbors.py

echo "[*] Not NN Analysis done."

source deactivate

echo "[Step 3 done.]"


