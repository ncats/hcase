# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

source activate my-rdkit-env

echo "[*] Cherry-Picking ChEMBL 24_1 Scaffolds ..."

python cherrypick_scaffolds.py

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 2, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 2 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_2_dim_2 > ../../log/log_sc_2_2.txt &

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 3, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 3 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_3_dim_2 > ../../log/log_sc_3_2.txt &

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 4, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 4 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_4_dim_2 > ../../log/log_sc_4_2.txt &

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 5, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 5 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_5_dim_2 > ../../log/log_sc_5_2.txt &

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 6, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 6 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_6_dim_2 > ../../log/log_sc_6_2.txt &

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 7, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 7 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_7_dim_2 > ../../log/log_sc_7_2.txt &

echo "[*] Cherry-Picked Scaffolds Analysis.  HC order: 8, dim: 2"
python embed_structure.py ../../data/cherrypicked_scaffolds.tab structure scaffold_id ../../data/hc_space.tab 8 2 ../../data/cp_scaffolds_chembl_24_1_bms_ord_8_dim_2 > ../../log/log_sc_8_2.txt &

source deactivate

echo "[Step 4 done.]"


