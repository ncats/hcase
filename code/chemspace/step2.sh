# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

source activate my-rdkit-env

echo "[*] HC order: 2, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 2 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2 > ../../log/log_2_2.txt &

echo "[*] HC order: 3, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 3 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2 > ../../log/log_3_2.txt &

echo "[*] HC order: 4, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 4 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2 > ../../log/log_4_2.txt &

echo "[*] HC order: 5, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 5 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2 > ../../log/log_5_2.txt &

echo "[*] HC order: 6, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 6 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2 > ../../log/log_6_2.txt &

echo "[*] HC order: 7, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 7 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2 > ../../log/log_7_2.txt &

echo "[*] HC order: 8, dim: 2"
python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 8 2 ../../data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2 > ../../log/log_8_2.txt &

source deactivate

echo "[Step 2 done.]"
