# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

# Placeholder for embed data into different spaces
# Of note, DrugBank data was embedded into ChEMBL 24.1 scsffold space in step2.sh, here CANVASS compounds will be embedded into the same space, and into ChEMBL NatProd space which will be created here.

source activate my-rdkit-env

echo "[*] HC order: 2, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 2 2 ../../data/canvass_chembl_24_1_bms_ord_2_dim_2 > ../../log/log_canvass_embedding_2_2.txt &

echo "[*] HC order: 3, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 3 2 ../../data/canvass_chembl_24_1_bms_ord_3_dim_2 > ../../log/log_canvass_embedding_3_2.txt &

echo "[*] HC order: 4, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 4 2 ../../data/canvass_chembl_24_1_bms_ord_4_dim_2 > ../../log/log_canvass_embedding_4_2.txt &

echo "[*] HC order: 5, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 5 2 ../../data/canvass_chembl_24_1_bms_ord_5_dim_2 > ../../log/log_canvass_embedding_5_2.txt &

echo "[*] HC order: 6, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 6 2 ../../data/canvass_chembl_24_1_bms_ord_6_dim_2 > ../../log/log_canvass_embedding_6_2.txt &

echo "[*] HC order: 7, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 7 2 ../../data/canvass_chembl_24_1_bms_ord_7_dim_2 > ../../log/log_canvass_embedding_7_2.txt &

echo "[*] HC order: 8, dim: 2"
python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/hc_space.tab 8 2 ../../data/canvass_chembl_24_1_bms_ord_8_dim_2 > ../../log/log_canvass_embedding_8_2.txt &

source deactivate

echo "[Step 6 done.]"
