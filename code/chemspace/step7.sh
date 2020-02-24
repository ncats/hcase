# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

# Generate NatProd HC-Space and embed Aprroved Drugs from DrugBank and CANVASS compounds into this space.
source activate my-rdkit-env

echo "[*] Generating NatPRod HC Space..."


python general_order_scaffolds.py ../../data/STD_ML_ChEMBL23_NatProd_10132017.txt Structure ID ../../data/natprod_hc_space.tab

echo "[Done.]"

### CANVASS into NatProd HC-Space


echo "[*] CANVASS into NatProd HC order: 2, dim: 2"

python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/natprod_hc_space.tab 2 2 ../../data/canvass_into_hc_natprod_bms_ord_2_dim_2 > ../../log/canvass_into_hc_natprod_log_2_2.txt &

echo "[*] CANVASS into NatProd HC order: 3, dim: 2"

python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/natprod_hc_space.tab 3 2 ../../data/canvass_into_hc_natprod_bms_ord_3_dim_2 > ../../log/canvass_into_hc_natprod_log_3_2.txt &

echo "[*] CANVASS into NatProd HC order: 4, dim: 2"

python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/natprod_hc_space.tab 4 2 ../../data/canvass_into_hc_natprod_bms_ord_4_dim_2 > ../../log/canvass_into_hc_natprod_log_4_2.txt &

echo "[*] CANVASS into NatProd HC order: 5, dim: 2"

python embed_structure.py ../../data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt Structure ID ../../data/natprod_hc_space.tab 5 2 ../../data/canvass_into_hc_natprod_bms_ord_5_dim_2 > ../../log/canvass_into_hc_natprod_log_5_2.txt &


### DrugBank into NatProd HC-Space

echo "[*] app_drugbank into NatProd HC order: 2, dim: 2"

python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/natprod_hc_space.tab 2 2 ../../data/app_drugbank_into_hc_natprod_bms_ord_2_dim_2 > ../../log/app_drugbank_into_hc_natprod_log_2_2.txt &

echo "[*] app_drugbank into NatProd HC order: 3, dim: 2"

python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/natprod_hc_space.tab 3 2 ../../data/app_drugbank_into_hc_natprod_bms_ord_3_dim_2 > ../../log/app_drugbank_into_hc_natprod_log_3_2.txt &

echo "[*] app_drugbank into NatProd HC order: 4, dim: 2"

python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/natprod_hc_space.tab 4 2 ../../data/app_drugbank_into_hc_natprod_bms_ord_4_dim_2 > ../../log/app_drugbank_into_hc_natprod_log_4_2.txt &

echo "[*] app_drugbank into NatProd HC order: 5, dim: 2"

python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/natprod_hc_space.tab 5 2 ../../data/app_drugbank_into_hc_natprod_bms_ord_5_dim_2 > ../../log/app_drugbank_into_hc_natprod_log_5_2.txt &


#### Analyzing Above Embeddings

echo "[*] Analyzing embeddings into NatProd HC Space..."

python analysis_embeddings.py

source deactivate

echo "[Step 7 done.]"



