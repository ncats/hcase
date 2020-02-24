# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

source activate my-rdkit-env

python aggregate_embedding.py ../../data/app_drugbank_into_hc_natprod_bms_ord_5_dim_2.tab ../../plots/aggregated_hm/hm_drugbank_natprod_ord_5_dim_2.png

python aggregate_embedding.py ../../data/canvass_into_hc_natprod_bms_ord_5_dim_2.tab ../../plots/aggregated_hm/hm_canvass_natprod_ord_5_dim_2.png

source deactivate


echo "[Step 12 done.]"
