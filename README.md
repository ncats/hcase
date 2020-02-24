#### LICENSE Related Section #######

* Links to licenses


- MIT License: [ https://opensource.org/licenses/MIT ]

- Creative Commons Attribution-ShareAlike 3.0 Unported License: [https://creativecommons.org/licenses/by-sa/3.0/]

- Creative Commons Attribution 4.0 International Public License (CC-BY 4.0 International): [ https://creativecommons.org/licenses/by/4.0/legalcode.txt ]

- GPL-2 License: [ https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html ]

- GPL-3 License: [ https://www.gnu.org/licenses/gpl-3.0.en.html ]

- Apache License 2.0: [ https://www.apache.org/licenses/LICENSE-2.0 ]

- 3-Clause BSD license: [https://opensource.org/licenses/BSD-3-Clause, https://github.com/scikit-learn/scikit-learn/blob/master/COPYING]



* Licenses of resources:

- ChEMBL DB version 24:	[ https://chembl.gitbook.io/chembl-interface-documentation/about ]	Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)	[ https://creativecommons.org/licenses/by-sa/3.0/ ]


- ChEMBL Natural Products Subset Data:	Creative Commons Attribution-ShareAlike 3.0 Unported License	[ https://creativecommons.org/licenses/by-sa/3.0/ ]


- DrugBank 5.0.9:	Creative Common’s Attribution-NonCommercial 4.0 International License	[ https://www.drugbank.ca/releases/5-0-9 , https://creativecommons.org/licenses/by-nc/4.0/legalcode ]


- CANVASS:	Creative Commons Attribution 4.0 International Public License	CC-BY 4.0 International	[ https://creativecommons.org/licenses/by/4.0/legalcode.txt ]


- SciKit-Learn Python module:	[ http://scikit-learn.org, https://github.com/scikit-learn/scikit-learn ]	3-Clause BSD license	[ https://github.com/scikit-learn/scikit-learn/blob/master/COPYING ]


- SciPy library:	[ https://github.com/scipy/scipy ]	BSD 3-Clause "New" or "Revised" License	[ https://github.com/scipy/scipy/blob/master/LICENSE.txt ] 


- RDKit	[ https://www.rdkit.org/ ]	BSD license	[ https://github.com/rdkit/rdkit/blob/master/license.txt ]


- hilbertcurve.hilbertcurve:	[ https://github.com/galtay/hilbertcurve ]	MIT License	[ https://github.com/galtay/hilbertcurve/blob/develop/LICENSE ]


- MatPlotLib:	[ https://github.com/matplotlib/matplotlib ]	Custom	[ https://github.com/matplotlib/matplotlib/blob/master/LICENSE/LICENSE ]


- NumPy:	Custom	[ https://numpy.org/license.html ]


- Pandas:	[ https://github.com/pandas-dev/pandas]	BSD 3-Clause License	[ https://github.com/pandas-dev/pandas/blob/master/LICENSE ] 


- Seaborn:	[ https://github.com/mwaskom/seaborn ]	BSD 3-Clause "New" or "Revised" License	[ https://github.com/mwaskom/seaborn/blob/master/LICENSE ] 


Derivative work of ChEMBL Natural Products dataset that is distributed under Creative Commons Attribution-ShareAlike 3.0 Unported License.




#### Workflow to replicate the analysis of the manuscript ######



* Step 1: Order reference Bemis-Murcko scaffold set according to Ert'l Scaffol-Key algorithm.

> bash step1.sh


* Step 2: Map structures onto a Hilbert-curve created from the ordered reference scaffold set.

> bash step2.sh

-- Make sure that the log files are exited with '[Done.]' in log/log_x_2.txt files where x = 2-8. 

* Step 3: HCASE Embedding of (k=5) KNNs of 5 Randomly Selected Molecules (Reference Scaffolds: ChEMBL, compounds: DrugBank)

> bash step3.sh

* Step 4: HCASE Embedding of Cherry-Picked Scaffolds and their 50+50 Immediate Neighbors Based on Scaffold-Key Ordering

> bash step4.sh

--Make sure that the log files are exited with '[Done.]' in log/sc_log_x_2.txt files where x = 2-8.


* Step 5: Analysis of HCASE Embedding of 25 Randomly Selected Molecules Randomly Selected Molecules (Making Plots)

> bash step5.sh


* Step 6: Embedding CANVASS library into ChEMBL 24 Scaffold Space.

> bash step6.sh

--Make sure that the log files are exited with '[Done.]' in  log/log_canvass_embedding_x_2.txt files where x = 2-8.

* Step 7: Embedding DrugBank into NatProd Scaffold Space.

> bash step7.sh

--Make sure that the log files are exited with '[Done.]' in log/app_drugbank_into_hc_natprod_log_x_2.txt files where x = 2-5.

* Step 8: Quantify Overlap.

> bash step8.sh

* Step 9: Chebyshev Distance and Rank Distance Correlation.

> bash step9.sh


* Step 10: Altered t-SNE Analysis of Cherry-Picked Scaffolds.

> bash step10.sh

* Step 11: Altered t-SNE Embedding of DrugBank Molecules into ChEMBL Space.

> bash step11.sh

* Step 12: Quantifying Overlap of HCASE-Embeddings of DrugBank/CANVASS Compounds in NatProd Space

> bash step12.sh

* Step 13: Original t-SNE Analysis of DrugBank Molecules.

> bash step13.sh



* Creating Hilbert-Curve Mapped Scaffold Embedding (HMSE)

[SYNTAX] python order_scaffolds.py <input> <output>
python order_scaffolds.py ../../data/scaffolds_chembl_24.tab ../../data/hc_space.tab

* Embedding structures into HMSE

[SYNTAX] python embed_structure.py <tab-separated file of structures> <Column name in structures file containing SMILES> <Column name of structure IDs> <Hilbert-Curve Mapped Scaffold Embeding>

python embed_structure.py ../../data/STD_drugbank_approved_structures_v5.txt Structure ID ../../data/hc_space.tab 8 2



* Annotation and license related remarks of input and generated datasets:

data/app_drugbank_all_embedding_coords_natprod_bms.tab	HCASE embedding of DrugBank approved drugs in NatProd scaffold space at every PHC-order value.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
data/app_drugbank_into_hc_natprod_bms_ord_2_dim_2.tab	HCASE embedding of DrugBank approved drugs in NatProd scaffold space using PHC-2.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
data/app_drugbank_into_hc_natprod_bms_ord_3_dim_2.tab	HCASE embedding of DrugBank approved drugs in NatProd scaffold space using PHC-3.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
data/app_drugbank_into_hc_natprod_bms_ord_4_dim_2.tab	HCASE embedding of DrugBank approved drugs in NatProd scaffold space using PHC-4.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
data/app_drugbank_into_hc_natprod_bms_ord_5_dim_2.tab	HCASE embedding of DrugBank approved drugs in NatProd scaffold space using PHC-5.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
data/app_drugs_drugbank_all_knn_5.tab	Similarity matrix of 5 randomly selected query molecules from DrugBank dataset others molecules DrugBank dataset, sorted by decreasing order of similarity.	Derivative work of DrugBank dataset.
data/app_drugs_drugbank_all_not_nn_5.tab	25 randomly selected molecules from DrugBank approved drugs dataset.	Derivative work of DrugBank dataset.
data/app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-2.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-3.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-4.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-5.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-6.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-7.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.tab	HCASE embedding of DrugBank approved drugs in ChEMBL (24.1) scaffold space using PHC-8.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_all_embedding_coords_chembl_24_1_bms.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space at every PHC-order value.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_all_embedding_coords_natprod_bms.tab	HCASE embedding of CANVASS compounds in NatProd scaffold space at every PHC-order value.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/canvass_chembl_24_1_bms_ord_2_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-2.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_chembl_24_1_bms_ord_3_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-3.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_chembl_24_1_bms_ord_4_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-4.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_chembl_24_1_bms_ord_5_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-5.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_chembl_24_1_bms_ord_6_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-6.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_chembl_24_1_bms_ord_7_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-7.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_chembl_24_1_bms_ord_8_dim_2.tab	HCASE embedding of CANVASS compounds in ChEMBL (24.1) scaffold space using PHC-8.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/canvass_into_hc_natprod_bms_ord_2_dim_2.tab	HCASE embedding of CANVASS compounds in NatProd scaffold space using PHC-2.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/canvass_into_hc_natprod_bms_ord_3_dim_2.tab	HCASE embedding of CANVASS compounds in NatProd scaffold space using PHC-3.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/canvass_into_hc_natprod_bms_ord_4_dim_2.tab	HCASE embedding of CANVASS compounds in NatProd scaffold space using PHC-4.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/canvass_into_hc_natprod_bms_ord_5_dim_2.tab	HCASE embedding of CANVASS compounds in NatProd scaffold space using PHC-5.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/chebyshev_stat_full.tab	Correlation statistics of rank-order and Chebyshev-distances of all embeddings (both compound libraries in both scaffolds spaces).	CC-BY 4.0 International
data/chebyshev_stat.tab	Correlation statistics of rank-order and Chebyshev-distances of DrugBank dataset embedding in the two scaffold spaces.	CC-BY 4.0 International
data/cherrypicked_scaffolds.tab	Cherry-picked scaffolds from ChEMBL 24.1 dataset.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_all_embedding_coords_chembl_24_1_bms.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds at every PHC order value.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_2_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-2.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_3_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-3.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_4_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-4.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_5_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-5.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_6_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-6.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_7_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-7.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/cp_scaffolds_chembl_24_1_bms_ord_8_dim_2.tab	HCASE embedding of cherry-picked ChEMBL 24.1 scaffolds using PHC-8.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/hc_space.tab	HCASE space defined by the coordinates of ChEMBL 24.1 scaffolds.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/knn_coords_app_drugs_drugbank_chembl_24_1_bms.tab	HCASE embedding of DrugBank molecules  in ChEMBL 24.1 scaffolds space at every PHC order value, highlighting the 5 randomly selected query molecules and their 5 nearest neighbors.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/lo_orig_tsne_perplexity_10.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=10. Coordinates of only the 5 randomly selected molecules and their 5 nearest neughbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/lo_orig_tsne_perplexity_20.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=20. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/lo_orig_tsne_perplexity_30.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=30. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/lo_orig_tsne_perplexity_35.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=35 (not used in analysis). Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/lo_orig_tsne_perplexity_40.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=40. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/lo_orig_tsne_perplexity_50.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=50. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/lo_orig_tsne_perplexity_5.tab	Original t-SNE embedding of a 90% subset of DrugBank molecules (merged with some other highlighted molecules) at perplexity=5. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/natprod_hc_space.tab	HCASE space defined by the coordinates of NatProd scaffolds.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/orig_tsne_perplexity_10.tab	Original t-SNE embedding of DrugBank molecules at perplexity=10. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/orig_tsne_perplexity_20.tab	Original t-SNE embedding of DrugBank molecules at perplexity=20. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/orig_tsne_perplexity_30.tab	Original t-SNE embedding of DrugBank molecules at perplexity=30. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/orig_tsne_perplexity_35.tab	Original t-SNE embedding of DrugBank molecules at perplexity=35 (not used in analysis). Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/orig_tsne_perplexity_40.tab	Original t-SNE embedding of DrugBank molecules at perplexity=40. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/orig_tsne_perplexity_50.tab	Original t-SNE embedding of DrugBank molecules at perplexity=50. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/orig_tsne_perplexity_5.tab	Original t-SNE embedding of DrugBank molecules at perplexity=5. Coordinates of only the 5 randomly selected molecules and their 5 nearest neighbors from the DrugBank dataset.	Derivative work of DrugBank dataset.
data/PUB_KNIME_PrepChEMBL_Workflow_Generated_Scaffold_Table.tab.csv	Input CheMBL 24.1 Bemis-Murcko Scaffolds.	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/quantified_overlaps.tab	Quantified overlap results of the embedding of CANVASS and DrugBank datasets in NatProd scaffolds spaces.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
data/rnd_5_app_drugs_drugbank_knn_5.tab	Randomly selected 5 query molecules and their 5 nearest neighbors from DrugBank approved drugs dataset.	Derivative work of DrugBank dataset.
data/scaffolds_chembl_24.tab	Symlink to data/PUB_KNIME_PrepChEMBL_Workflow_Generated_Scaffold_Table.tab.csv	Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/similarity_app_drugs_drugbank_all_not_nn_5.tab	Similarity matrix of 25 randomly selected query molecules from DrugBank dataset vs. their 5 nearest neighborsin DrugBank dataset, sorted by decreasing order of similarity.	Derivative work of DrugBank dataset.
data/STD_drugbank_approved_structures_v5.txt	Input datafile of DrugBank approved drugs.	Derivative work of DrugBank dataset.
data/STD_FINAL_20170525_ML_CANVASS_deduplicated.txt	Input datafile of CANVASS molecules.	CC-BY 4.0 International
data/STD_ML_ChEMBL23_NatProd_10132017.txt	Input datafile of ChEMBL 23 NatProd molecules.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
data/tsne_per_10.tab	Embedding of DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=10. Coordinates only of th 5 randomly selected molecules and their 5 nearest neighbors  from DrugBank dataset.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/tsne_per_20.tab	Embedding of DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=20. Coordinates only of th 5 randomly selected molecules and their 5 nearest neighbors  from DrugBank dataset.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/tsne_per_30.tab	Embedding of DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=30. Coordinates only of th 5 randomly selected molecules and their 5 nearest neighbors  from DrugBank dataset.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/tsne_per_40.tab	Embedding of DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=40. Coordinates only of th 5 randomly selected molecules and their 5 nearest neighbors  from DrugBank dataset.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/tsne_per_50.tab	Embedding of DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=50. Coordinates only of th 5 randomly selected molecules and their 5 nearest neighbors  from DrugBank dataset.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/tsne_per_5.tab	Embedding of DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=5. Coordinates only of th 5 randomly selected molecules and their 5 nearest neighbors  from DrugBank dataset.	Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
data/tsne_embedding_coords_all_cmpd_lo_orig_tsne_perplexity_10.tab	Embedding of 90% subset (merged with other highlighted molecules) DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=10.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_lo_orig_tsne_perplexity_20.tab	Embedding of 90% subset (merged with other highlighted molecules) DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=20.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_lo_orig_tsne_perplexity_30.tab	Embedding of 90% subset (merged with other highlighted molecules) DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=30.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_lo_orig_tsne_perplexity_40.tab	Embedding of 90% subset (merged with other highlighted molecules) DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=40.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_lo_orig_tsne_perplexity_50.tab	Embedding of 90% subset (merged with other highlighted molecules) DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=50.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_lo_orig_tsne_perplexity_5.tab	Embedding of 90% subset (merged with other highlighted molecules) DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=5.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_orig_tsne_perplexity_10.tab	Embedding of all DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=10.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_orig_tsne_perplexity_20.tab	Embedding of all DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=20.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_orig_tsne_perplexity_30.tab	Embedding of all DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=30.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_orig_tsne_perplexity_40.tab	Embedding of all DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=40.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_orig_tsne_perplexity_50.tab	Embedding of all DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=50.	Derivative work of DrugBank dataset.
data/tsne_embedding_coords_all_cmpd_orig_tsne_perplexity_5.tab	Embedding of all DrugBank molecules with the help of Scaffold-tSNE method in ChEMBL 24.1 scaffold space at perplexity=5.	Derivative work of DrugBank dataset.
		
plots/aggregated_hm/hm_canvass_natprod_ord_5_dim_2.png	Heatmap: number of molecules associated with a specific coordinate in the HCASE embedded space aggregated, cells colored accordingly. Embedding: CANVASS molecules in NatProd scaffolds spaces.	Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/aggregated_hm/hm_drugbank_natprod_ord_5_dim_2.png	Heatmap: number of molecules associated with a specific coordinate in the HCASE embedded space aggregated, cells colored accordingly. Embedding: DrugBank molecules in NatProd scaffolds spaces.	Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
		
		
plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_2_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_3_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_4_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/app_drugbank_into_hc_natprod_bms_ord_5_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_2_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_3_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_4_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_5_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_6_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_7_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_chembl_24_1_bms_ord_8_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_2_dim_2.png		Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_3_dim_2.png		Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_4_dim_2.png		Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/canvass_into_hc_natprod_bms_ord_5_dim_2.png		Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_2_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_3_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_4_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_5_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_6_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_7_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_chembl_bms_ord_8_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_natrpod_ord_2_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_natrpod_ord_3_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_natrpod_ord_4_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/comparative_embeddings/dual_canvass_drugbank_natrpod_ord_5_dim_2.png		Derivative work of DrugBank and ChEMBL Natural Products (ChEMBL v23) datasets.
		
		
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_2_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_3_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_4_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_5_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_6_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_7_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/cp_scaffolds/cp_scaffolds_chembl_24_1_bms_ord_8_dim_2.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
		
		
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_2_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_3_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_4_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_5_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_6_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_7_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/knn_5_app_drugs_drugbank_chembl_24_1_bms_ord_8_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_2_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_3_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_4_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_5_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_6_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_7_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/knn/not_nn_app_drugbank_chembl_24_1_bms_ord_8_dim_2.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
		
		
plots/knn/molecules/knn_molpanel_1.png		Derivative work of DrugBank dataset.
plots/knn/molecules/knn_molpanel_2.png		Derivative work of DrugBank dataset.
plots/knn/molecules/knn_molpanel_3.png		Derivative work of DrugBank dataset.
plots/knn/molecules/knn_molpanel_4.png		Derivative work of DrugBank dataset.
plots/knn/molecules/knn_molpanel_5.png		Derivative work of DrugBank dataset.
plots/knn/molecules/not_nn_molpanel.png		Derivative work of DrugBank dataset.
		
		
plots/tsne/lo_orig_tsne_perplexity_10.png		Derivative work of DrugBank dataset.
plots/tsne/lo_orig_tsne_perplexity_20.png		Derivative work of DrugBank dataset.
plots/tsne/lo_orig_tsne_perplexity_30.png		Derivative work of DrugBank dataset.
plots/tsne/lo_orig_tsne_perplexity_40.png		Derivative work of DrugBank dataset.
plots/tsne/lo_orig_tsne_perplexity_50.png		Derivative work of DrugBank dataset.
plots/tsne/lo_orig_tsne_perplexity_5.png		Derivative work of DrugBank dataset.
plots/tsne/orig_tsne_perplexity_10.png		Derivative work of DrugBank dataset.
plots/tsne/orig_tsne_perplexity_20.png		Derivative work of DrugBank dataset.
plots/tsne/orig_tsne_perplexity_30.png		Derivative work of DrugBank dataset.
plots/tsne/orig_tsne_perplexity_40.png		Derivative work of DrugBank dataset.
plots/tsne/orig_tsne_perplexity_50.png		Derivative work of DrugBank dataset.
plots/tsne/orig_tsne_perplexity_5.png		Derivative work of DrugBank dataset.
plots/tsne/tsne_10.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_20.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_30.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_40.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_50.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_5.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_knn_perplexity_10.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_knn_perplexity_20.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_knn_perplexity_30.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_knn_perplexity_40.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_knn_perplexity_50.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/tsne/tsne_knn_perplexity_5.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
		
		
plots/main_text_edited_figures/Figure1.png		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/Figure1.xcf		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/Figure2.png		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/Figure2.xcf		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/Figure3.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure3.xcf		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure4.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure4.xcf		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure5.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure5.xcf		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure6.png		Derivative work of DrugBank, ChEMBL scaffolds (ChEMBL v24_1) and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/main_text_edited_figures/Figure6.xcf		Derivative work of DrugBank, ChEMBL scaffolds (ChEMBL v24_1) and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/main_text_edited_figures/Figure7.png		Derivative work of DrugBank, ChEMBL scaffolds (ChEMBL v24_1) and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/main_text_edited_figures/Figure7.xcf		Derivative work of DrugBank, ChEMBL scaffolds (ChEMBL v24_1) and ChEMBL Natural Products (ChEMBL v23) datasets.
plots/main_text_edited_figures/Figure8.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/Figure8.xcf		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS1.png		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/FigureS1.xcf		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/FigureS2.png		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/FigureS2.xcf		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/FigureS3.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS3.xcf		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS4.png		Derivative work of DrugBank dataset.
plots/main_text_edited_figures/FigureS5.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS5.xcf		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS6.png		Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/main_text_edited_figures/FigureS6.xcf		Derivative work of ChEMBL Natural Products (ChEMBL v23) datasets.
plots/main_text_edited_figures/FigureS7.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS7.xcf		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS8.png		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS8.xcf		Derivative work of ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS9.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/FigureS9.xcf		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/GraphicalAbstract.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/main_text_edited_figures/GraphicalAbstract.xcf		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
		
		
plots/Algorithm_a.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/Algorithm_b.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.
plots/Algorithm_c.png		Derivative work of DrugBank and ChEMBL scaffolds (ChEMBL v24_1) datasets.



