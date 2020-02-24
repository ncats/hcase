# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#

source activate my-rdkit-env

echo "[*] t-SNE perplexity 5 started..."
python altered_tsne_analysis.py 5 ../../data/tsne_per_5.tab ../../plots/tsne/tsne_5.png &

echo "[*] t-SNE perplexity 10 started..."
python altered_tsne_analysis.py 10 ../../data/tsne_per_10.tab ../../plots/tsne/tsne_10.png &

echo "[*] t-SNE perplexity 20 started..."
python altered_tsne_analysis.py 20 ../../data/tsne_per_20.tab ../../plots/tsne/tsne_20.png &

echo "[*] t-SNE perplexity 30 started..."
python altered_tsne_analysis.py 30 ../../data/tsne_per_30.tab ../../plots/tsne/tsne_30.png &

echo "[*] t-SNE perplexity 40 started..."
python altered_tsne_analysis.py 40 ../../data/tsne_per_40.tab ../../plots/tsne/tsne_40.png &

echo "[*] t-SNE perplexity 50 started..."
python altered_tsne_analysis.py 50 ../../data/tsne_per_50.tab ../../plots/tsne/tsne_50.png &

source deactivate

echo "[Step 10 done.]"
