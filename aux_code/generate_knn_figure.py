# Author: Gergely Zahoranszky-Kohalmi, PhD
#
# Email: gergely.zahoranszky-kohalmi@nih.gov
#
#
# Organization: National Center for Advancing Translational Sciences (NCATS/NIH)
#
#
# Ref: https://stackoverflow.com/questions/39870642/matplotlib-how-to-plot-a-high-resolution-graph
# Ref: https://www.rdkit.org/docs/GettingStartedInPython.html
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html
# Ref: https://iwatobipen.wordpress.com/2017/11/03/draw-high-quality-molecular-image-in-rdkit-rdkit/
# Ref: https://sourceforge.net/p/rdkit/mailman/message/36477772/
# Ref: https://www.rdkit.org/docs/source/rdkit.Chem.Draw.MolDrawing.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/36080986/
# Ref: https://sourceforge.net/p/rdkit/mailman/message/36311747/
# Ref: https://gist.github.com/greglandrum/d5f12058682f6b336905450e278d3399
# Ref: http://rdkit.blogspot.com/2015/02/new-drawing-code.html
# Ref: https://sourceforge.net/p/rdkit/mailman/message/36011548/



import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import DrawingOption
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.Draw import rdMolDraw2D
#from rdkit.Chem.Draw import MolDrawing
#from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold


try:
    import Image
except ImportError:
    from PIL import Image
from io import BytesIO

def DrawMolsZoomed(mols, molsPerRow=3, subImgSize=(200, 200)):
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow: nRows += 1
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    full_image = Image.new('RGBA', fullSize )
    for ii, mol in enumerate(mols):
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
        column = ii % molsPerRow
        row = ii // molsPerRow
        offset = ( column*subImgSize[0], row * subImgSize[1] )
        d2d = rdMolDraw2D.MolDraw2DCairo(subImgSize[0], subImgSize[1])
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        sub = Image.open(BytesIO(d2d.GetDrawingText()))
        full_image.paste(sub,box=offset)
    return (full_image)



def smiles2bmscaffold (smiles):
    if smiles == 'NA':
        bms = 'NA'
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            #print (Chem.MolToSmiles(mol))
            bms = MurckoScaffold.GetScaffoldForMol(mol)
            bms = Chem.MolToSmiles(bms)
    
        except:
            bms = 'NA'
    
    if len(bms) == 0:
        bms = 'NA'
    
    return (bms)


"""
def save_mol_depiction (smiles, id, dataset_name, path):
    DrawingOptions.atomLabelFontSize = 55
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 3.0
    mol = Chem.MolFromSmiles(smiles)
    fname = path + dataset_name + '_' + id + '.png'
    Chem.Draw.MolToFile(mol, fname, size=(1000, 1000))
"""


df = pd.read_csv ('../data/rnd_5_app_drugs_drugbank_knn_5.tab', sep ='\t')
print (df.head())


query_ids = []

q = list(df['knn_query_id'])
for i in range(len(q)):
    id = q[i]
    if id not in query_ids:
        query_ids.append(id)

print (query_ids)

target_ids = list(df['knn_target_id'])

molpanel = []


for i in range(len(query_ids)):
#for i in range(1,2):
    molpanel = []
    
    target_mols = []
    qid = query_ids[i]
    x = df[df['knn_query_id'] == qid].copy()
    q_structure = list(x['knn_query_structure'])[0]
    q_mol = Chem.MolFromSmiles(q_structure)
    q_mol.SetProp('_Name', qid)
    q_mol.SetProp('sim', '')
    q_mol.SetProp('data_label', '')
    molpanel.append(q_mol)

    target_ids = list(x['knn_target_id'])
    target_structures = list(x['knn_target_structure'])
    target_sim = list(x['knn_sim'])
    target_data_labels = list(x['data_label'])
    
    for j in range(len(target_ids)):
    #for j in range(1,2):
        tid = target_ids[j]
        t_structure = target_structures[j]
        tansim = ("%.2f" % target_sim[j])
        data_label = target_data_labels[j]
        #print("tansim: %s" % (tansim))
        t_mol = Chem.MolFromSmiles (t_structure)
        t_mol.SetProp('_Name', tid)
        t_mol.SetProp('sim', tansim)
        t_mol.SetProp('data_label', data_label)

        molpanel.append(t_mol) 
    fname = '../plots/knn/molecules/knn_molpanel_' + str(i + 1) + '.png'

    molsPerRow = 6
    subImgSize= (800,800)
    nRows = len(molpanel) // molsPerRow
    if len(molpanel) % molsPerRow:
        nRows += 1
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
    legends = []
    for k in range (len(molpanel)):
        mo = molpanel[k]
        l = mo.GetProp('_Name') + ' ' + mo.GetProp('sim') + ' ' + mo.GetProp('data_label')
        legends.append(l)


    d2d = rdMolDraw2D.MolDraw2DCairo(fullSize[0],fullSize[1],subImgSize[0], subImgSize[1])
    
#    d2d.drawOptions().atomLabelFontSize = 164
#    d2d.DrawingOptions.dotsPerAngstrom = 1000
#    d2d.DrawingOptions.bondLineWidth = 25.0

    highlights_atom = []

    for l in range (len(molpanel)):
        mol = molpanel[l]
        bms = smiles2bmscaffold(Chem.MolToSmiles(mol))
        highlight = mol.GetSubstructMatch(Chem.MolFromSmarts(bms))
        highlights_atom.append(highlight)
    
    d2d.drawOptions().legendFontSize=64
    d2d.DrawMolecules(molpanel,legends=legends, highlightAtoms=highlights_atom)
    d2d.FinishDrawing()
    open(fname,'wb+').write(d2d.GetDrawingText())
    #img = Draw.MolsToGridImage(molpanel,molsPerRow=6,subImgSize=(800,800),legends=[x.GetProp("_Name") for x in molpanel])
    #img = DrawMolsZoomed (molpanel, molsPerRow = 6, subImgSize=(200, 200))
    #img.save(fname)



