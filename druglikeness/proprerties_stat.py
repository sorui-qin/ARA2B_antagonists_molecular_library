import pandas
import rdkit.Chem as rkc
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import properties


with open(r'D:\科研\大创\论文\Self\A2B antagonists\antagonists\smiles.txt', 'r') as f:
    sta_set = f.read()
    sta_set = sta_set.split(r'\n')


def mole_proper(mol):
    num_hdonors = Lipinski.NumHDonors(mol)
    num_hacceptors = Lipinski.NumHAcceptors(mol)
    num_rotatable = Lipinski.NumRotatableBonds(mol)
    num_aromatic = Lipinski.NumAromaticRings(mol)
    mol_weight = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)
    mol_TPSA = Descriptors.TPSA(mol)
    proper= [num_hdonors, num_hacceptors, num_rotatable, num_aromatic, mol_weight, mol_logp, mol_TPSA]
    return proper

def save_excel(proper_list):
    df = pandas.DataFrame(proper_list, columns=['HBD', 'HBA', 'Rotatable bones', 'Mol Weight', 'LogP', 'TPSA'])
    df.to_excel('Molecules Properties.xlsx', index=None)

def main():
    proper_list=[]
    for smi in sta_set:
        mol = rkc.MolFromSmiles(smi)
        proper = mole_proper(mol)
        proper_list.append(proper)

    save_excel(proper_list)

if __name__ == '__main__':
    main()
