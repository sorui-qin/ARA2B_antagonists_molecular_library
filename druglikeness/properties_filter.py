import pandas
import rdkit.Chem as rkc
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import properties

def mole_proper(mol):
    num_hdonors = Lipinski.NumHDonors(mol)
    num_hacceptors = Lipinski.NumHAcceptors(mol)
    num_rotatable = Lipinski.NumRotatableBonds(mol)
    mol_weight = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)
    mol_TPSA = Descriptors.TPSA(mol)
    proper= (num_hdonors, num_hacceptors, num_rotatable, mol_weight, mol_logp, mol_TPSA)
    return proper

def filter(proper):
    (num_hdonors, num_hacceptors, num_rotatable, mol_weight, mol_logp, mol_TPSA) = proper
    if num_hdonors < 5 and num_hacceptors > 4 and mol_weight > 250 and mol_logp > 1.5 and mol_TPSA >50 and mol_TPSA < 120:
        return True
    else:
        return False

def main():
    with open(r'D:\科研\微科研\结果\腺苷受体\20210515_R6R9_L8800\lipinski_filted.txt', 'r') as f:
        smi_list=f.read()
        smi_list=smi_list.split('\n')

    Filted_smi=[]
    for smi in smi_list:
        if filter(mole_proper(rkc.MolFromSmiles(smi))):
            Filted_smi.append(smi)
        else:
            pass
    
    for smi in Filted_smi:
        with open('Molecules_properties_filted.txt','a+') as fi:
            fi.write(smi+'\n')

if __name__ == '__main__':
    main()
