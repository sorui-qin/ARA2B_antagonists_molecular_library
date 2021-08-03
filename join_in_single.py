import utils.scaffold as usc
import rdkit.Chem as rkc
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import os
import random

def join(scaffold_smi, decoration_smi):
    joined_mol=usc.join_first_attachment(scaffold_smi,decoration_smi)
    if joined_mol:
        return rkc.MolToSmiles(joined_mol)
    else:
        pass

def smi_standard(smi):
    if  rkc.MolFromSmiles(smi):
        return smi
    else:
        return None

def read_decoration_smi(textfile):
    decorations=[line[:-1] for line in open(textfile) if smi_standard(line)]
    return decorations  

def write_results(scaffold_smi, decoration_file, file_name):
    '''
    Function for writing results.
    Read the scaffold in smiles type and decorated fragments document in txt format.
    Save combined smiles and molecule in the package "file_name".
    '''
    def print_results(smiles, sacffold, file_name):
        with open(file_name + "_result/smiles.txt", "wt") as f:
            all_smi_num=0
            for i, smile in enumerate(smiles):
                f.write(smile + "\n")
                #mol = rkc.MolFromSmiles(smile)
                #scf = rkc.MolFromSmiles(scaffold)
                #sub = mol.GetSubstructMatches(scf)
                #Draw.MolToFile(mol, file_name + "_result/img/" + str(i) + ".png", size=(600, 600), highlightAtoms=sub[0] , legend=smile)
                all_smi_num += 1
        print('{0} results saved in {1}'.format(all_smi_num, file_name))

    def make_dir(smiles,scaffold,file_name):
        os.mkdir(file_name + "_result/")
        os.mkdir(file_name + "_result/img")
        print_results(smiles,scaffold,file_name)       

    decorations=read_decoration_smi(decoration_file)
    smiles=[]
    for decoration_smi in decorations:
        smiles.append(join(scaffold_smi,decoration_smi))
    for i in smiles:
        if not i:
            smiles.remove(i)

    scaffold=scaffold_smi.replace('([*:0])','')

    try:
        make_dir(smiles,scaffold,file_name)
        print_results(smiles,scaffold,file_name)
    except FileExistsError:
        rand_num=str(random.randint(0,100))
        file_name=file_name +'_'+ rand_num
        print('Filename existed, results were saved in ./{}'.format(file_name))
        make_dir(smiles,scaffold,file_name)


    

