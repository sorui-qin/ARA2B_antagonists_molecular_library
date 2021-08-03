import utils.scaffold as usc
import rdkit.Chem as rkc
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import os
import random
import argparse


def scaffold_attachment_points_number(scaffold_smi):
    '''
    Judging the number of attachment points number of scaffold molecule.
    '''
    return len(usc.get_attachment_points(scaffold_smi))


def smi_standard(smi):
    '''
    To judge the input SMILES is proper or not.
    '''
    if  rkc.MolFromSmiles(smi):
        return smi
    else:
        pass


def read_decoration_smi(textfile):
    '''
    Read the decorations SMILES from the given address.
    The decoration fragments should be saved in txt format.
    Return with a list-type of filted decorations SMILES.
    '''
    decorations=[line.strip('\n') for line in open(textfile,'r',encoding='utf-8') if smi_standard(line)]

    return decorations

    
def single_joiner(scaffold_smi, decorations):
    '''
    The function to execute the single_join fuction.
    Join only one decoration to the scaffold.
    The scaffold has only one attachment incisions.
    Scaffold and decoration are both in SMILES type.
    '''
    def join(scaffold_smi, decoration_smi):
        joined_mol=usc.join_first_attachment(scaffold_smi,decoration_smi)
        if joined_mol:
            return rkc.MolToSmiles(joined_mol)
        else:
            pass
    
    smiles=[]
    for decoration_smi in decorations:
        smile=join(scaffold_smi,decoration_smi)
        smiles.append(smile)
    return smiles


def double_joiner(scaffold_smi, decorations):
    '''
    The function to execute the double_join fuction.
    '''
        
    def double_decorations_list(decorations): #Chosing two decorations SMILES into a new list.
        dec=[]
        for st in decorations:
            for nd in decorations:
                li=[]
                li.append(st)
                li.append(nd)
                dec.append(li)
        return dec
    
    def join_in_double(scaffold_smi, decorations_list):
        '''
        The function to join the scaffold which has two attachment incisions.
        There must be put in with the SMILES of scaffold molecule,
        and two decorations SMILES stored in list-type.
        '''
        if len(decorations_list) != 2:
            return None
        smi_1 = rkc.MolToSmiles(usc.join_first_attachment(scaffold_smi, decorations_list[0]))
        smile = rkc.MolToSmiles(usc.join_first_attachment(smi_1, decorations_list[1]))
        return smile

    
    dou_dec_list=double_decorations_list(decorations)
    smiles=[]
    for dec in dou_dec_list:
        smile=join_in_double(scaffold_smi, dec)
        smiles.append(smi_standard(smile))
    return smiles


def write_results(scaffold_smi, decoration_file, file_name):
    '''
    Main fuction. Function for writing results.
    Read the scaffold in SMILES type and decorated fragments document in txt format.
    Save combined SMILES and the images of molecule in the package "file_name".
    If the input filename exists, there will create a new filename which follows a random number follow the given filename as suffix.
    '''
    def print_results(smiles, file_name):
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

    def make_dir(smiles,file_name):
        os.mkdir(file_name + "_result/")
        os.mkdir(file_name + "_result/img")
        print_results(smiles,file_name)       

    decorations=read_decoration_smi(decoration_file)
    
    att_num = scaffold_attachment_points_number(scaffold_smi)
    if att_num == 1:
        smiles=single_joiner(scaffold_smi, decorations)
        scaffold=scaffold_smi.replace('([*:0])','')
    elif att_num == 2:
        smiles=double_joiner(scaffold_smi, decorations)
        scaffold=scaffold_smi.replace('([*:0])','').replace('([*:1])','')
    else:
        return None
    
    for i in smiles:
        if not i:
            smiles.remove(i)
    
    try:
        make_dir(smiles,file_name)
    except FileExistsError:
        rand_num=str(random.randint(0,100))
        file_name=file_name +'_'+ rand_num
        print('Filename existed, results were saved in ./{}'.format(file_name))
        make_dir(smiles,file_name)

def parse_args():
    """Parses input arguments."""
    parser = argparse.ArgumentParser(description="Join a scaffold with a given decoration fragments library.")
    parser.add_argument("--scaffold-smiles", "-s", help="SMILES of sacffold.", type=str, required=True)
    parser.add_argument("--input-decorations-path", "-i",
                        help="Path to the input file with decoration fragments library in txt format.", type=str, required=True)
    parser.add_argument("--output-path", "-o",
                        help="Path to the output file.",
                        type=str, required=True)
    return parser.parse_args()

def main():
    '''Main Fuction when program run in Linux server.'''
    args = parse_args()

    write_results(args.scaffold_smiles, args.input_decorations_path, args.output_path)
