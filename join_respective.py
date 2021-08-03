import argparse
import os
import rdkit.Chem as rkc
import utils.scaffold as usc

def read_fragments_library(textfile):
    with open(textfile, 'r+') as frli:
        frag_li = frli.read() 
        frag_li = frag_li.split('\n')
    return frag_li

def join(scaffold_smi, decoration_smi):
    joined_mol=usc.join_first_attachment(scaffold_smi,decoration_smi)
    if joined_mol:
        return rkc.MolToSmiles(joined_mol)
    else:
        pass

def join_respective(textfile1, textfile2, scaffold_smi):
    def sca_accuracy(scaffold_smi):
        if '([*:0])' and '([*:1])' in scaffold_smi:
            return True
        else:
            return False
    
    frag_li1 = read_fragments_library(textfile1) 
    frag_li2 = read_fragments_library(textfile2)
    if sca_accuracy(scaffold_smi):
        first_join_list=[]
        for first_decoration_smi in frag_li1:
            first_join = join(scaffold_smi, first_decoration_smi)
            if first_join:
                first_join = first_join.replace('[*:1]','[*:0]')
                first_join_list.append(first_join)
            else:
                print('SMILES format error.')
                break
        
        final_join_list=[]
        for decorated_smi in first_join_list:
            for decoration_smi in frag_li2:
                fin_join = join(decorated_smi, decoration_smi)
                if fin_join:
                    final_join_list.append(fin_join)
                else:
                    print('SMILES format error.')
                    break
    return final_join_list

def write_in_txt(li,file_name):
    os.mkdir(file_name)
    with open(file_name + './join.txt', 'a+') as fin:
        for smi in li:
            fin.write(smi+'\n')
    
def parse_args():
    """Parses input arguments."""
    parser = argparse.ArgumentParser(description="Join a scaffold with two decoration fragments libraries respectively.")
    parser.add_argument("--scaffold-smiles", "-s", help="SMILES of sacffold.", type=str, required=True)
    parser.add_argument("--first-decorations-path", "--st",
                        help="Path to the input file of the first decoration fragments library in txt format.", type=str, required=True)
    parser.add_argument("--second-decorations-path", "--nd",
                        help="Path to the input file of the second decoration fragments library in txt format.", type=str, required=True)
    parser.add_argument("--filename", "-f",
                        help="The name of the final output file.",
                        type=str, required=True)
    return parser.parse_args()

def main():
    args = parse_args()
    final_join_list = join_respective(args.first_decorations_path, args.second_decorations_path, args.scaffold_smiles)
    write_in_txt(final_join_list, args.filename)

if __name__ == "__main__":
    final_join_list = join_respective(r'D:\Anaconda\envs\my-rdkit-env\utils\1.txt',r'D:\Anaconda\envs\my-rdkit-env\utils\2.txt','NC1=NC(CC2=C([*:0])C=CC=C2([*:1]))=CN3C1=NC4=C3C=CC=C4')
    write_in_txt(final_join_list, 'test')
    
