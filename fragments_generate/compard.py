import argparse
import rdkit.Chem as rkc
import re

def standard_smi(smile):
    if rkc.MolFromSmiles(smile):
        return rkc.MolToSmiles(rkc.MolFromSmiles(smile))

def delete_attachment(smile):
    atta=rkc.MolFromSmiles('*')
    mol=rkc.MolFromSmiles(smile)
    rm=rkc.DeleteSubstructs(mol,atta)
    return rkc.MolToSmiles(rm)

def compare(att_txt,frag_txt,savefile):
    def read_txt(filename):
        with open(filename,'r') as fi:
            li=fi.read()
            li=re.split(r'[;\n\t]',li)
            li=list(set(li))
        return li

    def compared_smi(attached_smilist,frag_smilist):
        target_smi=[]
        for smile in attached_smilist:
            if delete_attachment(smile) in frag_smilist:
                target_smi.append(smile)
        return target_smi
    
    def saved_file(filename,target_smi):
        with open('{}.txt'.format(filename),'a+') as f:
            for smile in target_smi:
                f.write(smile+'\n')
    
    att=read_txt(att_txt)
    attached_smilist=[]
    for smi in att:
        if "[*]"  in smi:
            attached_smilist.append(smi)
    
    frag_smili = read_txt(frag_txt)
    frag_smilist = [standard_smi(smi) for smi in frag_smili]
    target_smi=compared_smi(attached_smilist,frag_smilist)
    saved_file(savefile,target_smi)

def parse_args():
    """Parses input arguments."""
    parser = argparse.ArgumentParser(description="Compared the attached SMILES with an given fragments library.")
    parser.add_argument("--input-decorations-path", "-i",
                        help="Path to the input file with sliced fragments library in txt format.", type=str, required=True)
    parser.add_argument("--compared-decorations-path", "-c",
                        help="Path to the input file with specific fragments library in txt format.", type=str, required=True)
    parser.add_argument("--output-path-name", "-o",
                        help="Name of the output file.",
                        type=str, required=True)
    return parser.parse_args()

def main():
    arg = parse_args()
    compare(arg.input_decorations_path,arg.compared_decorations_path,arg.output_path_name)
    
if __name__ == '__main__':
    main()