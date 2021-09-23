'''
This program is in order to delete the same molecules in a molecular library in txt/smi format.
Molecules should be presented in SMILES.
RDKit is needed to run this program.
'''

import rdkit
import rdkit.Chem as rkc

def standrad_smi(smile):
    return rkc.MolToSmiles(rkc.MolFromSmiles(smile))

def delete_same(li):
    fin_li=[]
    for smi in li:
        fin_li.append(standrad_smi(smi))
        del_li = sorted(set(fin_li),key=fin_li.index)
    return del_li

def write(filename):
    with open(filename,'r') as f:
        w=f.read()
        li=w.split('\n')
    fin_li=delete_same(li)
    with open(filename,'w') as fi:
        for i in fin_li:
            fi.write(i+'\n')


filename=(r'D:\Research\A2B\Pharmacophore\nonxan_All.smi')
write(filename)
