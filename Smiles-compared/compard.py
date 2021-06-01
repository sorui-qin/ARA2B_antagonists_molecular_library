import rdkit
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
    
    frag_smilist=read_txt(frag_txt)
    target_smi=compared_smi(attached_smilist,frag_smilist)
    saved_file(savefile,target_smi)

compare(r'D:\科研\微科研\smiles比对\drd2.recap.txt',r'D:\科研\微科研\片段库\L8800_standrad.txt','drd2_frag')


