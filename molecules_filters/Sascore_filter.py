from numpy.core.fromnumeric import sort
import SA_Score.sascorer as sa
import rdkit.Chem as rkc
import pandas as pd

def read_molecules(filename):
    with open(filename, 'r+') as f:
        li = f.read().split('\n')
    return li

def calculate(smis,n):
    def sort_di(dict_, n):
        sorted_li=sorted(dict_.items(), key = lambda kv:(kv[1], kv[0]))
        sorted_li= sorted_li[:n]
        return sorted_li

    sa_value={}
    for smi in smis:
        try:
            sa_value[smi] = sa.calculateScore(rkc.MolFromSmiles(smi))
        except ZeroDivisionError:
            pass
    return sort_di(sa_value,n)

def write(filename, sorted_li):
    def filename_check(filename):
        if '.txt' or '.smi' in filename:
            return filename.replace('.smi','').replace('.txt','')
        else:
            pass
    
    def to_smi(filename, sorted_li):
        with open('{}_SAfilted.smi'.format(filename_check(filename)),'a+') as f:
            for i in sorted_li:
                f.write(i[0] + '\n')
    
    def to_csv(filename, sorted_li):
        sorted_li = [list(i) for i in sorted_li]
        df = pd.DataFrame(sorted_li, columns=('SMILES','SA_score'))
        df.to_csv('{}_SAfilted.csv'.format(filename_check(filename)), index=None)

    to_smi(filename, sorted_li), to_csv(filename, sorted_li)
def main():
    filename = r"D:\科研\微科研\结果\R6R10_filted.smi"
    sorted_li = calculate(read_molecules(filename),5000)
    write(filename, sorted_li)

if __name__ == '__main__':
    main()