import random

def write(filename, n):
    def read(filename,n):
        with open(filename, 'r+') as f:
            li = f.read().split('\n')
            random_list = random.sample(li, n)
        return random_list

    def filename_check(filename):
        if '.txt' or '.smi' in filename:
            return filename.replace('.smi','').replace('.txt','')
        else:
            pass
    
    def to_smi(filename, li, n):
        with open('{0}_random_{1}.smi'.format(filename_check(filename), n),'a+') as f:
            for i in li:
                f.write(i + '\n')
    to_smi(filename, read(filename,n), n)

if __name__ == '__main__':
    filename= r"D:\Research\A2B\molecules\R8R10_filted.smi"
    write(filename, 2000)