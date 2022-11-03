import glob
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib

path = os.path.abspath(os.getcwd())
files = glob.glob(path+'/data_csv/'+'*.txt') # Get all .txt files in the folder

def substring(string,start,end):
    lstart = len(start)
    idx_start = string.find(start)+lstart
    idx_end = string.find(end)
    return string[idx_start:idx_end]

# criando vetor para acomodar os valores de cis do dataset
cis = np.zeros((len(files),4))

for i,file in enumerate(files):
    cis[i,0] = float(substring(files[i],'x1i_','x2i_'))
    cis[i,1] = float(substring(files[i],'x2i_','x3i_'))
    cis[i,2] = float(substring(files[i],'x3i_','x4i_'))
    cis[i,3] = float(substring(files[i],'x4i_','.txt'))

matplotlib.rcParams['text.usetex'] = True # habilitando tex 

fig, axs = plt.subplots(4,3,figsize=(6,4), tight_layout=True)

ii, jj = 0, 0
for i in range(4):
    for j in range(4):
        if i != j:
            print(i,j,ii,jj)
            axs[ii,jj].plot(cis[:,i],cis[:,j],'.b')
            axs[ii,jj].set_xlabel('$x_'+str(i+1)+'$')
            axs[ii,jj].set_ylabel('$x_'+str(j+1)+'$')
            jj += 1
    ii += 1
    jj = 0

plt.show()