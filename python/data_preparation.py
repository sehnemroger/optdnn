import numpy as np
import glob
import os

# Path to the data directory
path = '/mnt/e0e91479-681a-494e-8442-ea1d3ae8c6e7/Documentos/2-Mestrado/2-semestre-1/inteligencia_computacional/Trabalhos/9-Trab_29112021-TrabalhoFinal/python/data_csv/'
files = glob.glob(path+'*.txt') # Get all .txt files in the folder

for file in files:
    with open('data_csv/' + os.path.basename(file) , 'w') as f:
        data = np.genfromtxt(file)
        np.savetxt(f, data, delimiter=',')
