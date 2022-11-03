import numpy as np
import glob
import os

# Path to the data directory
path = '/mnt/e0e91479-681a-494e-8442-ea1d3ae8c6e7/Documentos/2-Mestrado/2-semestre-1/inteligencia_computacional/Trabalhos/9-Trab_29112021-TrabalhoFinal/python/data_csv/'
files = glob.glob(path+'*.txt') # Get all .txt files in the folder

file = files[3]
print(file)
data = np.loadtxt(file,delimiter=',')

print(data[0,0:4])

# Indice de desempenho
r = 1
q = 1
Q = q*np.eye(4)
def desemp(x,u,t):
    J_em_t = x.dot(Q).dot(x)+u**2+q*t
    return J_em_t

t=data[:,8]
Jt = np.zeros_like(t)
for i in range(len(t)):
    Jt[i] = desemp(sol.y[:,i],F[i],t[i])


J = np.trapz(Jt,x=t)
print(J)
