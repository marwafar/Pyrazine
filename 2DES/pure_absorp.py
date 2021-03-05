import numpy as np

#-------------------------------------
# Read the file for GB and SE

R=np.genfromtxt("spectra_rephase_T0.txt", delimiter=" ")
NR=np.genfromtxt("spectra_non_rephase_T0.txt", delimiter=" ")

n_row,n_colm=R.shape
n_step=int(np.sqrt(n_row))
dt=0.5

# compute the pure absorption spectra
S_R=R[:,2]+R[:,3]*1j
S_NR=NR[:,2]+NR[:,3]*1j

pur_abs=S_R+S_NR
print(pur_abs.shape)

dw=3.335*10000/(n_step*dt)
spect=open("2DES_T0.txt","w+")
k=0
for i in range(n_step):
 for j in range(n_step):
  spect.write(str(i*dw/1000)+" "+str(j*dw/1000)+" "+str(np.real(pur_abs[k]))+" "+str(np.imag(pur_abs[k]))+" "+str(abs(pur_abs[k]))+"\n")
  k+=1
 spect.write("\n")

spect.close()
