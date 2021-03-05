import numpy as np

#-------------------------------------
# Read the file for GB and SE

bleach_01=np.genfromtxt("response_GB_NR_01.txt", delimiter=" ")
emission_01=np.genfromtxt("response_SE_NR_01.txt", delimiter=" ")
bleach_02=np.genfromtxt("response_GB_NR_02.txt", delimiter=" ")
emission_02=np.genfromtxt("response_SE_NR_02.txt", delimiter=" ")
emission_02_01=np.genfromtxt("response_SE_NR_02_01.txt", delimiter=" ")
absorp_22_02=np.genfromtxt("response_ESA_NR_22_02.txt", delimiter=" ")
absorp_22_01=np.genfromtxt("response_ESA_NR_22_01.txt", delimiter=" ")
absorp_11_02=np.genfromtxt("response_ESA_NR_11_02.txt", delimiter=" ")
absorp_11_01=np.genfromtxt("response_ESA_NR_11_01.txt", delimiter=" ")

print(bleach_02.shape)

n_row,n_colm=bleach_02.shape
n_step=int(np.sqrt(n_row))

dephas=70
dt=0.5
# compute the response function
S_GB_01=(bleach_01[:,2]+bleach_01[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_SE_01=(emission_01[:,2]+emission_01[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_GB_02=(bleach_02[:,2]+bleach_02[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_SE_02=(emission_02[:,2]+emission_02[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_SE_02_01=(emission_02_01[:,2]+emission_02_01[:,3]*1j)*complex(0.0,1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_ESA_22_02=(absorp_22_02[:,2]+absorp_22_02[:,3]*1j)*complex(0.0,-1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_ESA_22_01=(absorp_22_01[:,2]+absorp_22_01[:,3]*1j)*complex(0.0,-1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_ESA_11_02=(absorp_11_02[:,2]+absorp_11_02[:,3]*1j)*complex(0.0,-1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))
S_ESA_11_01=(absorp_11_01[:,2]+absorp_11_01[:,3]*1j)*complex(0.0,-1.0)*np.exp(-(bleach_01[:,0]+bleach_01[:,1])/(dephas))

response=S_GB_01+S_SE_01+S_GB_02+S_SE_02+S_SE_02_01+S_ESA_22_02+S_ESA_22_01+S_ESA_11_02+S_ESA_11_01

#print(S_GB[1],bleach[1,2],bleach[1,3])
#print(S_GB.shape)

# Compute the fourier transform.
temp=np.zeros((n_step),dtype=complex)
BW=np.zeros((n_step,n_step), dtype=complex)
k=0

for i in range(n_step):
 for j in range(n_step):
  temp[j]=response[k]
  k+=1
 BW[i,:]=np.fft.ifft(temp)

FW=np.zeros((n_step,n_step), dtype=complex)

for j in range(n_step):
 for i in range(n_step):
  temp[i]=BW[i,j]
 FW[:,j]=np.fft.ifft(temp)

spect=open("spectra_non_rephase_T0.txt","w+")
for i in range(n_step):
 for j in range(n_step):
  spect.write(str(i)+" "+str(j)+" "+str(np.real(FW[i,j]))+" "+str(np.imag(FW[i,j]))+" "+str(abs(FW[i,j]))+"\n")
# spect.write("\n")

spect.close()

