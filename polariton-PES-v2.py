import numpy as np
from numpy import linalg as LA
#--------------------------------------
def diag(M):
 #Diagonalize matrix
 E,V=LA.eigh(M)
 return E,V
#--------------------------------------
def x_Hermit(tot_grid):
 # Compute the grid points using Hermite polynomial

 # Here we employ mass*freq=1.0

 x_eq=0.0
 x_ij=np.zeros((tot_grid,tot_grid))
 for row in range(tot_grid):
  for colm in range(tot_grid):
   x_ij[row,colm]=np.sqrt((row+1.0)/(2.0))*float(row==colm-1)\
                  +x_eq*float(row==colm)+np.sqrt(row/2.0)*float(row==colm+1)

 x_i,vect=diag(x_ij)
 return x_i,vect
#----------------------------
def molecular_diabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2):

 # Compute the diabatic Hamiltonian

 x_i,vect=x_Hermit(tot_grid)

 H_t=0.5*freq_t*x_i*x_i
 H_c=0.5*freq_c*x_i*x_i
 H_0= H_t[:,None]+H_c

 H_t=E_1+0.5*freq_t*x_i*x_i+grad_1*x_i
 H_1=H_t[:,None]+H_c

 H_t=E_2+0.5*freq_t*x_i*x_i+grad_2*x_i
 H_2=H_t[:,None]+H_c

 H_12=lamda*x_i

 return H_0,H_1,H_2,H_12
#----------------------------------------------------------------------------------------------
def molecular_adiabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2):

 # Compute the adiabatic potential and the transformation matrix from diabatic to adiabatic.

 H_0,H_1,H_2,H_12=molecular_diabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2)
 x_i,vect=x_Hermit(tot_grid)

 ad_pot=open("ad_pot.csv","w+")
 di_pot=open("di_pot.txt","w+")
 cpl=open("diabatic_coupling.txt","w+")

 n_state=3
 ED=np.zeros((n_state,n_state))
 H_ad=np.zeros((tot_grid*tot_grid,n_state,n_state))
# di_to_ad=np.zeros((tot_grid*tot_grid,n_state,n_state))
 k=0
 for i in range(tot_grid):
  for j in range(tot_grid):
   ED[0,0]=H_0[i,j]
   ED[1,1]=H_1[i,j]
   ED[2,2]=H_2[i,j]
   ED[1,2]=H_12[j]
   ED[2,1]=H_12[j]
   E_ad,U=diag(ED)
#   di_to_ad[k,:,:]=U
   H_ad[k,0,0]=E_ad[0]
   H_ad[k,1,1]=E_ad[1]
   H_ad[k,2,2]=E_ad[2]

   k+=1

   ad_pot.write(str(x_i[i])+","+str(x_i[j])+","+str(E_ad[0])+","+ str(E_ad[1])+","+str(E_ad[2])+"\n")
   di_pot.write(str(x_i[i])+" "+str(x_i[j])+" "+str(H_0[i,j])+" "+ str(H_1[i,j])+" "+str(H_2[i,j])+"\n")
   cpl.write(str(x_i[i])+" "+str(x_i[j])+" "+str(H_12[j])+"\n")
#  ad_pot.write("\n")
  di_pot.write("\n")
  cpl.write("\n")

 ad_pot.close()
 di_pot.close()
 cpl.close()

 return H_ad
#----------------------------------------------------------------------------------------
def polariton(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity,g_c_01,g_c_12,g_c_02):

 # Compute the polariton PES.
 H_0,H_1,H_2,H_12=molecular_diabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2)
 x_i,vect=x_Hermit(tot_grid) 
 H_ad=molecular_adiabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2)

 n_state=3
 ED=np.zeros((n_state,n_state))
 Hij=np.zeros((3*n_fock,3*n_fock))
 
 PES=open("Polariton-PES-g01-240.csv","w+")
 PES_di=open("di-dress-PES-g01-240.csv","w+")
 
 for x in range(tot_grid):
  for y in range(tot_grid):
   ED[0,0]=H_0[x,y]
   ED[1,1]=H_1[x,y]
   ED[2,2]=H_2[x,y]
   ED[1,2]=H_12[y]
   ED[2,1]=H_12[y]

   for row in range(3*n_fock):
    a=int(row/n_fock)
    m=row%n_fock

    for colm in range(3*n_fock):
     b=int(colm/n_fock)
     n=colm%n_fock
     
     Hij[row,colm]=ED[a,b]*float(m==n)+ (n+0.5)*freq_cavity*float(a==b)*float(m==n)
     Hij[row,colm]+=g_c_01*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))*\
	float((a==1 and b==0) or (a==0 and b==1))
     Hij[row,colm]+=g_c_12*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))*\
        float((a==1 and b==2) or (a==2 and b==1))
     Hij[row,colm]+=g_c_02*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))*\
        float((a==0 and b==2) or (a==2 and b==0))
#     Hij[row,colm]+=g_c_01*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))*float(a!=b)

   E_pol,U=diag(Hij)
   PES.write(str(x_i[x])+ "," + str(x_i[y]) + "," + ",".join(E_pol.astype(str))+"\n")
   PES_di.write(str(x_i[x])+ "," + str(x_i[y]) + "," +str(Hij[0,0])+","+str(Hij[1,1])+","+str(Hij[4,4])+","+str(Hij[8,8])+"\n")
#  PES.write("\n")
#  PES_di.write("\n")
 PES.close()
 return
#---------------------------
if __name__== "__main__":

 # The parameters of the molecular Hamiltonian
 tot_grid=21
 cm_to_h=219474.63
 freq_t=597.0/cm_to_h
 freq_c=952.0/cm_to_h
 grad_1=-847.0/cm_to_h
 grad_2=1202.0/cm_to_h
 lamda=2110.0/cm_to_h
 E_1=31800/cm_to_h
 E_2=39000/cm_to_h

 # The parameters of the cavity
 n_fock=4
 freq_cavity=4.3/27.21138386
# g_c_12=0.24/27.21138386
 g_c_12=0.0
 g_c_01=0.24/27.21138386
 g_c_02=0.0
 mu_ge1=1.0
 mu_ge2=1.0
#-------------------------------
 # Compute the PES
 H_ad=molecular_adiabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2)
 polariton(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity,g_c_01,g_c_12,g_c_02)
