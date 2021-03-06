import numpy as np
from numpy import linalg as LA
import math
#---------------------------------
def diag(M):

 # Diagonalize matrix
 E,V=LA.eigh(M)
 return E,V

#-----------------------
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

#---------------------------
def weight_Hermit(tot_grid):
 # Compute the weight for each grid point.
 
 x_eq=0.0

 x_i,vect=x_Hermit(tot_grid)
 w_i= ((1.0/np.pi)**(-0.25)*np.exp(0.5*x_i*x_i)*vect[0,:])**2

 return w_i
#-------------------------------------------
def second_derivavtive_Hermit(tot_grid):

 # Compute the second derivavtive  

 x_i,vect=x_Hermit(tot_grid)

 k=np.array([i+0.5 for i in range(tot_grid)])

 dif2mat=-2.0*np.matmul(np.matmul(vect.T,np.diag(k)),vect)
 dif2mat+=(np.diag(x_i))**2

 return dif2mat
#-------------------------------------------------------
def init_wf(tot_grid,freq):
 
 x_i,vect=x_Hermit(tot_grid)
 dif2mat=second_derivavtive_Hermit(tot_grid)

 pot=0.5*freq*x_i*x_i
 H_ij=-0.5*freq*dif2mat+np.diag(pot)
 energy,coef=diag(H_ij)

 return energy,coef

#----------------------------------------------------------------------------------
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

 ad_pot=open("ad_pot.txt","w+")
 di_pot=open("di_pot.txt","w+")
 cpl=open("diabatic_coupling.txt","w+")

 n_state=3
 ED=np.zeros((n_state,n_state))
 di_to_ad=np.zeros((tot_grid*tot_grid,n_state,n_state))
 k=0
 for i in range(tot_grid):
  for j in range(tot_grid):
   ED[0,0]=H_0[i,j]
   ED[1,1]=H_1[i,j]
   ED[2,2]=H_2[i,j]
   ED[1,2]=H_12[j]
   ED[2,1]=H_12[j]
   E_ad,U=diag(ED)
   di_to_ad[k,:,:]=U
   k+=1

   ad_pot.write(str(x_i[i])+" "+str(x_i[j])+" "+str(E_ad[0])+" "+ str(E_ad[1])+" "+str(E_ad[2])+"\n")
   di_pot.write(str(x_i[i])+" "+str(x_i[j])+" "+str(H_0[i,j])+" "+ str(H_1[i,j])+" "+str(H_2[i,j])+"\n")
   cpl.write(str(x_i[i])+" "+str(x_i[j])+" "+str(H_12[j])+"\n")
  ad_pot.write("\n")
  di_pot.write("\n")
  cpl.write("\n")

 ad_pot.close() 
 di_pot.close()
 cpl.close()   

 return di_to_ad
#----------------------------------------------------------------------------------------
def kinetic_energy_operator(tot_grid,coef_gs,coef_ex1,coef_ex2,freq_t,freq_c,n_fock):

 # compute the kinetic energy operator for the propagation.
 dif2mat=second_derivavtive_Hermit(tot_grid)
 
 KEO_gs=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 KEO_ex1=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 KEO_ex2=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 for m in range(n_fock):
  KEO_gs[:,:,m]=-0.5*freq_t*(dif2mat.dot(coef_gs[:,:,m]))-0.5*freq_c*(coef_gs[:,:,m].dot(dif2mat.T))
  KEO_ex1[:,:,m]=-0.5*freq_t*(dif2mat.dot(coef_ex1[:,:,m]))-0.5*freq_c*(coef_ex1[:,:,m].dot(dif2mat.T))
  KEO_ex2[:,:,m]=-0.5*freq_t*(dif2mat.dot(coef_ex2[:,:,m]))-0.5*freq_c*(coef_ex2[:,:,m].dot(dif2mat.T))


 return KEO_gs,KEO_ex1,KEO_ex2
#---------------------------------------------------------------------------------------
def potential_energy_operator(tot_grid,coef_gs,coef_ex1,coef_ex2,freq_t,freq_c,\
                             grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity):

 # compute the potential energy operator for the propagation.
 H_0,H_1,H_2,H_12=molecular_diabatic_pot(tot_grid,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2) 

 PEO_gs=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 PEO_ex1=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 PEO_ex2=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 diab_cpl21=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 diab_cpl12=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 for m in range(n_fock):
  PEO_gs[:,:,m]=coef_gs[:,:,m]*(H_0+(m+0.5)*freq_cavity)
  PEO_ex1[:,:,m]=coef_ex1[:,:,m]*(H_1+(m+0.5)*freq_cavity)
  PEO_ex2[:,:,m]=coef_ex2[:,:,m]*(H_2+(m+0.5)*freq_cavity)
  diab_cpl21[:,:,m]=H_12[:]*coef_ex1[:,:,m]
  diab_cpl12[:,:,m]=H_12[:]*coef_ex2[:,:,m]    
 
 return PEO_gs,PEO_ex1,PEO_ex2,diab_cpl12,diab_cpl21
#---------------------------------------------------------------------
def interaction(coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12):

 Hcm_01=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 Hcm_10=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 Hcm_12=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 Hcm_21=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 # Compute the interaction operator between the molecule and the cavity.  
 for n in range(n_fock):
  for m in range(n_fock):
   Hcm_01[:,:,n]+=g_c_01*coef_ex1[:,:,m]*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))
   Hcm_10[:,:,n]+=g_c_01*coef_gs[:,:,m]*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n)) 
   Hcm_12[:,:,n]+=g_c_12*coef_ex2[:,:,m]*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))
   Hcm_21[:,:,n]+=g_c_12*coef_ex1[:,:,m]*(np.sqrt(m)*float(m-1==n)+np.sqrt(m+1)*float(m+1==n))
 
 return Hcm_01,Hcm_10,Hcm_12,Hcm_21
#-------------------------------------------------------
def equation_of_motion(dt,coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12,tot_grid,freq_t,freq_c,\
                       grad_1,grad_2,lamda,E_1,E_2,freq_cavity):

 # Compute the numerical integration using RK4 method.
 
 old_coef_gs=np.copy(coef_gs)
 old_coef_ex1=np.copy(coef_ex1)
 old_coef_ex2=np.copy(coef_ex2)

 # RK1
 KEO_gs,KEO_ex1,KEO_ex2=kinetic_energy_operator(tot_grid,coef_gs,coef_ex1,coef_ex2,freq_t,freq_c,n_fock)
 PEO_gs,PEO_ex1,PEO_ex2,diab_cpl12,diab_cpl21=potential_energy_operator(tot_grid,coef_gs,coef_ex1,\
                                   coef_ex2,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity)
 Hcm_01,Hcm_10,Hcm_12,Hcm_21=interaction(coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12)
 
 RK1_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK1_ex1=complex(0.0,-1.0)*(KEO_ex1+PEO_ex1+diab_cpl12+Hcm_10+Hcm_12)
 RK1_ex2=complex(0.0,-1.0)*(KEO_ex2+PEO_ex2+diab_cpl21+Hcm_21)

 coef_gs=old_coef_gs+dt/2.0*RK1_gs
 coef_ex1=old_coef_ex1+dt/2.0*RK1_ex1
 coef_ex2=old_coef_ex2+dt/2.0*RK1_ex2

 #RK2
 KEO_gs,KEO_ex1,KEO_ex2=kinetic_energy_operator(tot_grid,coef_gs,coef_ex1,coef_ex2,freq_t,freq_c,n_fock)
 PEO_gs,PEO_ex1,PEO_ex2,diab_cpl12,diab_cpl21=potential_energy_operator(tot_grid,coef_gs,coef_ex1,\
                                   coef_ex2,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity)
 Hcm_01,Hcm_10,Hcm_12,Hcm_21=interaction(coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12)

 RK2_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK2_ex1=complex(0.0,-1.0)*(KEO_ex1+PEO_ex1+diab_cpl12+Hcm_10+Hcm_12)
 RK2_ex2=complex(0.0,-1.0)*(KEO_ex2+PEO_ex2+diab_cpl21+Hcm_21)

 coef_gs=old_coef_gs+dt/2.0*RK2_gs
 coef_ex1=old_coef_ex1+dt/2.0*RK2_ex1
 coef_ex2=old_coef_ex2+dt/2.0*RK2_ex2

 #RK3
 KEO_gs,KEO_ex1,KEO_ex2=kinetic_energy_operator(tot_grid,coef_gs,coef_ex1,coef_ex2,freq_t,freq_c,n_fock)
 PEO_gs,PEO_ex1,PEO_ex2,diab_cpl12,diab_cpl21=potential_energy_operator(tot_grid,coef_gs,coef_ex1,\
                                   coef_ex2,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity)
 Hcm_01,Hcm_10,Hcm_12,Hcm_21=interaction(coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12)

 RK3_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK3_ex1=complex(0.0,-1.0)*(KEO_ex1+PEO_ex1+diab_cpl12+Hcm_10+Hcm_12)
 RK3_ex2=complex(0.0,-1.0)*(KEO_ex2+PEO_ex2+diab_cpl21+Hcm_21)

 coef_gs=old_coef_gs+dt*RK3_gs
 coef_ex1=old_coef_ex1+dt*RK3_ex1
 coef_ex2=old_coef_ex2+dt*RK3_ex2

 #RK4
 KEO_gs,KEO_ex1,KEO_ex2=kinetic_energy_operator(tot_grid,coef_gs,coef_ex1,coef_ex2,freq_t,freq_c,n_fock)
 PEO_gs,PEO_ex1,PEO_ex2,diab_cpl12,diab_cpl21=potential_energy_operator(tot_grid,coef_gs,coef_ex1,\
                                   coef_ex2,freq_t,freq_c,grad_1,grad_2,lamda,E_1,E_2,n_fock,freq_cavity)
 Hcm_01,Hcm_10,Hcm_12,Hcm_21=interaction(coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12)

 RK4_gs=complex(0.0,-1.0)*(KEO_gs+PEO_gs+Hcm_01)
 RK4_ex1=complex(0.0,-1.0)*(KEO_ex1+PEO_ex1+diab_cpl12+Hcm_10+Hcm_12)
 RK4_ex2=complex(0.0,-1.0)*(KEO_ex2+PEO_ex2+diab_cpl21+Hcm_21)

 coef_gs=old_coef_gs+dt/6.0*(RK1_gs+2.0*RK2_gs+2.0*RK3_gs+RK4_gs)
 coef_ex1=old_coef_ex1+dt/6.0*(RK1_ex1+2.0*RK2_ex1+2.0*RK3_ex1+RK4_ex1)
 coef_ex2=old_coef_ex2+dt/6.0*(RK1_ex2+2.0*RK2_ex2+2.0*RK3_ex2+RK4_ex2)

 return coef_gs,coef_ex1,coef_ex2

#-----------------------------------------------------------
if __name__== "__main__":

 # The parameters of the molecular Hamiltonian
 tot_grid=21
 cm_to_h=219474.63
 freq_tuning=597.0/cm_to_h
 freq_coupling=952.0/cm_to_h
 gradient_1=-847.0/cm_to_h
 gradient_2=1202.0/cm_to_h
 lambda_c=2110.0/cm_to_h
 E_1=31800/cm_to_h
 E_2=39000/cm_to_h

 # The parameters of the cavity
 n_fock=3
 freq_cavity=4.3/27.21138386
# g_c_12=0.24/27.21138386
 g_c_12=0.0
 g_c_01=0.0
 mu_ge1=1.0
 mu_ge2=1.0

 # Propagation inputs
 n_steps=10
 output=2
 dt=0.01*41.34137333656
#--------------------------------------------------------
# compute the initial coef
 erg_t,coef_t=init_wf(tot_grid,freq_tuning)
 erg_c,coef_c=init_wf(tot_grid,freq_coupling)

 coef=np.outer(coef_t[:,0],coef_c[:,0])
 # check normalization
# norm=(coef**2).sum()
# print(norm)
#-------- 
 coef_gs=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 coef_ex1=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 coef_ex2=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 coef_ex2[:,:,0]=coef.astype(complex)

 # transform the initial wavefunction.
# di_to_ad=molecular_adiabatic_pot(tot_grid,freq_tuning,freq_coupling,gradient_1,gradient_2,lambda_c,E_1,E_2)
 
# coef_gs[:,:,0]=(coef_gs[:,:,0].flatten()*di_to_ad[:,0,0]+coef_ex1[:,:,0].flatten()\
#                      *di_to_ad[:,0,1]+coef_ex2[:,:,0].flatten()*di_to_ad[:,0,2]).reshape(tot_grid,tot_grid) 
# coef_ex1[:,:,0]=(coef_gs[:,:,0].flatten()*di_to_ad[:,1,0]+coef_ex1[:,:,0].flatten()\
#                      *di_to_ad[:,1,1]+coef_ex2[:,:,0].flatten()*di_to_ad[:,1,2]).reshape(tot_grid,tot_grid)
# coef_ex2[:,:,0]=(coef_gs[:,:,0].flatten()*di_to_ad[:,2,0]+coef_ex1[:,:,0].flatten()\
#                      *di_to_ad[:,2,1]+coef_ex2[:,:,0].flatten()*di_to_ad[:,2,2]).reshape(tot_grid,tot_grid)

#----------------------------------------------------------
 # Run the propagation.
 
 di_to_ad=molecular_adiabatic_pot(tot_grid,freq_tuning,freq_coupling,gradient_1,gradient_2,lambda_c,E_1,E_2)

 # Initialize arrays:
 pop_gs=np.zeros((n_fock),dtype=complex)
 pop_ex1=np.zeros((n_fock),dtype=complex)
 pop_ex2=np.zeros((n_fock),dtype=complex)

 coef_gs_ad=np.zeros((tot_grid*tot_grid,n_fock),dtype=complex)
 coef_ex1_ad=np.zeros((tot_grid*tot_grid,n_fock),dtype=complex)
 coef_ex2_ad=np.zeros((tot_grid*tot_grid,n_fock),dtype=complex)

 # file name
 pop_ad_gs=open("ad_pop_gs.txt","w+")
 pop_ad_ex1=open("ad_pop_ex1.txt","w+")
 pop_ad_ex2=open("ad_pop_ex2.txt","w+")
 pop_di_gs=open("di_pop_gs.txt","w+")
 pop_di_ex1=open("di_pop_ex1.txt","w+")
 pop_di_ex2=open("di_pop_ex2.txt","w+")

 for step in range(n_steps):

  time_fs=dt*step/41.34137333656
  out=int(step/output)*output

  # compute the diabatic population and export it.
  if (step==out):
   for n in range(n_fock):
    pop_gs[n]=np.matmul(coef_gs[:,:,n].flatten(),np.conj(coef_gs[:,:,n].flatten()))
    pop_ex1[n]=np.matmul(coef_ex1[:,:,n].flatten(),np.conj(coef_ex1[:,:,n].flatten()))     
    pop_ex2[n]=np.matmul(coef_ex2[:,:,n].flatten(),np.conj(coef_ex2[:,:,n].flatten()))
   pop_di_gs.write(str(time_fs)+" "+ " ".join(np.real(pop_gs).astype(str))+"\n")
   pop_di_ex1.write(str(time_fs)+" "+ " ".join(np.real(pop_ex1).astype(str))+"\n")
   pop_di_ex2.write(str(time_fs)+" "+ " ".join(np.real(pop_ex2).astype(str))+"\n")

   #compute the adiabatic population and export it.
   for n in range(n_fock):
    coef_gs_ad[:,n]=coef_gs[:,:,n].flatten()*di_to_ad[:,0,0]+coef_ex1[:,:,n].flatten()\
                      *di_to_ad[:,1,0]+coef_ex2[:,:,n].flatten()*di_to_ad[:,2,0]
    coef_ex1_ad[:,n]=coef_gs[:,:,n].flatten()*di_to_ad[:,0,1]+coef_ex1[:,:,n].flatten()\
                      *di_to_ad[:,1,1]+coef_ex2[:,:,n].flatten()*di_to_ad[:,2,1]
    coef_ex2_ad[:,n]=coef_gs[:,:,n].flatten()*di_to_ad[:,0,2]+coef_ex1[:,:,n].flatten()\
                      *di_to_ad[:,1,2]+coef_ex2[:,:,n].flatten()*di_to_ad[:,2,2]
    
   for n in range(n_fock):
    pop_gs[n]=np.matmul(coef_gs_ad[:,n],np.conj(coef_gs_ad[:,n]))
    pop_ex1[n]=np.matmul(coef_ex1_ad[:,n],np.conj(coef_ex1_ad[:,n]))     
    pop_ex2[n]=np.matmul(coef_ex2_ad[:,n],np.conj(coef_ex2_ad[:,n]))
  
   pop_ad_gs.write(str(time_fs)+" "+ " ".join(np.real(pop_gs).astype(str))+"\n")
   pop_ad_ex1.write(str(time_fs)+" "+ " ".join(np.real(pop_ex1).astype(str))+"\n")
   pop_ad_ex2.write(str(time_fs)+" "+ " ".join(np.real(pop_ex2).astype(str))+"\n")
    
  # Compute the coef at time t
  coef_gs,coef_ex1,coef_ex2=equation_of_motion(dt,coef_gs,coef_ex1,coef_ex2,n_fock,g_c_01,g_c_12,\
                                        tot_grid,freq_tuning,freq_coupling,gradient_1,gradient_2,\
                                         lambda_c,E_1,E_2,freq_cavity)
 
 pop_ad_gs.close()
 pop_ad_ex1.close()
 pop_ad_ex2.close()
 pop_di_gs.close()
 pop_di_ex1.close()
 pop_di_ex2.close()

