import propagation_pyrazine as propag
import numpy as np
#---------------------------------------------------------
def response(ket_coef,bra_coef):
 # Compute the response function.

 R_3_t=complex(0.0,1.0)*(np.dot(ket_coef.flatten(),np.conj(bra_coef.flatten())))

 return R_3_t
#---------------------------------------------------------
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
 # Makesure if you change the parameter you excite the required states.

 n_fock=5
 freq_cavity=4.3/27.21138386
 # g_c_12=0.24/27.21138386
 g_c_12=0.0
 g_c_01=0.24/27.21138386
 mu_ge1=1.0
 mu_ge2=1.0

 output=50
 dt_fs=0.01
 dt=dt_fs*41.34137333656
 
#--------------------------------------------------------
# compute the initial coef
 erg_t,coef_t=propag.init_wf(tot_grid,freq_tuning)
 erg_c,coef_c=propag.init_wf(tot_grid,freq_coupling)

 coef=np.outer(coef_t[:,0],coef_c[:,0])

#--------------------------------------------
 bra_coef_gs=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 bra_coef_ex1=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 bra_coef_ex2=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 ket_coef_gs=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 ket_coef_ex1=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)
 ket_coef_ex2=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 temp_coef=np.zeros((tot_grid,tot_grid,n_fock),dtype=complex)

 # Compute the initial wavefunction at time t1=0
 bra_coef_pol=np.zeros((tot_grid*tot_grid,3*n_fock),dtype=complex)
 ket_coef_pol=np.zeros((tot_grid*tot_grid,3*n_fock),dtype=complex) 

 # Transformation matrix from diabatic to polariton.
 di_to_pol=propag.polariton_pot(tot_grid,freq_tuning,freq_coupling,gradient_1,gradient_2,lambda_c,E_1,E_2,n_fock,freq_cavity,g_c_01,g_c_12)

#--------------------------------------------------------------------------------------------------------
# Run the propagation
 
 GB_file=open("response_SE_NR_01.txt", "w+")
 dipol_di10=np.zeros((n_fock,n_fock))
 for i in range(n_fock):
  dipol_di10[i,i]=0.276

 dipol_di20=np.zeros((n_fock,n_fock))
 for i in range(n_fock):
  dipol_di20[i,i]=0.979

 t2_step=5000
 
 for t in range(0,t2_step,output):

  # t1,t2,t3 for feynman diagram
  # wait is the waiting time T in fs.

  wait=0
  t1=0
  t2=dt_fs*t
  t3=wait+t2
  
  # Propagation inputs
  # total time in fs
  t_time=50
  total=t_time+t3

  n_step=int(total/dt_fs)

  for step in range(n_step):
   time_fs=dt_fs*step
   out=int(step/output)*output
  
   # 1- Compute the WF after interaction with laser at time t=0 .
   if time_fs==t1:
    # Initially we excite the fourth state (it corresponds to |2,0> state in diabatic representation).
    # !!!!!Note that this polariton state will change if we change the cavity parameters.
    # Therefore makesure to excite the required state. 

    ket_coef_ex1[:,:,0]=(coef*0.276).astype(complex) 
    bra_coef_gs[:,:,0]=coef.astype(complex)

#    print((abs(bra_coef_ex2[:,:,0])**2).sum())
#    print((abs(ket_coef_gs[:,:,0])**2).sum())
   # 2- Compute WF after interaction with laser at time=t2.
   if time_fs==t2:
    for i in range(tot_grid):
     for j in range(tot_grid):
      bra_coef_ex1[i,j,:]=np.dot(dipol_di10[:,:],bra_coef_gs[i,j,:])
      bra_coef_gs[i,j,:]=complex(0.0,0.0)

#    print((abs(bra_coef_gs[:,:,0])**2).sum())

   # 3- Compute WF after interaction with laser at time=t3.
   if time_fs==t3:
    for i in range(tot_grid):
     for j in range(tot_grid):
      ket_coef_ex1[i,j,:]=np.dot(dipol_di10[:,:],ket_coef_gs[i,j,:])        
      ket_coef_gs[i,j,:]=complex(0.0,0.0)

#    print((abs(ket_coef_ex2[:,:,0])**2).sum())

   # 4- Compute the third order response function.
   if (time_fs>=t3) and (step==out):
    for i in range(tot_grid):
     for j in range(tot_grid):
      temp_coef[i,j,:]=np.dot(dipol_di10[:,:],ket_coef_ex1[i,j,:])

    correlation=response(temp_coef[:,:,1],bra_coef_gs[:,:,1])

    GB_file.write(str(t2)+" "+str(time_fs)+" "+str(np.real(correlation))+ " "+str(np.imag(correlation))+" "+ str(abs(correlation))+ "\n")

   # Compute the bra coef at time t during the propagation.
   bra_coef_gs,bra_coef_ex1,bra_coef_ex2=propag.equation_of_motion(dt,bra_coef_gs,bra_coef_ex1,bra_coef_ex2,n_fock,g_c_01,g_c_12,\
                                        tot_grid,freq_tuning,freq_coupling,gradient_1,gradient_2,\
                                         lambda_c,E_1,E_2,freq_cavity)

   # Compute the ket coef at time t during the propagation.
   ket_coef_gs,ket_coef_ex1,ket_coef_ex2=propag.equation_of_motion(dt,ket_coef_gs,ket_coef_ex1,ket_coef_ex2,n_fock,g_c_01,g_c_12,\
                                        tot_grid,freq_tuning,freq_coupling,gradient_1,gradient_2,\
                                         lambda_c,E_1,E_2,freq_cavity)

 GB_file.close()
