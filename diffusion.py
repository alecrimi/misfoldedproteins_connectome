###################################################################################
########## Spreading model based on Heat-kernel diffusion #########################
########## Warining: It requires: Python 3.4 due to the decorators#################
##########                                                               ##########   
########## Code initially developped by Jimmy Jackson (AIMS-Ghana) and   ##########
########## then further improved by Alessandro Crimi and Matteo Frigo    ##########
###################################################################################   

from numpy import genfromtxt
import os
import numpy as np
import matplotlib.pyplot as plt


############# Parameters
beta = 1.5  #As in the Raj et al. papers
iterations = 1000
rois = 184 #Shen atlas has 184 rois
tstar = 10.0
timestep = tstar / iterations

# Increment done with decorators
def integration_step(x0, t):
        xt = eigvecs.T @ x0
        xt = np.diag(np.exp(-beta * t * eigvals)) @ xt
        return eigvecs @ xt
 
############# Load data
# Baseline values from PET
list1 = os.listdir("b_pet")
baseline_data =  np.empty((rois,len(list1) ) ) 
for index in range(len(list1)):
#for file_name in list1:    
    baseline_data[:,index] = genfromtxt('b_pet/' + list1[index], delimiter=',')
'''
    #Store initial misfolded proteins
    X0 = np.zeros(len( control_data[:,:,index]))

    #Seed region for Pankinson is the brainstem
    #  8,13,97,182 
    X0[8] = 1
    X0[13] = 1
    X0[97] = 1
    X0[182] = 1
'''
# Structural connectomes
list2 = os.listdir("c_connectomes")
con_data =  np.empty((rois,rois,len(list2) ) )   
for index in range(len(list2)):
#for file_name in list2:      
    con_data[:,:,index] = genfromtxt('c_connectomes/' + list2[index], delimiter=',')

# compute diffusion
for index in range(len(list1)):
    print(index)
    #Store initial misfolded proteins
    X0 =  baseline_data[:,index] 
    #Compute Laplacian 
    adjacency = con_data[:,:,index]
    degrees = np.sum(adjacency, axis=1)
    dm12 = np.diag(1.0 / np.sqrt(degrees))
    laplacian = np.eye(adjacency.shape[0]) - (dm12 @ adjacency) @ dm12
    eigvals, eigvecs = np.linalg.eig(laplacian)


    #Iterate spreading    
    all_steps = [X0]  #List containing all timepoints

    for _ in range(iterations):
        next_step = integration_step(all_steps[-1], timestep)
        all_steps.append(next_step)

    A = np.asarray(all_steps)
    plt.imshow(A.T) #, interpolation='nearest'
    #plt.plot(tot)    
    plt.xlabel('Iteration' )
    plt.ylabel('ROIs' )
    plt.show() 
    #Save results
    np.savetxt(list1[index] + "_sim.csv",A.T, delimiter=" ")
