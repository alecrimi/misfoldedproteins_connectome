from numpy import genfromtxt
import os
import numpy as np

#Parameters
beta = 1.5
iterations = 10
rois = 184
 
#Load data
list = os.listdir("HC")
control_data =  np.empty((rois,rois,len(list) ) ) 
for index in range(len(list)):
#for file_name in list:    
    control_data[:,:,index] = genfromtxt('HC/' + list[index], delimiter=',')

 
list2 = os.listdir("PD")
pd_data =  np.empty((rois,rois,len(list2) ) )   
for index in range(len(list2)):
#for file_name in list:      
    pd_data[:,:,index] = genfromtxt('PD/' + list2[index], delimiter=',')

# compute diffusion
for index in range(len(list)):
    print(index)
    #Compute Laplacian 
    D = np.diag(sum(control_data[:,:,index]))
    # This is not normalized
    L = D - control_data[:,:,index]

    #Store initial misfolded proteins
    X0 = np.zeros(len( control_data[:,:,index]))

    #Seed region for Pankinson is the brainstem
    #  8,13,97,182 
    X0[8] = 1
    X0[13] = 1
    X0[97] = 1
    X0[182] = 1

    #Compute first spreading
    increment =  np.exp(-L*beta)
    Xt = increment.dot(X0)

    #Iterate spreading    
    for iter in range(iterations):
 
       Xt = increment.dot(Xt) + Xt ;
 
     
    #Save results
    np.savetxt(list[index] + "sim.csv", Xt, delimiter=" ")
