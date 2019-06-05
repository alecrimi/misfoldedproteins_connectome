import numpy as np
import os
import math as mt
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import zscore
#from sklearn.preprocessing import normalize

# this function compute the correlation matrix
def static_connectivity_matrix(subject):
  corr= np.corrcoef(subject)
  return corr


def Simulation(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate):
	
 # A function to simulate the spread of misfolded beta_amyloid

 ##input parameters (inside parenthesis are values used in the paper)	
 #N_regions: number of regions 
 #v: speed (1)
 # dt: time step (0.01)
 # T_total: total time steps (10000)
 #GBA: GBA gene expression (zscore, N_regions * 1 vector) (empirical GBA expression)
 #SNCA: SNCA gene expression after normalization (zscore, N_regions * 1 vector) (empirical SNCA expression)
 # sconnLen: structural connectivity matrix (length) (estimated from HCP data)
 # sconnDen: structural connectivity matrix (strength) (estimated from HCP data)
 # ROIsize: region sizes (voxel counts)
 # seed: seed region of misfolded beta-amyloid injection (choose as you like? )
 # syn_control: a parameter to control the number of voxels in which beta-amyloid may get synthesized (region size, i.e., ROIsize)
 # init_number: number of injected misfolded beta-amyloid (1) 
 # prob_stay: the probability of staying in the same region per unit time (0.5)
 # trans_rate: a scalar value, controlling the baseline infectivity

 ## output parameters
 # Rnor_all: A N_regions * T_total matrix, recording the number of normal beta-amyloid in regions
 # Rmis_all: A N_regions * T_total matrix, recording the number of misfolded beta-amyloid in regions
 # Rnor0: a N_Regions * 1 vector, the population of normal agents in regions before pathogenic spreading 
 # Pnor0: a N_Regions * 1 vecotr, the population of normal agents in edges before pathogenic spreading 
 #Pnor_all: a N_regions * N_regions * T_total matrix, recording the number of normal beta_amyloid in paths could be memory-consuming
 #Pmis_all: a N_regions * N_regions * T_total matrix, recording the number of misfolded beta_amyloid in paths could be memoryconsuming
 

 
 #make sure the diag is zero
 sconnDen = sconnDen - np.diag(np.diag(sconnDen))
 sconnLen = sconnLen - np.diag(np.diag(sconnLen))
 
 #set the mobility pattern
 weights = sconnDen
 delta0 = 1 * trans_rate /ROIsize 
 g = 0.5 #global tuning variable that  quantifies the temporal MP deposition inequality among the different brain regions
 mu_noise = 0.2 #mean of the additive noise
 sigma_noise = 0.1 # standard deviation of the additive noise
 Ki = np.random.normal(mu_noise, sigma_noise, (N_regions,N_regions))

 #model
 #regional probability receiving MP infectous-like agents

 Epsilon = np.zeros((N_regions, N_regions))
 for i in range(N_regions):
  t = 0
  for j in range(N_regions):
   if i != j:
    t = t +  (sconnDen[i][j] * (g*(1 - np.exp(-delta0 *prob_stay))) * prob_stay + sconnDen[i][i] * (1 - g)*(1 - np.exp( - delta0*prob_stay))* prob_stay)
    Epsilon[i][j] = t
 weights = (1 - prob_stay)* Epsilon - prob_stay * np.exp(- init_number* prob_stay) +  Ki



 #The probability of moving from region i to edge (i,j)
 sum_line= [sum(weights[i]) for i in range(len(weights))]
 Total_el_col= np.tile(np.transpose(sum_line), (1,1))
 weights = weights / Total_el_col 
 
 #convert gene expression scores to probabilities
 clearance_rate = norm.cdf(zscore(GBA)) 
 synthesis_rate = norm.cdf(zscore(SNCA))

 #store the number of normal/misfolded beta-amyloid at each time step
 Rnor_all = np.zeros((N_regions, T_total))
 Rmis_all = np.zeros((N_regions, T_total))
 Pnor_all = np.zeros((N_regions, T_total))
 Pmis_all = np.zeros((N_regions, T_total))

 #Rnor, Rmis, Pnor, Pmis store results of single simulation at each time
 Rnor = np.zeros(N_regions)# number of normal beta-amyloid in regions
 Rmis = np.zeros(N_regions) #number of misfolded beta-amyloid in regions
 Pnor = np.zeros((N_regions, N_regions)) # number of normal beta-amyloid in paths
 Pmis= np.zeros((N_regions, N_regions)) # number of misfolded beta-amyloid in paths

 ##normal alpha-syn growth 
 # fill the network with normal proteins

 iter_max = 10000
 
 #normal alpha synuclein growth
 for t in range(iter_max):  
  ###moving process
  # regions towards paths
  # movDrt stores the number of proteins towards each region. i.e. moving towards l
  #movDrt = np.kron(np.ones((1, N_regions)), Rnor) * weights
  movDrt = Rnor * weights
  movDrt = movDrt * dt 
  movDrt = movDrt - np.diag(np.diag(movDrt))
    
    
  # paths towards regions
  # update moving
	  

  movOut = Pnor * v / sconnLen  #longer path & smaller v = lower probability of moving out of paths
  movOut = movOut - np.diag(np.diag(movOut))

  Pnor = Pnor - movOut * dt + movDrt
  Pnor = Pnor -  np.diag(np.diag(Pnor))
  Sum_rows_movOut = [sum(movOut[i])for i in range(len(movOut))]
  Sum_cols_movDrt = [sum(movDrt[:,i]) for i in range(len(movDrt))]

  Rtmp = Rnor
  Rnor = Rnor +  np.transpose(Sum_rows_movOut) * dt -  Sum_cols_movDrt

  #growth process	
  Rnor = Rnor -  Rnor * (1 - np.exp(- clearance_rate * dt))   + (synthesis_rate * syn_control) * dt

  if np.absolute(Rnor - Rtmp) < (1e-7 * Rtmp):
   break
	
 
 Pnor0 = Pnor
 Rnor0 = Rnor

 # misfolded protein spreading process
 
 #inject misfolded beat_amyloid
 Rmis[seed] = init_number;
 
 for t  in range (T_total):
  #moving process
  # normal proteins: region -->> paths
  #movDrt_nor = np.kron(np.ones((1, N_regions)), Rnor) * weights * dt
  movDrt_nor = Rnor * weights * dt
  movDrt_nor = movDrt_nor - np.diag(np.diag(movDrt_nor))     

  movOut_nor = Pnor * v  / sconnLen 
  movOut_nor = movOut_nor - np.diag(np.diag(movOut_nor))
    
    
  #misfolded proteins: region -->> paths
  movDrt_mis =  Rnor * weights * dt
  #movDrt_mis = np.kron(np.ones((1, N_regions)), Rnor) * weights * dt
  movDrt_mis = movDrt_mis - np.diag(np.diag(movDrt_mis))
     
  #misfolded proteins: paths -->> regions
  movOut_mis = Pmis * v / sconnLen
  movOut_mis = movOut_mis - np.diag(np.diag(movOut_mis))
    
  #update regions and paths
  Pnor = Pnor - movOut_nor * dt + movDrt_nor 
  Pnor = Pnor - np.diag(np.diag(Pnor))

  Sum_rows_movOut_nor = [sum(movOut_nor[i]) for i in range(len(movOut_nor))]
  Sum_cols_movDrt_nor = [sum(movDrt_nor[:,i]) for i in range(len(movDrt_nor))]

  Rnor = Rnor + np.transpose(Sum_rows_movOut_nor ) * dt -  Sum_cols_movDrt_nor

  Pmis = Pmis - movOut_mis*dt + movDrt_mis; 
  Pmis = Pmis - np.diag(np.diag(Pmis))
     
  Sum_rows_movOut_mis = [sum(movOut_mis[i]) for i in range(len(movOut_mis))]
  Sum_cols_movDrt_mis = [sum(movDrt_mis[:,i]) for i in range(len(movDrt_mis))]

  Rmis = Rmis + np.transpose (Sum_rows_movOut_mis)*dt - Sum_cols_movDrt_mis    
       
  Rnor_cleared = Rnor * (1 - np.exp(-clearance_rate * dt))
  Rmis_cleared = Rmis * (1 - np.exp(-clearance_rate * dt))
  #the probability of getting misfolded
  delta0 = 1 * trans_rate /ROIsize 
  misProb = 1 - np.exp( - Rmis * delta0 * dt ) # trans_rate: default
  #number of newly infected
  N_misfolded = Rnor * np.exp(- clearance_rate) * misProb 
  #update
  Rnor = Rnor - Rnor_cleared - N_misfolded + (synthesis_rate * syn_control) *dt
  Rmis = Rmis - Rmis_cleared + N_misfolded

  Rnor_all[: , t]= Rnor 
  Rmis_all[: , t] = Rmis 

  #uncomment the following lines if you want outputs of alpha-syn in
  #paths
  #Pnor_ave(:, :, t) = Pnor
  #Pmis_ave(:, :, t) = Pmis
  #additive noise


 return Rnor_all, Rmis_all, Pnor_all, Pmis_all, Rnor0, Pnor0



if __name__=="__main__":
  #data control healthy
 filename_control = os.listdir("Healthy")
 for filename in filename_control:
  subject = np.loadtxt(open("Healthy/"+filename, "rb"), delimiter=",")
  static_CM = static_connectivity_matrix(subject)
 
 #data alzheimer's disease
 filename_ad = os.listdir("Sick")
 i = 1
 for filename in filename_ad:
  subject = np.loadtxt(open("Sick/"+filename, "rb"), delimiter=",")
  static_CM = static_connectivity_matrix(subject)
 
  GNA = np.diag(static_CM)
  SNCA = GNA / np.linalg.norm(GNA)
  Rnor_all, Rmis_all, Pnor_all, Pmis_all, Rnor0, Pnor0 = Simulation(96, 1, 0.01, 10000 , GNA, SNCA , subject, subject, 4, 2, 1 , 1 , 0.5, 2)
  print('the number of normal beta_amyloid in regions for subject ',i,'is\n', Rnor_all,) 
  print('the number of misfolded beta_amyloid in regions for subject ',i,'is\n', Rmis_all,) 
  print('the number of normal beta_amyloid in paths could be memory-consuming for subject',i,'is\n', Pnor_all,) 
  print('the number of misfolded beta_amyloid in paths could be memory-consuming for subject',i,'is\n', Pmis_all_all,) 
  print('the population of normal agents in regions before pathogenic spreading\n',Rnor0)
  print('the population of normal agents in regions before pathogenic spreading\n',Pnor0)

