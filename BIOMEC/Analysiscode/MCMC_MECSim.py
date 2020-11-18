# -*- coding: utf-8 -*-
"""
Created on Tue May 29 14:46:49 2018

@author: Python
"""

from scipy.fftpack import rfft, irfft, rfftfreq
from pandas import read_csv, read_table, DataFrame
import numpy as np
import time
import os
from functools import partial
from multiprocessing import Pool
# Custom function imports
from Script_generator import *     # writing modules
from ML_signal_processing import *   #General ML_processing funcion
from MCMC_modules import *

t1 = time.time()

method_settings, data, Exp_data = globalinreader('global.txt') #Takes the input data and seperates to its sections

var, bandwidth, harm_weights, op_settings, datatype, MCMC_settings, scalvar = MCMC_Varibles(method_settings) # changes CMA_settings to usable data

curr, Exp_t = Exp_data_spliter(Exp_data, datatype) # splits the experimental input depending on input type

harm_weights = harm_weights_trans(harm_weights)  # reshapes harm_inputs into row

"""Add a varience propigator for the wieghts"""

spaces = data_finder(data)

AC_freq, AC_amp = AC_info(data,spaces)

DCAC_method = ACDC_method(AC_amp,op_settings)

if DCAC_method[0] == 1: # AC method {DC = 0, AC = 1}
    
    Ex_fft_res = rfft(curr)
    Ex_freq = rfftfreq(len(curr), d = Exp_t)   #defines frequency 
    n_window, f_window, g = ifft_windowing(AC_freq, Ex_freq, bandwidth, spaces)
    EX_hil_store, EX_N_ifft = Harmonic_gen(Ex_fft_res, n_window, spaces, bandwidth, g)   
    
    # adjusts the std in error for wieghts
    MCMC_settings[3] = wieghted_err_std(MCMC_settings[3],harm_weights)
     
else:
    
    pass  

# for all the messy stuff
space_holder = var[[4,5]].values   # retains values for use in printed values in np.array
var, scalvar = line_assigner(var,spaces,scalvar)
op_settings[1] = MCMC_settings[4]    # sets number of cores
parallelgen(op_settings)     # sets up the working files

# set up the prior and the initial values
phi0, covar, prior_range = prior_man(var)

if DCAC_method[0] == 1: # AC method
    
    # gets the MECSim settings set up to iterativeMECSIM through baysian code 
    MEC_set = {'data':data,'var':var,'spaces':spaces, 'DCAC_method':DCAC_method, 'Exp_data':EX_hil_store, 'harm_weights':harm_weights, 'bandwidth':bandwidth, 'scalvar':scalvar} # here to allow functional ACDC functionality
    
    # parallelization starts here
    # collects the MCMC settings into a dictionary
    MCMC_set = {'phi0':phi0,'covar':covar,'prior_range':prior_range,'MCMC_settings':MCMC_settings,'MEC_set':MEC_set}
    
    # defines a working directory # and a random number seed
    listin = listin_def(MCMC_settings[4])
    
    # Attaches kwargs to wrapper function
    func = partial(MCMC_Chain,**MCMC_set)
    
    # Multiprocessing call around Iterative MECSim
    with Pool(processes=MCMC_settings[4]) as p:
            
        multi_res = p.map(func, listin)     # only args we pass in is the counter for parralelization
        # note that results is a list of lists of lists
    
    # creates and generates the paralized output
    resultM, accep_list, Err = Para_out_setter(MCMC_settings[4],multi_res)
    
    #Compines the paralized outputs into a single output
    dist, accept_rate, Error = Para_combiner(resultM, accep_list, Err)
    
    MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = covarinence_stat(dist, MCMC_settings[2])
    
    MCMC_mean = np.append(MCMC_mean,1)      # sets up for working directory
    
    Perr = Iterative_MECSim(MCMC_mean, **MEC_set)
    
    Perr = Perr_ACtran(Perr, spaces) # translates Perr into a matrix form
      
else:    
    
    # gets the MECSim settings set up
    MEC_set = {'data':data,'var':var,'spaces':spaces, 'DCAC_method':DCAC_method, 'Exp_data':curr, 'scalvar':scalvar}
    
    # set settings to pass into Iterative MECSim
    
    # parallelization starts here
    # collects the MCMC settings into a dictionary
    MCMC_set = {'phi0':phi0,'covar':covar,'prior_range':prior_range,'MCMC_settings':MCMC_settings,'MEC_set':MEC_set}
    
    # defines a working directory # and a random number seed
    listin = listin_def(MCMC_settings[4])
    
    # Attaches kwargs to wrapper function
    func = partial(MCMC_Chain,**MCMC_set)
    
    # Multiprocessing call around Iterative MECSim
    with Pool(processes=MCMC_settings[4]) as p:
            
        multi_res = p.map(func, listin)     # only args we pass in is the counter for parralelization
        # note that results is a list of lists of lists
    
    # creates and generates the paralized output
    resultM, accep_list, Err = Para_out_setter(MCMC_settings[4],multi_res)
    
    #Compines the paralized outputs into a single output
    dist, accept_rate, Error = Para_combiner(resultM, accep_list, Err)
    
    MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = covarinence_stat(dist, MCMC_settings[2])
    
    MCMC_mean = np.append(MCMC_mean,1)      # sets up for working directory
    
    Perr = Iterative_MECSim(MCMC_mean, **MEC_set)
    

#Collects run time
MCMC_time = (time.time() - t1 )/60 # Time in minutes

MCMC_para_logger('MCMC_parallel_logger.txt', resultM,Err) # prints the multi chain output
MCMC_dist_logger('global_logger.txt', dist,Error)    # saves the distrabution to a text file 

# set up for printing
Niter = MCMC_settings[1] # gets the number of iterations

# ints CMA results in a .txt file THIS DOESN'T need to be like this desion inside a signular function
if DCAC_method[0] == 1: # AC method
    
    MCMC_output_printAC('global_out.txt',MCMC_time,Niter,accept_rate,MCMC_mean,colapsed_std,phi0,Perr,covar_std,cal_covar,corr_matrix)
    
else:
    
    MCMC_output_printAC('global_out.txt',MCMC_time,Niter,accept_rate,MCMC_mean,colapsed_std,phi0,Perr,covar_std,cal_covar,corr_matrix)
