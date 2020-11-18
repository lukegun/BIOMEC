# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 14:55:49 2018

@author: Luke Gundry
"""
"""
LIST OF ALL MODGULS USED
scipy
pandas
numpy
time
functools
cma
os
subprocess
"""

from scipy.fftpack import rfft, irfft, rfftfreq
from pandas import read_csv, read_table, DataFrame
import numpy as np
import time
# Custom function imports
from Script_generator import *     # writing modules
from ML_signal_processing import *   #General ML_processing funciont
from CMA_modules import *
import os
#val_in = [0,1.0e-6,1.0e-5,1000,.5,0]    # draft input for CMA calling


t1 = time.time()

#stuff to move inside the directory
cwd = os.getcwd()
s = '%s/' %(cwd)
os.chdir(s) # changes directory

# Collects input data
CMA_settings, data, Exp_data = globalinreader('global.txt') #Takes the input data and seperates to its sections

var, bandwidth, harm_weights, op_settings, datatype, scalvar = CMA_Varibles(CMA_settings) # changes CMA_settings to usable data

curr, Exp_t = Exp_data_spliter(Exp_data, datatype) # splits the experimental input depending on input type

harm_weights = harm_weights_trans(harm_weights)  # reshapes harm_inputs into row

# pre-allocates the initial starting conditions
val_in = list(var.iloc[:][0])

spaces = data_finder(data)

AC_freq, AC_amp = AC_info(data,spaces)

DCAC_method = ACDC_method(AC_amp,op_settings)

if DCAC_method[0] == 1: # AC method {DC = 0, AC = 1}
    
    Ex_fft_res = rfft(curr)
    Ex_freq = rfftfreq(len(curr), d = Exp_t)   #defines frequency 
    n_window, f_window, g = ifft_windowing(AC_freq, Ex_freq, bandwidth, spaces)
    EX_hil_store, EX_N_ifft = Harmonic_gen(Ex_fft_res, n_window, spaces, bandwidth, g)   

else: 
    pass  

space_holder = var[[4,5]].values   # retains values for use in printed values in np.array
var, scalvar = line_assigner(var,spaces,scalvar)

#sets up for the 
logger("global_logger.txt",val_in) # creates a blank file for recording varibles
parallelgen(op_settings) # Generates the paralell folders

if DCAC_method[0] == 1: # AC method
    
    freq_set = {'harm_weights':harm_weights,'bandwidth':bandwidth}  # here to allow functional ACDC functionality
    
    res = CMA_MEC(data, var, val_in, spaces, DCAC_method, EX_hil_store, op_settings, scalvar, **freq_set) # ,CMA_var_handler
   
    var_out, c = CMA_output(res, data, var, spaces, DCAC_method, EX_hil_store, op_settings, scalvar, **freq_set)
    
    DCAC_method[1] = 1 # sets up to use the Perr method    
    var_out, Perr = CMA_output(res, data, var, spaces, DCAC_method, EX_hil_store, op_settings, scalvar, **freq_set) # here Perr is a row
    
    #Perr = Perr[1] #extracts the Harmonic wieght for AC
    
    Perr = Perr_ACtran(Perr, spaces) # translates Perr into a matrix form
    
    
else:    
    
    # Exp_data input has to be current only
    res = CMA_MEC(data, var, val_in, spaces, DCAC_method, curr, op_settings, scalvar) #,CMA_var_handler
    
    # gets a meaningfull output
    var_out, c = CMA_output(res, data, var, spaces, DCAC_method, curr, op_settings, scalvar)

    DCAC_method[1] = 1    # sets up to use the Perr method 
    var_out, Perr = CMA_output(res, data, var, spaces, DCAC_method, curr, op_settings, scalvar)
    
    #Perr = Perr[1] #extracts the fit for dc
    
#Collects run time  
CMA_time = (time.time() - t1 )/60 # Time in minutes

# ints CMA results in a .txt file THIS DOESN'T need to be like this desion inside a signular function
if DCAC_method[0] == 1: # AC method
    
    values = res[5]
    mean_var_out = CMA_meanval(res, data, var, spaces, DCAC_method, curr, harm_weights)
    s = CMA_output_printAC(CMA_time,space_holder, var_out,Perr,res,c,mean_var_out)
    
else:
    
    values = res[5]
    std = values = res[6]
    mean_var_out = CMA_meanval(res, data, var, spaces, DCAC_method, curr, harm_weights)
    s = CMA_output_printDC(CMA_time,space_holder, var_out,Perr,res,c,mean_var_out)

#Might want something to delete workingdir#