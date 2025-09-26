# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 15:07:45 2018

@author: Luke Gundry

"""

# {PUT CMA_var_inserter, MECSim, fitter} in here

import cma
from numpy import random, zeros, arctanh, log10, log, dot, append, array, asarray, append
from pandas import DataFrame
from functools import wraps, partial
from multiprocessing import Pool
from ML_signal_processing import Iterative_MECSim
from Script_generator import iter_logger
from timeseries_modules import *
   
def originalCMA_var_handler(Iterative_MECSim0): # AC is needed so args can be used in the wrapper
    def CMA_Var(*args,**kwargs): #"""THIS IS GOING TO BE A PAIN # data,var,spaces,Ex_data """
        
        DCAC_method = kwargs.get('DCAC_method')
        var = kwargs.get('var')
        val_in = kwargs.pop('val_in') # remove val_in fro keywords 
        
        # Will need to put an if loop to not do a parallel
        
        Npop = len(args) # gets population iteration size 
        Ncore = DCAC_method[2] # gets the number of usable cores
        
        val_log = [] # creates a empty list for log values
        listin = [] # creates a blank list
        i = 0
        for input_val in args:
            
            Nallo = i % Ncore + 1  # allocation for parallel file
                
            input_val = append(input_val,Nallo)  # adds job reference number
            listin.append(input_val)    # adds the modified array to new list 
            i += 1
        
        Nvar = len(listin[0]) - 1   # number of varibles being modified
        
        Nin = len(listin)
        """parallel starts here FUCK FUCK FUCK"""
             
        i = 0 # CMA input counter
        while i != Npop: # this sets the CMA MODEL INPUT forif DCAC_method == 1: # AC method   
           j = 0
           while j != Nvar:
               if var.iloc[j][3] == 0:  #Non-logaritmic scalable constant
                   if listin[i][j] > 0:
                       listin[i][j] = (var.iloc[j][2] - var.iloc[j][0])*listin[i][j] + var.iloc[j][0] #Parameter scaler
                   else:
                       listin[i][j] = (var.iloc[j][0] - var.iloc[j][1])*listin[i][j] + var.iloc[j][0] #Parameter scaler
               elif var.iloc[j][3] == 1: #logaritmic scalable constant
                   if listin[i][j] > 0:
                       listin[i][j] = 10**((log10(var.iloc[j][0]/var.iloc[j][1]))*listin[i][j] + log10(var.iloc[j][0])) #Parameter scaler
                   else:
                       listin[i][j] = 10**((log10(var.iloc[j][2]/var.iloc[j][0]))*listin[i][j] + log10(var.iloc[j][0])) #Parameter scaler
               else:
                   print('please put a 0 or 1 in the right input')
               j += 1
           i += 1
        
        # multicore processsing
        
        # Attaches kwargs to wrapper function
        func = partial(Iterative_MECSim0,**kwargs)
        
        # Multiprocessing call around Iterative MECSim
        with Pool(processes=Ncore) as p:
            
            multi = p.map(func, listin)
               
        # prints the multivarible output in order
        val_out = []
        [val_out.append(s) for s in multi] # CMA-ES output is a list of NP-arrays
                       
            #p = Process(target=Iterative_MECSim0, args=(listin), kwargs = kwargs)
        
        # c = Iterative_MECSim0(x,**kwargs) # this needs to be a more rounded 
        
        if DCAC_method[1] == 1 and DCAC_method[0] == 1:   # use percentage fitter and AC method
            x = []  # space holder for logger
            out = [] # space holder for output c
            
            for output in val_out:
                
                c = output
                harm_weights = kwargs.get('harm_weights')
                
                c = dot(c,harm_weights)/len(harm_weights)
                
                out.append(c)  # adds the c value to a list
            
        else:
            out = [] # space holder for output c
            for output in val_out:
                  out.append(output[-1])
        
        iter_logger(listin, out) # append in loop then print to txt at c point redesign
    
        return out
    return CMA_Var


def CMA_MEC(data, var, val_in, spaces, DCAC_method, Exp_data, op_settings, scalvar, **kwargs):
    
    harm_weights = kwargs.get('harm_weights')
    
    if DCAC_method[0] == 1:   # AC method
        
        harm_weights = kwargs.get('harm_weights') # need that flexability function
        bandwidth = kwargs.get('bandwidth')
        MEC_set = {'data':data,'var':var,'val_in':val_in,'spaces':spaces, 'DCAC_method':DCAC_method, 'Exp_data':Exp_data, 'harm_weights':harm_weights, 'bandwidth':bandwidth, 'scalvar':scalvar}
    
    else:   # DC method
        
        MEC_set = {'data':data,'var':var,'val_in':val_in,'spaces':spaces, 'DCAC_method':DCAC_method, 'Exp_data':Exp_data, 'scalvar':scalvar}
    
    x = 1 
    N = random.randint(1000,9999)
    
    Dim = var.shape[0]  # Extracts the number of varibles being used
    popsize = int(3*log(Dim))  # calculates number of function evals per iteration (popsize)
    options = {'seed':N,'verb_disp':1,'tolx': op_settings[3],'bounds': [-x, x],'tolstagnation': int(100 + 100 * Dim**1.5 / (popsize*2)), 'maxiter': 1000} #Keyword assignment for CMA, 'tolfun':100 'tolx':0.05
   
    #Allocates wrapper  around iterative MECSim function
    CMA_Iterative_MECSim = CMA_var_handler(Iterative_MECSim)
    
    # initial std (Needs to be more personalized) output would be a list so just needs modification from there
    sigma = op_settings[4]   #33 probs better but sebs
    
    #logger = cma.CMAESDataLogger()      # records whats going on iterativly
    es = cma.CMAEvolutionStrategy(Dim*[0], sigma,options)# optimizes to a minimal point
    #logger = cma.CMADataLogger().register(es)
    "Might need an if loop for single core systems"
    #func = CMA_Iterative_MECSim(*num,**MEC_set)
    #with cma.fitness_transformations.EvalParallel(DCAC_method[2]*2) as eval_all:
    while not es.stop():
        # the multiprocessing needs to go in here
        X = es.ask()
        
        #input output list to do parallel in wrapper
        FIT = CMA_Iterative_MECSim(*X,**MEC_set) # Input a list of input lists for the wrapper to work on, as well as keyword settings
        
        es.tell(X, FIT) # Does the CMA-ES fitting
        
        #logger.add()
        #es.logger.disp([-1]) 
    #cma.plot()
    res = es.result
    
    return res

# Final translation 
def CMA_output(res, **kwargs):

    var = kwargs.get('var')
    DCAC_method = kwargs.get('DCAC_method')

    N = var.shape[0]
    var_out = list(zeros(N))
    
    res1 = res[0] # reassigns to best solution [5] better with noise
    
    j = 0
    while j != N:
        if var.iloc[j][3] == 0:  #Non-logaritmic scalable constant
            if res1[j] > 0:
                var_out[j] = (var.iloc[j][2] - var.iloc[j][0])*res1[j] + var.iloc[j][0] #Parameter scaler
            else:
                var_out[j] = (var.iloc[j][0] - var.iloc[j][1])*res1[j] + var.iloc[j][0] #Parameter scaler
            
            j += 1
                # function added in to adjust logarithmic values for cma-es
        elif var.iloc[j][3] == 1: #logaritmic scalable constant
                    
            if res1[j] > 0:
                var_out[j] = 10**((log10(var.iloc[j][0]/var.iloc[j][1]))*res1[j] + log10(var.iloc[j][0])) #Parameter scaler
            else:
                var_out[j] = 10**((log10(var.iloc[j][2]/var.iloc[j][0]))*res1[j] + log10(var.iloc[j][0])) #Parameter scaler
                
            j += 1
        else:
            print('please put a 0 or 1 in the right input')
            
    # here so we can allocate mean out to run on a processor

    EXcurr = kwargs.get('Exp_data')
    if DCAC_method[0] == 1:   # AC method
        
        harm_weights = kwargs.get('harm_weights') # need that flexability function
        bandwidth = kwargs.get('bandwidth')

        Scurr = Iterative_MECSim_Curr(var_out, **kwargs)
        Loss = HE_percentagefit(EXcurr, Scurr,**kwargs)

    else:   # DC method

        Scurr = Iterative_MECSim_Curr(var_out, **kwargs)
        Loss = TotCurr_diffsquare(EXcurr, Scurr)
    
    return var_out, Loss

# NEEDS WORK
def CMA_meanval(res, var):
    
    N = var.shape[0]
    var_out = list(zeros(N))
    
    #gets mean and std from CMA
    values = res[5]
    std = res[6]
    
    j = 0
    while j != N:
        if var.iloc[j][3] == 0:  #Non-logaritmic scalable constant
            if values[j] > 0:
                var_out[j] = (var.iloc[j][2] - var.iloc[j][0])*values[j] + var.iloc[j][0] #Parameter scaler
            else:
                var_out[j] = (var.iloc[j][0] - var.iloc[j][1])*values[j] + var.iloc[j][0] #Parameter scaler
            
            j += 1
                # function added in to adjust logarithmic values for cma-es
        elif var.iloc[j][3] == 1: #logaritmic scalable constant

            if values[j] > 0:
                var_out[j] = 10**((log10(var.iloc[j][0]/var.iloc[j][1]))*values[j] + log10(var.iloc[j][0])) #Parameter scaler
            else:
                var_out[j] = 10**((log10(var.iloc[j][2]/var.iloc[j][0]))*values[j] + log10(var.iloc[j][0])) #Parameter scaler
                
            j += 1
        else:
            print('please put a 0 or 1 in the right input')
    
    means = var_out
    #Insert STD scaler here
    var_out = list(zeros(N))
    j = 0
    while j != N:
        if var.iloc[j][3] == 0:  #Non-logaritmic scalable constant
            if values[j] > 0:
                var_out[j] = (var.iloc[j][2] - var.iloc[j][0])*std[j] #Parameter scaler
            else:
                var_out[j] = (var.iloc[j][0] - var.iloc[j][1])*std[j] #Parameter scaler
            
            j += 1
                # function added in to adjust logarithmic values for cma-es
        
        elif var.iloc[j][3] == 1: #logaritmic scalable constant ( NOT NESSACARY BUT WHAT I HAD BEFORE WAS WRONG)
                    
            if values[j] > 0:
                var_out[j] = (var.iloc[j][2] - var.iloc[j][0])*std[j] #Parameter scaler
            else:
                var_out[j] = (var.iloc[j][0] - var.iloc[j][1])*std[j]#Parameter scaler
            
            j += 1
        else:
            print('please put a 0 or 1 in the right input')
        
    standev = var_out
    meanprop = [means,standev]    
       
    return meanprop
    