# -*- coding: utf-8 -*-
"""
Created on Tue May 22 00:38:20 2018

MCMC modules for numerical analysis of MECSim data

@author: Luke Gundry
"""

import numpy as np
from pandas import DataFrame
from functools import wraps
import datetime
import timeseries_modules as tpseries
import time
from ML_signal_processing import Iterative_MECSim
from Script_generator import format_e

def prior_man(var):
    
    N = var.shape[0]    # gets the # of varibles
    
    # preallocation
    phi0 = np.zeros(N)
    covar = np.zeros((N,N))
    prior_range = np.zeros((N,2))
    
    i = 0
    while i != N:
        
        phi0[i] = var.values[i][0]
        prior_range[i,0] = var.iloc[i][1]
        prior_range[i,1] = var.iloc[i][2]
        covar[i,i] = (prior_range[i,1] - prior_range[i,0])**2/12 # varience of a uniform distrabution
        
        i += 1
        
    return phi0, covar, prior_range

# function for handling DCAC_methods for iterative mecsim
def sumerr_MECS(val_in,counter,**kwargs):
    
    DCAC_method = kwargs.get('DCAC_method')
    
    val_in = np.append(val_in,counter) # sets the working directory (works for now)
    
    c = Iterative_MECSim(val_in,**kwargs)
    
    if DCAC_method[1] == 1 and  DCAC_method[0] == 1:   # use percentage fitter and AC method
        harm_weights = kwargs.get('harm_weights')
        c = np.dot(c,harm_weights)/len(harm_weights) # gets adverage per harmonic wieght
    else:
        pass
        
    return c

def prior_test(val_in, prior_range):
    
    N = len(val_in)
    prob = 1
    i = 0
    while i != N:
        if val_in[i] >= prior_range[i,0] and  val_in[i] <= prior_range[i,1]:
            prob *= 1
        else:
            prob = 0
            break  # cancels while loop with prob = 0 
        i += 1
    
    return prob

# the calculation step between the priori and final sampler
def MCMC_burnin_sample(phi0, covar, prior_range, noise, Nsamp,counter, **kwargs):
    
    # allocates error
    Err = []
    
    ao = 1
    t = 1 # iteration counter
    N = Nsamp*len(phi0) # max counter
    phit = [np.array(phi0)]    # sets the initial mean {MAY NOT BE RIGHT}
    Naccept = 0         # acceptance counter
    
    curr_fit = sumerr_MECS(phi0, counter,**kwargs) # experimental data gets passed in kwargs
    Err.append(curr_fit)
    
    
    
    while t != N + 1:
        
        phi_trail = np.random.multivariate_normal(phit[-1],ao*covar)
        prob = prior_test(phi_trail, prior_range) # use this to fitler experimental results with respect to prior
    
        if prob != 0: # checks to see if position is possible
            
            phi_test_fit = sumerr_MECS(phi_trail, counter,**kwargs) # tests the possible solution
            logtrial = -(t+1)*np.log(noise) - (1/(2*noise**2))*phi_test_fit
            
            # Probability comparison code
            comp = np.exp(logtrial - (-(t)*np.log(noise) - (1/(2*noise**2))*curr_fit))
            r = min(1,comp)
            
            u = np.random.uniform(0,1)  # generates random uniform distrabution 
            if u < r:
                
                phit.append(phi_trail) 
                curr_fit = phi_test_fit  # resets the sum
                Err.append(curr_fit)
                Naccept += 1 # accept rate
                
            else:
                phit.append(phit[-1])
                Err.append(curr_fit)

        else:
            phit.append(phit[-1]) # accept the current state
            Err.append(curr_fit)
        
        t += 1
    
    return phit, Naccept, Err

#MCMC metropolis hastings algoritm with adaptive covarience 
def MCMC_metro_hast(phit, covar, prior_range, noise, Nsamp, NMax, Naccept, counter,**kwargs):
    
    # Error record
    Err = []
    
    # set up for the loop
    N = Nsamp*len(phit[0])
    t = N + 1   # iteration start point
    
    # define important parameters
    mean = np.array(phit[0])  # watch THIS 
    a = 1
    curr_fit = sumerr_MECS(phit[-1], counter,**kwargs)   # sets up s0 point
    
    while t != NMax:
        
        s = t - N
        gamma = (s + 1)**(-0.6)
        
        phi_trail = np.random.multivariate_normal(phit[-1],a*covar)  # gets the trail solution
        
        prob = prior_test(phi_trail, prior_range)   # Tests if solution is possible assuming uniform distrabution
        
        if prob != 0: # checks to see if position is possible
            
            phi_test_fit = sumerr_MECS(phi_trail, counter,**kwargs) # tests the possible solution
            logtrial = -(t+1)*np.log(noise) - (1/(2*noise**2))*phi_test_fit
            
            # Probability comparison code
            comp = np.exp(logtrial - (-(t)*np.log(noise) - (1/(2*noise**2))*curr_fit))
            r = min(1,comp)
            
            u = np.random.uniform(0,1)  # generates random uniform distrabution 
            if u < r:
                
                phit.append(phi_trail) 
                curr_fit = phi_test_fit  # resets the sum
                Err.append(curr_fit)
                accepted = 1
                Naccept += 1    # counter for num accepted
                
            else:
                phit.append(phit[-1])
                Err.append(curr_fit)
                accepted = 0
                
        else:
            phit.append(phit[-1]) # accept the current state
            Err.append(curr_fit)
            accepted = 0
        
        # iteration probability modifiers
        covar = ((1 -gamma)*covar + gamma*np.outer((phit[-1] - mean),(phit[-1] - mean)))
        mean = (1 -gamma)*mean + gamma*phit[-1]
        a = np.exp(np.log(a) + gamma*(accepted-0.25)) # 0.25 will need to be a varible
        
        t += 1
    
    rate = Naccept/t    # Calculates chain acceptance rate 
    
    return phit, rate, Err

# Extracts Covarience statistics from MCMC output
def covarinence_stat(phit, burnin):
    
    N = len(phit[0,:])    # extracts the cma best fit values
    
    Nmax = phit.shape[0]    # Gets the overall number of iterations
    
    MCMC_mean = []
    i = 0
    while i != N:
        
        MCMC_mean.append(np.mean(phit[int(Nmax*burnin):,i]))    # calculates and adds the mean value
        i += 1
    
    # single value standard deviations
    colapsed_std = []
    i = 0
    while i != N:
        
        colapsed_std.append(np.std(phit[int(Nmax*burnin):,i]))    # calculates and adds the mean value
        i += 1
        
    # multivarible dist stats
    cal_covar = np.cov(np.transpose(phit[int(Nmax*burnin):,:]))    # gets the covarence matrix of the 

    covar_std = np.sqrt(np.diag(cal_covar)) # extracts the sigmas from covar matrix
    
    sig_sq = np.outer(covar_std,covar_std)      # gets a sigma squared matrix
    
    corr_matrix = cal_covar/sig_sq  # Gets the coleation matrix
    
    return MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix

# Extracts the settings for mecsim from the block
def MCMC_Varibles(MCMC_settings):
    
    nvar = int(MCMC_settings.iloc[0][0])
    n = MCMC_settings[0].shape
    n = int(n[0])
    nscal = int(MCMC_settings.iloc[nvar + 1][0])
    
    # Gets the scaling 
    scalvarpre = MCMC_settings[nvar+ 1:nvar+ 2+nscal] # gets the num of scaling to
    
    nband = int((n-9-nvar - 1 - nscal)/2)
    
    scalvar = [[int(scalvarpre.iloc[0][0])]]
    i = 0 
    x = scalvarpre[0].str.split(',',expand=True).astype('float')
    x = x.values
    
    while i != scalvar[0][0]:
        scalvar.append(x[i + 1,:])
        i += 1
    
    #Seperates the CMA input data
    var = MCMC_settings[1:nvar + 1]
    bandwidth = MCMC_settings[nvar+ 1 + 1 + nscal:nvar+ 1 +nband + 1 + nscal]
    harm_weights = MCMC_settings[nvar+ 1 +nband + 1 + nscal:nvar+ 1 + 2*nband + 1 + nscal]
    
    #Changes the seperated input values into 
    var = var[0].str.split(',',expand=True).astype('float')
    bandwidth = bandwidth[0].str.split(',',expand=True).astype('float')
    bandwidth = bandwidth.values    # changes from a df to np.array
    harm_weights = harm_weights[0].str.split(',',expand=True).astype('float')
    harm_weights = harm_weights.values
    
    #Optimization settings
    Np = nvar + 2*nband + nscal +1
    op_settings = [0,0,0]     #prealocation
    op_settings[0] = int(MCMC_settings.iloc[Np + 1][0])   # This one setts the optimization method
    datatype = int(MCMC_settings.iloc[Np + 2][0]) #Gets the type of data were compairing against
    op_settings[1] = int(MCMC_settings.iloc[Np + 3][0]) # this one does the num cores
    op_settings[2] = 2**float(MCMC_settings.iloc[Np + 4][0])  # sets the number of datapoints to compair in the current

    #MCMC settings
    bay_settings = [0,0,0,0,0]
    bay_settings[0] = int(MCMC_settings.iloc[Np + 5][0])     # MCMC initial algoritm trail's per varible
    bay_settings[1] = int(MCMC_settings.iloc[Np + 6][0])     # number of trail for overall chain
    bay_settings[2] = float(MCMC_settings.iloc[Np + 7][0])     # burnin period as ratio of chain legth
    # treatment of the multi noise 
    bay_settings[3] = MCMC_settings.iloc[Np + 8]
    bay_settings[3] = bay_settings[3].str.split(',',expand=True).astype('float')
    bay_settings[3] = bay_settings[3].values[0,:]    # noise of fit
    
    bay_settings[4] = int(MCMC_settings.iloc[Np + 9][0])     # number of chain to run at once
    
    return var, bandwidth, harm_weights, op_settings, datatype, bay_settings, scalvar 

def MCMC_Chain(*args,**kwargs):
    
    # gets the input arguments from kwargs
    phi0 = kwargs.pop('phi0')
    covar = kwargs.pop('covar')
    prior_range = kwargs.pop('prior_range')
    MCMC_settings = kwargs.pop('MCMC_settings')
    MEC_set = kwargs.pop('MEC_set')
    
    #sets chain specific settings
    counter = args[0][0]
    np.random.seed(args[0][1]) # sets the random number seed   
    
    # reset the lengths of the chain with respect to the parallelization    
    Nlength = int(MCMC_settings[1]/MCMC_settings[4])    # length of overall chain
    Nintvar = int(MCMC_settings[0]/MCMC_settings[4])    # legth of initial sampling
    
    # intial sampling of the MCMC algoritm
    burn_dist, Naccept, Err1 = MCMC_burnin_sample(phi0, covar, prior_range, MCMC_settings[3], Nintvar, counter,**MEC_set) # ,CMA_var_handler
    # bulk adaptive sampling
    dist, accept_rate, Err2 = MCMC_metro_hast(burn_dist, covar, prior_range, MCMC_settings[3], Nintvar, Nlength, Naccept, counter,**MEC_set)    
    
    Err = Err1 + Err2
    
    return dist, accept_rate,Err

# changes the noise from an unscaled approximation to what we use 
def wieghted_err_std(noise,harmweights):
    
    N = len(harmweights)    # sets number of signals
    noise_new = 0 # sets the adjusted std
    
    # gets the number of wieghts per harmonic
    NHz = (N-1)/(len(noise)-1)
    
    # DC signals
    noise_new += (harmweights[0]*noise[0])**2

    # AC signals
    j = 1
    k = 1
    while j != N:   # goes through freq
        i = 1
        while i != NHz+1:     # goes through harm set
            
            noise_new += (harmweights[int(NHz*(k-1) + i)]*noise[k])**2
            i += 1
            j += 1
        k += 1
    noise_new = np.sqrt(noise_new)/N    # final adjustment
    
    return noise_new

# prints the DC output
def MCMC_output_printDC(filename,MCMC_time,Niter,accept_rate,MCMC_mean,colapsed_std,phi0,Perr,covar_std,cal_covar,corr_matrix):
   
    N = len(phi0)
    
    f = open(filename,"w+")    # creates a blank .inp to be written
    
    #writes the input 
    f.write('Input file name: %s\n' %(filename))
    # date of the colour
    f.write('Date: %s\n' %(datetime.datetime.today().strftime('%d-%m-%Y')))
    # inserts time
    f.write('Completetion time: %f\n\n' %(MCMC_time))
    
    f.write('number of function Iterations: %i\n' %(Niter))
    
    f.write('Acceptance rate of MCMC: %f\n\n' %(accept_rate))
    
    # Statistics
    f.write('Statistical values\n')
    
    # mean_var_out = str(mean_var_out)
    f.write('MCMC 1D Values Mean parametervalues (mean,std):\n')
    i = 0
    while i != N:
        
        m = format_e(MCMC_mean[i])
        s = format_e(colapsed_std[i])
        
        f.write('Var %d: %s,%s\n' %(int(i+1),m,s))
        i +=1
        
    f.write('\n')
    
    f.write('Comparison between input mean and MCMC mean:\n')
    i = 0
    while i != N:
        
        m = MCMC_mean[i]/phi0[i]
        f.write('Var %d: %.5f\n' %(int(i+1),m))
        i +=1
        
    f.write('\n')
    
    # Perr different to AC
    x = format_e(Perr)
    f.write('Percentage Error: %s\n\n' %(x))
    
    f.write('Multivariate normal distribution statistics\n\nCovarience standard deviation\n')
    i = 0
    while i != N:
        f.write('Var %d: %s\n' %(int(i+1),format_e(covar_std[i])))
        i +=1
        
    f.write('\n')    #seperator
    
    f.write('Covarience Matrix\n')
    
    i = 0
    while i != N:
        j = 0
        while j != N:
            f.write('%s\t' %(format_e(cal_covar[i,j])))
            j += 1
        f.write('\n')
        i += 1
    f.write('\n')    #seperator
    
    f.write('Correlation Matrix\n')
    i = 0
    while i != N:
        j = 0
        while j != N:
            f.write('%s\t' %(format_e(corr_matrix[i,j])))
            j += 1
        f.write('\n')
        i += 1
    f.write('\n')
    
    f.close()
    
    return

# prints the AC output to text 
def MCMC_output_printAC(filename,MCMC_time,Niter,accept_rate,MCMC_mean,colapsed_std,phi0,Perr,covar_std,cal_covar,corr_matrix):
   
    N = len(phi0)
    
    f = open(filename,"w+")    # creates a blank .inp to be written
    
    #writes the input 
    f.write('Input file name: %s\n' %(filename))
    # date of the colour
    f.write('Date: %s\n' %(datetime.datetime.today().strftime('%d-%m-%Y')))
    # inserts time
    f.write('Completetion time: %f\n\n' %(MCMC_time))
    
    f.write('number of function Iterations: %i\n' %(Niter))
    
    f.write('Acceptance rate of MCMC: %f\n\n' %(accept_rate))
    
    # Statistics
    f.write('Statistical values\n')
    
    # mean_var_out = str(mean_var_out)
    f.write('MCMC 1D Values Mean parametervalues (mean,std):\n')
    i = 0
    while i != N:
        
        m = format_e(MCMC_mean[i])
        s = format_e(colapsed_std[i])
        
        f.write('Var %d: %s,%s\n' %(int(i+1),m,s))
        i +=1
        
    f.write('\n')
    
    f.write('Comparison between input mean and MCMC mean:\n')
    i = 0
    while i != N:
        
        m = MCMC_mean[i]/phi0[i]
        f.write('Var %d: %.5f\n' %(int(i+1),m))
        i +=1
        
    f.write('\n')
    
    f.write('Percentage Error Values for each harmonic (lowest-highest freq):\n')
    H = np.count_nonzero(Perr)   #Gets the number of fittered haarmonics
    Sum_err = np.sum(Perr)/H  #Retreaves Sum error in all harmonics
    x = np.array_str(Perr)  #Converts Perr to a list 
    f.write(x)       #Might need to be fixed
    f.write('\n')    #seperator
    f.write('Total prcentage error of all harmonics: %.5f\n\n' %(Sum_err))
    
    f.write('Multivariate normal distribution statistics\n\nCovarience standard deviation\n')
    
    i = 0
    while i != N:
        f.write('Var %d: %s\n' %(int(i+1),format_e(covar_std[i])))
        i +=1
        
    f.write('\n')    #seperator
    
    f.write('Covarience Matrix\n')
    
    i = 0
    while i != N:
        j = 0
        while j != N:
            f.write('%s\t' %(format_e(cal_covar[i,j])))
            j += 1
        f.write('\n')
        i += 1
    f.write('\n')    #seperator
    
    f.write('Correlation Matrix\n')
    i = 0
    while i != N:
        j = 0
        while j != N:
            f.write('%s\t' %(format_e(corr_matrix[i,j])))
            j += 1
        f.write('\n')
        i += 1
    f.write('\n')
    
    f.close()
    
    return

# Prints the intergrated values of phit to a text file
def MCMC_dist_logger(filename, phit,Error):
    
    f = open(filename,"w+")
    
    Nv = phit.shape[0]      # number of points 
    Nd = phit.shape[1]      # number of varibles
    
    i = 0
    while i != Nv:
        j = 0
        while j != Nd:
            
            x = format_e(phit[i,j])
            f.write('%s\t' %(x))
            j += 1
        
        f.write('%s\n' %(format_e(Error[i])))   # prints error
        
        i += 1
    
    f.close()
    
    return

# prints the multichain output to another txt file
def MCMC_para_logger(filename, multi_res,Err):
    
    # remember multi_res is a list of lists of lists
    
    f = open(filename,"w+")
    
    Nchain = len(multi_res)
    Npoint = len(multi_res[0])
    Nvar = len(multi_res[0][0])
    
    i = 0
    while i != Npoint:
        j = 0
        while j != Nchain:
            k = 0
            while k != Nvar:
                
                f.write('%s\t' %(format_e(multi_res[j][i][k])))  
                k +=1 
            
            f.write('%s\t|\t' %(format_e(Err[j][i]))) # chain seperatof
            j += 1
        f.write('\n')
        i += 1
    
    f.close()
    
    return


# Prints the intergrated values of phit to a text file
def MCMC_dist_logger2(filename, phit):
    f = open(filename, "w+")

    Nv = phit.shape[0]  # number of points
    Nd = phit.shape[1]  # number of varibles

    i = 0
    while i != Nv:
        j = 0
        while j != Nd:
            x = format_e(phit[i, j])
            f.write('%s\t' % (x))
            j += 1

        f.write('\n')  # prints error

        i += 1

    f.close()

    return


# prints the multichain output to another txt file minus Error
def MCMC_para_logger2(filename, multi_res):
    # remember multi_res is a list of lists of lists

    f = open(filename, "w+")

    Nchain = len(multi_res)
    Npoint = len(multi_res[0])
    Nvar = len(multi_res[0][0])

    i = 0
    while i != Npoint:
        j = 0
        while j != Nchain:
            k = 0
            while k != Nvar:
                f.write('%s\t' % (format_e(multi_res[j][i][k])))
                k += 1

            f.write('|\t')  # chain seperatof
            j += 1
        f.write('\n')
        i += 1

    f.close()

    return

# creates listin of working file allocation and random seed
def listin_def(Ncores):
    
    listin = []
    i = 1
    while i != Ncores + 1:
        listin.append([i,np.random.randint(1,10000)]) # creates a list of working directorys
        i += 1
    
    return listin

# creates and generates the paralized output
def Para_out_setter(Ncores,para_out):
    
    accep_list = [] 
    resultM = []
    Err = []
    i = 0
    while i != Ncores:
        
        resultM.append(para_out[i][0])     # paralized varibles
        accep_list.append(para_out[i][1])  # acceptance rate
        Err.append(para_out[i][2])     # para Error
        i += 1   
    
    return resultM, accep_list, Err

#Compines the paralized outputs into a single output
def Para_combiner(resultM, accep_list, para_Err):
    
    Ncores = len(resultM)
    
    # Nvar and Nchain arn't needed
    dist = np.zeros((Ncores*len(resultM[0]),len(resultM[0][0])))    # Nchain = resultM[0], Nvar = resultM[0][0]
    Error = np.zeros(Ncores*len(resultM[0]))
    # combines the acceptance/MCMC rate
    i = 0
    while i != Ncores:
    
        dist[i::Ncores,:] = np.array(resultM[i][:])
        Error[i::Ncores] = para_Err[i]
        i += 1
    
    accept_rate = sum(accep_list)/Ncores
    
    return dist, accept_rate, Error


def Gel_rub_test(filename,resultM, Err, burnin):
    # make everthing an array
    R_hat = []
    NN = []

    Nvar = len(resultM[0][0])
    N = len(resultM[0])
    Nx = int(0.1 * N)
    NN.append(Nx)
    Nchain = len(resultM)

    mean_var = np.zeros(Nchain)
    mean_collect = []
    totmean = []

    # apples burning
    i = 0
    while i != Nchain:
        resultM[i] = np.array(resultM[i][int(burnin * N):N])
        i += 1

    # update burnin legth
    N = len(resultM[0][:, 0])

    k = 0
    j = 0
    while NN[k] < N:

        # calc chain means
        mean_collect = []
        totmean = np.zeros(Nvar)
        j = 0
        while j != Nvar:
            i = 0
            mean_var = np.zeros(Nchain)
            while i != Nchain:
                mean_var[i] = np.mean(resultM[i][0:NN[k], j])
                i += 1
            mean_collect.append(mean_var)
            totmean[j] = np.mean(mean_var)
            j += 1

        # Calculate B
        i = 0
        Bh = np.zeros(Nvar)
        while i != Nvar:
            m = 0
            while m != Nchain:
                Bh[i] += (totmean[i] - mean_collect[i][m]) ** 2
                m += 1
            i += 1

        Bh = NN[k] / (Nchain - 1) * Bh  # adds degrees of freedom modifier

        # calculate  W
        Wh = np.zeros(Nvar)
        m = 0
        while m != Nvar:
            j = 0
            Win2 = 0
            while j != Nchain:
                i = 0
                Win1 = 0
                while i != NN[k]:
                    Win1 += (resultM[j][i, m] - mean_collect[m][j]) ** 2
                    i += 1
                Win2 += Win1 / (NN[k] - 1)
                j += 1
            Wh[m] = Win2 / Nchain
            m += 1

        # calc var
        Var = (NN[k] - 1) / NN[k] * Wh + Bh / NN[k]
        # calc R_hat
        R_hat.append(np.sqrt(Var / Wh))
        # store R-hat with respective iteration

        NN.append(NN[-1] + Nx)
        k += 1
    # put R_hat vector in a txt file

    # sets the last point in the array
    NN[-1] = N

    # calc chain means
    mean_collect = []
    totmean = np.zeros(Nvar)
    j = 0
    while j != Nvar:
        i = 0
        mean_var = np.zeros(Nchain)
        while i != Nchain:
            mean_var[i] = np.mean(resultM[i][0:N, j])
            i += 1
        mean_collect.append(mean_var)
        totmean[j] = np.mean(mean_var)
        j += 1

    # Calculate B
    i = 0
    Bh = np.zeros(Nvar)
    while i != Nvar:
        m = 0
        while m != Nchain:
            Bh[i] += (totmean[i] - mean_collect[i][m]) ** 2
            m += 1
        i += 1

    Bh = NN[-1] / (Nchain - 1) * Bh  # adds degrees of freedom modifier

    # calculate  W
    Wh = np.zeros(Nvar)
    m = 0
    while m != Nvar:
        j = 0
        Win2 = 0
        while j != Nchain:
            i = 0
            Win1 = 0
            while i != NN[-1]:
                Win1 += (resultM[j][i, m] - mean_collect[m][j]) ** 2
                i += 1
            Win2 += 1 / (NN[-1] - 1) * Win1
            j += 1
        Wh[m] = Win2 / Nchain
        m += 1

    # calc var
    Var = (NN[-1] - 1) / NN[-1] * Wh + Bh / NN[-1]
    # calc R_hat
    R_hat.append(np.sqrt(Var / Wh))

    # prints the R_hat statistics
    R_hat_print(filename,NN, R_hat)

    return R_hat[-1]


# Prints R-Hat values
def R_hat_print(filename,NN, R_hat):
    # pas final to function
    f = open(filename, "w")
    N = len(R_hat[0])

    f.write('N-burn\t')

    i = 0
    while i != N - 1:
        s = 'Var %i\t' % i
        f.write(s)
        i += 1
    f.write('Err\n')

    i = 0
    while i != len(R_hat):
        s = '%i\t' % NN[i]
        f.write(s)
        j = 0
        while j != N:
            s = '%f\t' % (R_hat[i][j])
            f.write(s)
            j += 1
        f.write('\n')
        i += 1
    f.close()

    return

# sets up the priors for MCMC
def prior_set(var, noise,Noisebool):

    N = len(var.iloc[:][0])

    # priors = [[min/mean,max/std,Uniform = 0/ normal = 1]]
    L = []  # lower bounds
    H = []  # higher bounds
    type = []
    prior = []
    i = 0
    while i != N:
        if var.iloc[i][1] == var.iloc[i][2]: # std case for normal distrabution priors
            L.append(var.iloc[i][1])
            H.append(var.iloc[i][2])
            type.append(1)
        else:
            L.append(var.iloc[i][1])
            H.append(var.iloc[i][2])
            type.append(0)
        i += 1

    if Noisebool:
        L.append(noise[1])
        H.append(noise[2])
        type.append(1)
    else:
        pass

    prior = [L,H,type]

    return prior

# checks if the input to MCMC is possible
def prior_check(val_in, prior):

    #priors = [[min/mean,max/std,Uniform = 0/ normal = 1]]

    N = len(val_in)
    out = True # passes the test

    i = 0
    while i != N:
        if prior[i][2] == 1:
            print('hasn,t been set up yet')
            exit()
        else:
            if prior[i][0] < val_in[i] < prior[i][1]:
                pass
            else: # fails the test
                out = False
                i = N - 1   # just skip the rest
        i += 1



    return out


# Final translation
def MCMC_harmfitpercentage(MCMC_mean,EXcurr,AC_method, **kwargs):

    if AC_method:  # AC method

        harm_weights = kwargs.get('harm_weights')  # need that flexability function
        bandwidth = kwargs.get('bandwidth')

        Scurr = tpseries.Iterative_MECSim_Curr(MCMC_mean, **kwargs)
        Loss = tpseries.HE_percentagefit(EXcurr, Scurr, **kwargs)

    else:  # DC method

        Scurr = tpseries.Iterative_MECSim_Curr(MCMC_mean, **kwargs)
        Loss = tpseries.TotCurr_diffsquare(EXcurr, Scurr)

    return Loss
