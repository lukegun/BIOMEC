"""
Modules for PINTS MECSim interface

Aurthor: Luke Gundry
Date: 16/7/19


"""

import cma
import pandas as pd
import numpy as np
import pints
from multiprocessing import Pool
import ML_signal_processing as MLsp
#Iterative_MECSim, Iterative_MECSim_Curr
from Script_generator import iter_logger
from MCMC_modules import *
import matplotlib.pyplot as plt

# all the functions within the system
from CMAES import *
from ADMCMC import *

# reads the initial input file
def globalinreader(filename):
    Input_file = filename  # input("Name of the textfile input:") # will need to be modified to be used with the
    globalin = pd.read_csv(Input_file, sep='    ', index_col=False, header=None, comment='!')

    nsettings = globalin[globalin[0].str.match('MECSim settings')].index.values  # gets the CMA settings seperator
    nExp = globalin[globalin[0].str.match('Experimental settings')].index.values

    # splits the inpfile into its important sections
    settings = globalin[0:nsettings[0]]
    data = globalin[nsettings[0] + 1:nExp[0]]  # will need to be fixed in future updates to truncate at exdata
    Exp_data = globalin[nExp[0] + 1::]  # experimental data

    # This needs to be moved till after CMA-Settings
    # Exp_data = Exp_data[0].str.split('  ',expand=True).astype('float')# \t for tab sep files # '  ' for mecsim

    return settings, data, Exp_data

# Gets out the starting capacitance
def capacatance0(data,spaces,var):

    capacatance = []
    cap = False
    x = var.iloc[:][4].values

    i = 0
    while i != len(x):
        if 60 > int(x[i]) > 50:
            cap = True
        else:
            pass
        i += 1

    Csol = int(data.iloc[spaces[1]-1][0])

    i = 0
    while i != int(data.iloc[spaces[1]+Csol][0]) + 1:
        capacatance.append(float(data.iloc[spaces[1]+Csol +2+i][0]))
        i += 1

    capacatance1 = [cap, capacatance]

    return capacatance1

def varible_scaler(val_in,var):

    var_scal = []
    Nvar = len(var.iloc[:][0])

    j = 0
    while j != Nvar:
        if var.iloc[j][3] == 0:  # Non-logaritmic scalable constant

            if val_in[j] > 0:
                var_scal.append((var.iloc[j][2] - var.iloc[j][0]) * val_in[j] + var.iloc[j][0])  # Parameter scaler
            else:
                var_scal.append((var.iloc[j][0] - var.iloc[j][1]) * val_in[j] + var.iloc[j][0])  # Parameter scaler
        elif var.iloc[j][3] == 1:  # logaritmic scalable constant
            if val_in[j] > 0:
                var_scal.append( 10 ** ((np.log10(var.iloc[j][0] / var.iloc[j][1])) * val_in[j] + np.log10(
                    var.iloc[j][0])))  # Parameter scaler
            else:
                var_scal.append(10 ** ((np.log10(var.iloc[j][2] / var.iloc[j][0])) * val_in[j] + np.log10(
                    var.iloc[j][0])))  # Parameter scaler
        else:
            print('please put a 0 or 1 in the right input')
        j += 1

    return var_scal

# Para_combiner_pints(chains)
def Para_combiner_pints(chains):

    Ncores = len(chains)

    # Nvar and Nchain arn't needed
    dist = np.zeros((Ncores * len(chains[0]), len(chains[0][0])))  # Nchain = resultM[0], Nvar = resultM[0][0]

    # combines the acceptance/MCMC rate
    i = 0
    while i != Ncores:
        dist[i::Ncores, :] = np.array(chains[i][:])
        i += 1

    return dist

# extracts the PINTS settings from the block data
def PINTS_Varibles(settings):

    # extracts header settings
    header = settings.iloc[0][0]
    header = header.split(' ')

    # clean up the header file
    y = []
    for x in header:
        x = x.strip()
        x = x.strip('\t')
        if x != '':
            y.append(x)

    # cleans up the header varible
    x = []
    i = 0
    while i != len(header):
        if header[i] == '\t':
            pass
        else:
            x.append(header[i])
        i += 1
    header = x

    # method selection input
    if header[2] == 'CMAES': # and (header[1] == 'TCDS' or header[1] == 'HEPWsigma' header[1] == 'HarmPerFit' or header[1] == 'Baye_HarmPerFit'):
        var, bandwidth, harm_weights, op_settings, datatype,\
            scalvar, funcvar, truntime = STAND_CMAES_TOTCURR_settings(settings,header)
    elif header[2] == 'ADMCMC': # and (header[1] == 'TCDS' or header[1] == 'HEPWsigma' or header[1] == 'Baye_HarmPerFit'):
        var, bandwidth, harm_weights, op_settings, datatype, \
        scalvar, funcvar, truntime = PINT_ADMCMC_TOTCURR_settings(settings)
    else:
        print('incorrect Header parameters or method not avalable')

    return header, var, bandwidth, harm_weights, op_settings, datatype, scalvar, funcvar, truntime   # Modified CMA VALUES var

# Reads experimental input files of the form MECSim/FTACV
def output_reader_FTMC(name):

    results = pd.read_csv(name, delim_whitespace=True, skiprows=19, names=["v", "i", "t"])
    results = results.values  # changes results fro df object to np array

    Exp_t = results[1, 2]
    Extime = results[-1, 2]
    curr = results[:, 1]

    return curr, Exp_t, Extime

# splits exp data depending on data type
def Exp_data_spliter(Exp_data, datatype,freq, bandwidth,spaces,header,op_settings):
    exp = Exp_data.values   # converts to list
    curr_col = []

    if datatype == 0 or datatype == 1:  # MECsim/FTACV simulation data (of form {v,i,t})
        # extracts Exp file names
        i = 0
        while i != len(Exp_data.iloc[:][0]):

            curr, Exp_t,Extime = output_reader_FTMC(exp[i,0])
            curr_col.append(np.array(curr))
            i += 1

    elif datatype == 2:  # CHI data type (of form {v,i}) (NEED SMETHING TO get time)
        Exp_data = Exp_data[0].str.split(', ', expand=True).astype('float')

        # sets up for usable information
        curr = np.array(Exp_data.values)
        curr_col.append(np.array(curr[:, 1]))
        Exp_t = None  # (CHI IS DC ONLY SO TIME DOMAIN DOEN't matter)
        Extime = None

    else:
        print('need a data type to compair to')


    Np = len(curr_col[0])
    deci = int(Np/op_settings[2])

    # exception for output reader
    if len(Exp_data.iloc[:][0]) == 1:
        if header[1] == 'TCDS':
            curr = curr_col[0][::deci]
            sigma = False
        elif header[1] == 'FTC' or header[1] == 'Log10FTC':
            curr = curr_col[0][:]
            sigma = False

        elif header[1] == 'HarmPerFit':
            curr = MLsp.harm_gen(curr_col[0], Exp_t, freq, bandwidth, spaces)
            curr = curr[:,::deci]
            sigma = False
        else:
            print('incorrect header parameters used, in loading experimental files')
            exit()
    else:

        if header[1] == 'HarmPerFit':   # for this no sigma in needed as just
            curr, sigma = Volt_Sigma(curr_col, Exp_t, freq, bandwidth,spaces,deci)
            sigma = False
        #elif header[1] == 'Baye_HarmPerFit':
        elif header[1] == 'TCDS':
            print('This is for emperical analysis (Warning: NOT FOR CALCULATIONS)')
            curr = curr_col[0][::deci]
            sigma = False
        #    curr, sigma = Volt_Sigma(curr_col, Exp_t, freq, bandwidth, spaces, deci)
        elif header[1] == 'HEPWsigma' or header[1] == 'Baye_HarmPerFit':
            curr, sigma = Volt_Sigma_Pwise(curr_col, Exp_t, freq, bandwidth, spaces,deci)

    # correction to Exp_t for time steps
    Exp_t *= deci

    return curr, Exp_t,Extime, sigma

# truncation method for Simulated total current method
def tot_current_truncation(Excurr, Scurr,**kwargs):

    truncation = kwargs.get('Nsimdeci')
    op_settings = kwargs.get('op_settings')

    # gets decimication required for simulation
    deci = int(len(Scurr) / op_settings[2])

    # Truncates Simulation harmonics at time interval
    Scurr = Scurr[::deci]
    Scurr = Scurr[truncation[0]:truncation[1]]

    if Scurr.shape == Excurr.shape:
        pass
    else:
        print('Experimental current where not of same dimensions as simulated current')
        exit()

    return  Scurr

# Truncation and preprossing for simulated harmonics for a Harmonic envolope method
def Harm_envolop_truncation(EXHarmonics, Scurr,**kwargs):
    # extracts kwargs
    truncation = kwargs.get('Nsimdeci')
    MECSettings = kwargs.get('data')
    spaces = kwargs.get('spaces')
    bandwidth = kwargs.get('bandwidth')
    AC_freq = kwargs.get('AC_freq')
    op_settings = kwargs.get('op_settings')

    dt = abs((2 * (float(MECSettings.iloc[3]) - float(MECSettings.iloc[2])) / float(MECSettings.iloc[5])) / 2 ** int(
        MECSettings.iloc[6]))

    # Extracts Harmonics from simulation
    hil_store = MLsp.harm_gen(Scurr, dt, AC_freq, bandwidth, spaces)

    # gets decimication required for simulation
    deci = int(len(Scurr) / op_settings[2])

    # Truncates Simulation harmonics at time interval
    hil_store = hil_store[:, ::deci]
    hil_store = hil_store[:, truncation[0]:truncation[1]]

    Nharm = len(bandwidth[0]) + 1
    Np = truncation[1] - truncation[0]

    return hil_store, Nharm, Np

########################################################################################################################
############################################### METHOD FITTING FUNCTIONS ###############################################
########################################################################################################################
""" All the functions for fitting simulated Curr to ExCurr"""
# decimate simulation current in the comparison

# Total Current difference squared
def TotCurr_diffsquare(Excurr, Scurr,**kwargs):

    # truncates total current from simulation
    Scurr = tot_current_truncation(Excurr, Scurr,**kwargs)

    Loss = sum((Excurr - Scurr)**2)

    return Loss

def BayeTotCurr_diffsquare(Excurr,Scurr,sigma,**kwargs):

    Scurr = tot_current_truncation(Excurr, Scurr, **kwargs)

    N = len(Scurr)

    diff = sum((Excurr - Scurr) ** 2)

    Loss = -N*np.log(sigma) - diff/(2*sigma**2)

    return Loss

# Fourier Transform difference squared
def FTC_diffsquare(Excurr, Scurr,**kwargs):
    Nsimdeci = kwargs.get('Nsimdeci')

    # normalization
    Scurr = np.fft.rfft(Scurr)/len(Scurr)
    Scurr = Scurr[Nsimdeci[0]:Nsimdeci[1]]

    # need the absolute as functions are imaginary
    Loss = sum(np.abs((Excurr - Scurr))**2)

    return Loss

# Log10 Fourier Transform difference squared
def Log10FTC_diffsquare(Excurr, Scurr,**kwargs):

    Nsimdeci = kwargs.get('Nsimdeci')

    # truncates total current from simulation
    Scurr = np.fft.rfft(Scurr)
    Scurr = np.abs(Scurr)/len(Scurr)

    Current = np.empty([0])#Scurr[Nsimdeci[0][0]:Nsimdeci[0][1]]
    #del Nsimdeci[0]

    for X in Nsimdeci:
        Current = np.concatenate((Current, Scurr[X[0]:X[1]]))

    Scurr = np.log10(Current)

    # note that for this case the Excurr is calculated and turned to base ten early in calculations
    Loss = sum(abs((Excurr.real - Scurr.real)**2))

    return Loss

# Harmonic Perctage best fit
def HE_percentagefit(EXHarmonics, Scurr,**kwargs):

    harm_weights = kwargs.get('harm_weights')

    # extracts processes and truncates simulated harmonics to be of same dimensions as EXharm
    hil_store, Nharm, Np = Harm_envolop_truncation(EXHarmonics, Scurr, **kwargs)

    y = 0
    for j in range(Nharm):
        x = sum((EXHarmonics[j,:] - hil_store[j,:])**2)/sum(EXHarmonics[j,:]**2)
        y += harm_weights[j]*np.sqrt(x)

    Loss = y/Nharm

    return Loss

# Harmonic Perctage best fit
def Baye_HarmPerFit(EXHarmonics, Scurr,Exsigma,**kwargs):

    harm_weights = kwargs.get('harm_weights')

    # extracts processes and truncates simulated harmonics to be of same dimensions as EXharm
    hil_store, Nharm, Np = Harm_envolop_truncation(EXHarmonics, Scurr, **kwargs)

    conc = -np.log(2*np.pi)*0.5
    z = 0
    for j in range(Nharm):
        x = sum((EXHarmonics[j,:] - hil_store[j,:])**2)/sum(EXHarmonics[j,:]**2)
        y = conc - np.log(Exsigma[j]) - (np.sqrt(x)/(2 * Exsigma[j] ** 2))
        z += harm_weights[j]*y      # using the WLB method (Newton and Raftery 1994)

    Loss = z

    return Loss


# Harmonic Envolope Point wise calculated sigma value
def HEPWsig_Bdiffsquare(EXHarmonics, Scurr, Exsigma,**kwargs):
    harm_weights = kwargs.get('harm_weights')
    # extracts processes and truncates simulated harmonics to be of same dimensions as EXharm
    hil_store, Nharm, Np = Harm_envolop_truncation(EXHarmonics, Scurr, **kwargs)

    error = []
    conc = -np.log(2*np.pi)*0.5      # log likilyhood constant

    x = 0
    #print('need to check if this is a fair equalizing parameter for log likilyhood')
    for j in range(Nharm):

        y = sum((conc - np.log(Exsigma[j,:]) - (((EXHarmonics[j,:] - hil_store[j,:])**2)/(2*Exsigma[j,:]**2))))
        x += harm_weights[j]*y      # using the WLB method (Newton and Raftery 1994)

    Loss = x

    return Loss

########################################################################################################################
###############################################      REST OF THE CRAP    ###############################################
########################################################################################################################



# Modal handler for MECSim in Pints
class CMAES_Model(pints.ForwardModel):

    def __int__(self,**kwargs):

        self.kwargs = kwargs

    def simulate(self, val_in, time):
        # Run a simulation with the given parameters for the
        # given times
        # and return the simulated value

        var = self.kwargs.get('var')
        val_in = varible_scaler(val_in,var)
        Scurr = MLsp.Iterative_MECSim_Curr(val_in, **self.kwargs)

        return Scurr


    def n_parameters(self):

        var = self.kwargs.get('var')

        return len(var[:][0])


# handles noise for the MCMC so that it is not passed to the MECSim
def Noise_handler(MEC_set):

    var = MEC_set.pop('var')    # takes var out of dic

    i = 0
    N = len(var.iloc[:][0])
    noisebool = False
    holder = []
    # checks to see if noise para is there
    while i != N:
        if var.iloc[i][4] == 13:    # noise para
            noise = var.iloc[i][:].values
            noisebool = True
        else:
            holder.append(var.iloc[i][:].values)
        i += 1

    if noisebool:
        var = pd.DataFrame(holder)
    else:
        noise = None

    MEC_set.update({'var':var})

    return var, MEC_set, noise, noisebool

# function to take in the priors and outputs a sigma function
def sigma_set(priors):

    # preallocation
    N = len(priors[0])  # gets the # of varibles
    covar = np.zeros((N, N))

    i = 0
    while i != N:

        covar[i, i] = ((priors[1][i] - priors[0][i])**2)/ 12  # varience of a uniform distrabution

        i += 1

    return covar

def CMA_var_handler(Iterative_MECSim0):  # AC is needed so args can be used in the wrapper
    def CMA_Var(*args, **kwargs):  # """THIS IS GOING TO BE A PAIN # data,var,spaces,Ex_data """

        DCAC_method = kwargs.get('DCAC_method')
        var = kwargs.get('var')

        # Will need to put an if loop to not do a parallel
        Ncore = DCAC_method[2]  # gets the number of usable cores
        listin = []  # creates a blank list
        i = 0
        for input_val in args:
            listin.append(input_val)  # adds the modified array to new list
            i += 1

        Nvar = len(listin[0])  # number of varibles being modified
        """parallel starts here FUCK FUCK FUCK"""
        Npop = len(args)  # gets population iteration size
        i = 0  # CMA input counter
        while i != Npop:  # this sets the CMA MODEL INPUT forif DCAC_method == 1: # AC method
            j = 0
            while j != Nvar:
                if var.iloc[j][3] == 0:  # Non-logaritmic scalable constant
                    if listin[i][j] > 0:
                        listin[i][j] = (var.iloc[j][2] - var.iloc[j][0]) * listin[i][j] + var.iloc[j][0]  # Parameter scaler
                    else:
                        listin[i][j] = (var.iloc[j][0] - var.iloc[j][1]) * listin[i][j] + var.iloc[j][0]  # Parameter scaler
                elif var.iloc[j][3] == 1:  # logaritmic scalable constant
                    if listin[i][j] > 0:
                        listin[i][j] = 10 ** ((np.log10(var.iloc[j][0] / var.iloc[j][1])) * listin[i][j] + np.log10(
                            var.iloc[j][0]))  # Parameter scaler
                    else:
                        listin[i][j] = 10 ** ((np.log10(var.iloc[j][2] / var.iloc[j][0])) * listin[i][j] + np.log10(
                            var.iloc[j][0]))  # Parameter scaler
                else:
                    print('please put a 0 or 1 in the right input')
                j += 1
            i += 1

        # multicore processsing

        # Multiprocessing call around Iterative MECSim
        with Pool(processes=Ncore) as p:
            multi = p.map(Iterative_MECSim0, listin)

        # prints the multivarible output in order
        val_out = []
        [val_out.append(s) for s in multi]  # CMA-ES output is a list of NP-arrays

        # iter_logger(listin, out)  # append in loop then print to txt at c point redesign

        return val_out
    return CMA_Var

# scaler for the cma between predefined varibles
def CMA_varible_scaler(listin,**kwargs):

    var = kwargs.get('var')

    Nvar = len(listin[0])  # number of varibles being modified
    """parallel starts here FUCK FUCK FUCK"""
    Npop = len(listin)  # gets population iteration size
    i = 0  # CMA input counter
    while i != Npop:  # this sets the CMA MODEL INPUT forif DCAC_method == 1: # AC method
        j = 0
        while j != Nvar:
            if var.iloc[j][3] == 0:  # Non-logaritmic scalable constant
                if listin[i][j] > 0:
                    listin[i][j] = (var.iloc[j][2] - var.iloc[j][0]) * listin[i][j] + var.iloc[j][
                        0]  # Parameter scaler
                else:
                    listin[i][j] = (var.iloc[j][0] - var.iloc[j][1]) * listin[i][j] + var.iloc[j][
                        0]  # Parameter scaler
            elif var.iloc[j][3] == 1:  # logaritmic scalable constant
                if listin[i][j] < 0:
                    listin[i][j] = 10 ** ((np.log10(var.iloc[j][0] / var.iloc[j][1])) * listin[i][j] + np.log10(
                        var.iloc[j][0]))  # Parameter scaler
                else:
                    listin[i][j] = 10 ** ((np.log10(var.iloc[j][2] / var.iloc[j][0])) * listin[i][j] + np.log10(
                        var.iloc[j][0]))  # Parameter scaler
            else:
                print('please put a 0 or 1 in the right input')
            j += 1
        i += 1

    return listin


###################################### Experimental sigma calculator ###################################################

# superived percentage calculator per harmonic
def Volt_Sigma(Curr,Exp_t, freq, bandwidth,spaces,deci):

    Ndata = len(Curr)

    Nharm = len(bandwidth[0,:]) + 1
    Nac = spaces[4] # for now only do one

    # gets the length of the harmonics
    Np = int(len(Curr[0]) / deci)

    # preallocate the harmonics in a numpy files
    harmdic = {}
    i = 0
    while i != Nharm:
        s = '%i' % (i)  # sets harmonic

        harmdic[s] = np.zeros((Ndata, Np))

        i += 1

    # generate each harmonic from the current file then adds to a dictionary
    i = 0
    while i != Ndata:

        # gets harmonics
        harmstore = MLsp.harm_gen(Curr[i], Exp_t, freq, bandwidth, spaces)

        j = 0
        while j != Nharm:
            sn = '%i' % (j)
            harmcol = harmdic.pop(sn)  # this is really unefficent but im lazy

            # add harmonic from harmstore to dictionary store
            harmcol[i, :] = harmstore[j, 0::deci]

            # add back to dictionary
            harmdic[sn] = harmcol

            j += 1

        i += 1

    # pre-allocates mean and std
    harmmean = np.zeros((Nharm, Np))
    harmstd = np.zeros((Nharm, Np))

    # calculate mean and std of each point
    i = 0
    while i != Nharm:

        sn = '%i' % (i)
        harmcol = harmdic.get(sn)  # this is really unefficent but im lazy

        # go through the data points
        for j in range(Np):
            harmmean[i, j] = np.mean(harmcol[:, j])
            harmstd[i, j] = np.std(harmcol[:, j],ddof=1)    # / np.sqrt(Ndata)  # Calculates the std of the mean of the current value.

        i += 1

    # propigate error throughout
    sumerr = np.zeros(Nharm)
    i = 0
    while i != Nharm:

        x = sum((harmstd[i, :] * harmmean[i, :]) ** 2)
        sumerr[i] = 2 * np.sqrt(x) / (np.dot(harmmean[i, :], harmmean[i, :]))  # normalizes the std

        i += 1

    return harmmean, np.sqrt(sumerr)

def Harm_sigma_normal(harmmean,harmstd):

    Nharm = harmmean.shape[0]

    # propigate error throughout
    sumerr = np.zeros(Nharm)
    i = 0
    while i != Nharm:
        x = sum((harmstd[i, :] * harmmean[i, :]) ** 2)
        sumerr[i] = 2 * np.sqrt(x) / (np.dot(harmmean[i, :], harmmean[i, :]))  # normalizes the std

        i += 1

    return np.sqrt(sumerr)

# point wise sigma for eac point in the harmonics
def Volt_Sigma_Pwise(Curr,Exp_t, freq, bandwidth,spaces,deci):

    Ndata = len(Curr)
    Nharm = len(bandwidth[0,:]) + 1
    Nac = spaces[4] # for now only do one

    # gets the length of the harmonics
    Np = int(len(Curr[0]) / deci)

    # preallocate the harmonics in a numpy files
    harmdic = {}
    i = 0
    while i != Nharm:
        s = '%i' % (i)  # sets harmonic

        harmdic[s] = np.zeros((Ndata, Np))

        i += 1

    # generate each harmonic from the current file then adds to a dictionary
    i = 0
    while i != Ndata:

        # gets harmonics
        harmstore = MLsp.harm_gen(Curr[i], Exp_t, freq, bandwidth, spaces)

        j = 0
        while j != Nharm:
            sn = '%i' % (j)
            harmcol = harmdic.pop(sn)  # this is really unefficent but im lazy

            # add harmonic from harmstore to dictionary store
            harmcol[i, :] = harmstore[j, 0::deci]

            # add back to dictionary
            harmdic[sn] = harmcol

            j += 1

        i += 1

    # pre-allocates mean and std
    harmmean = np.zeros((Nharm, Np))
    harmstd = np.zeros((Nharm, Np))

    # calculate mean and std of each point
    i = 0
    while i != Nharm:
        sn = '%i' % (i)
        harmcol = harmdic.get(sn)  # this is really unefficent but im lazy

        # go through the data points
        j = 0
        for j in range(Np):
            harmmean[i, j] = np.mean(harmcol[:, j])
            harmstd[i, j] = np.std(harmcol[:, j],ddof=1)    #/ np.sqrt(Ndata)#Calculates the std of the mean of the curr value.

            j += 1

        i += 1

    return harmmean, harmstd


# Harmonic Perctage best fit
def HARM_percentagefit(EXHarmonics, Scurr,**kwargs):

    # extracts processes and truncates simulated harmonics to be of same dimensions as EXharm
    hil_store, Nharm, Np = Harm_envolop_truncation(EXHarmonics, Scurr, **kwargs)


    j = 0
    y = []
    while j != Nharm:
        x = sum((EXHarmonics[j,:] - hil_store[j,:])**2)/sum(EXHarmonics[j,:]**2)
        y.append(np.sqrt(x))
        j += 1

    return y
