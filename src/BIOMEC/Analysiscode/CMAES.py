# sub module stuff for CMA-ES using a Nickol Hansens python interface
# Aurthor: Luke Gundry
# Date: 2/9/19

import pandas as pd
import numpy as np
import cma
import pints
from multiprocessing import Pool
from ML_signal_processing import Iterative_MECSim, Iterative_MECSim_Curr
from Script_generator import iter_logger
from MCMC_modules import *
import timeseries_modules as tpseries
import matplotlib.pyplot as plt
from copy import deepcopy

#Stadarnd class for case where non pints modal
class standard_CMAES_Model(object):
    # preallocation of self
    kwargs = {}
    EXcurr = None
    Method = None
    EXSigma = None

    def __init__(self,*args, **kwargs):

        # args = [Excurr, Method]
        self.kwargs = kwargs

        #need someshit to specify which METHOD TO USE

        # defines the EXperimental current
        self.EXcurr = args[0]
        self.Method = args[1]
        self.Bayesfit = False       # for the log likily hoodversions of the univariate versions

        if self.Method == 'HEPWsigma' or self.Method == 'Baye_HarmPerFit':
            self.EXSigma = kwargs.pop('Exsigma')
            self.filters = kwargs.get('filters')

        elif self.Method == 'TCDS' or self.Method == 'FTC' or self.Method == 'Log10FTC':
            var = kwargs.get('var')
            Nvar = len(var.iloc[:][0])
            for i in range(Nvar):
                if var.iloc[i][4] == 13:
                    self.EXSigma = i
                    self.Bayesfit = True

        elif self.Method == 'Bayes_ExpHarmPerFit':
            self.EXPperrErr = kwargs.pop("EXPperrErr")
            self.EXSigma = kwargs.pop('Exsigma')
            self.filters = kwargs.get('filters')


    def simulate(self, *args):

        Scurr = Iterative_MECSim_Curr(args[0], **self.kwargs)
        # fitting function
        if self.Method == 'TCDS':

            if self.Bayesfit:   # use univariate loglikilyhood
                Loss = tpseries.TotCurr_diffsquare(self.EXcurr, Scurr,**self.kwargs)
                conc = -len(self.EXcurr) * np.log(2 * np.pi)*0.5
                Loss = conc - len(self.EXcurr)*np.log(args[0][self.EXSigma]) - Loss/(2*(args[0][self.EXSigma]**2))
                Loss = abs(Loss)
            else:
                Loss = tpseries.TotCurr_diffsquare(self.EXcurr, Scurr, **self.kwargs)

        #elif self.Method == 'Baye_TCDS':

        #   Loss = tpseries.BayeTotCurr_diffsquare(self.EXcurr, Scurr, self.EXSigma,**self.kwargs)
        #    Loss = abs(Loss)

        elif self.Method == 'FTC':
            if self.Bayesfit:   # use univariate loglikilyhood
                Loss = tpseries.FTC_diffsquare(self.EXcurr, Scurr,**self.kwargs)
                conc = -len(self.EXcurr)*np.log(2*np.pi)/2
                Loss = conc -len(self.EXcurr)*np.log(args[0][self.EXSigma]) - Loss/(2*(args[0][self.EXSigma]**2))
                Loss = abs(Loss)
            else:
                Loss = tpseries.FTC_diffsquare(self.EXcurr, Scurr,**self.kwargs)

        elif self.Method == 'Log10FTC':
            if self.Bayesfit:   # use univariate loglikilyhood
                Loss = tpseries.Log10FTC_diffsquare(self.EXcurr, Scurr ,**self.kwargs)
                Loss = -len(self.EXcurr)*np.log(args[0][self.EXSigma]) - Loss/(2*(args[0][self.EXSigma]**2))
                Loss = abs(Loss)
            else:
                Loss = tpseries.Log10FTC_diffsquare(self.EXcurr, Scurr, **self.kwargs)

        elif self.Method == 'HEPWsigma':

            Loss = tpseries.HEPWsig_Bdiffsquare(self.EXcurr, Scurr,self.EXSigma,**self.kwargs)
            Loss = abs(Loss)    # minimization for CMA-ES algorithim

        elif self.Method == 'HarmPerFit':
            # Calculates percentage best fit
            Loss = tpseries.HE_percentagefit(self.EXcurr, Scurr, **self.kwargs)

        elif self.Method == 'Baye_HarmPerFit':
            # Calculates percentage best fit
            # harmonic percentage traking hasn't been set up for CMA-ES yet
            Loss, harmperfit = tpseries.Baye_HarmPerFit(self.EXcurr, Scurr, self.EXSigma,**self.kwargs)
            Loss = abs(Loss)    # minimization for CMA-ES algorithim

        elif self.Method == 'Bayes_ExpHarmPerFit':
            # calculates harmonic percentage error based on simulation
            Loss, harmperfit = tpseries.Baye_ExpHarmPerFit(self.EXcurr, Scurr, self.EXSigma, self.EXPperrErr,**self.kwargs)
            Loss = abs(Loss)

        else:
            print("incorrect fitting function header used")

        return Loss
	
	

# Modal handler for MECSim in Pints
class CMAES_Model(pints.ForwardModel):

    def __int__(self,**kwargs):

        self.kwargs = kwargs

    def simulate(self, val_in, time):
        # Run a simulation with the given parameters for the
        # given times
        # and return the simulated value

        var = self.kwargs.get('var')
        val_in = tpseries.varible_scaler(val_in,var)
        Scurr = Iterative_MECSim_Curr(val_in, **self.kwargs)

        return Scurr


    def n_parameters(self):

        var = self.kwargs.get('var')

        return len(var[:][0])

    def n_outputs(self):

            # Return the dimension of the parameter vector
        return 1


################################################################################################
######################################   CMAES_TOTCURR_METHOD  #################################
################################################################################################
# input reader for this method
def STAND_CMAES_TOTCURR_settings(settings,header):

    nvar = int(settings.iloc[1][0])

    # stores the varibles
    var = settings[2:nvar + 2]
    var = var[0].str.split(',', expand=True).astype('float')
    n = settings[0].shape
    n = int(n[0])
    nscal = int(settings.iloc[nvar + 2][0])  # gets number of scaling parameters
    Nfunc = int(settings.iloc[nvar + nscal + 3][0])

    # Gets the scaling
    scalvarpre = settings[nvar + 2:nvar + 3 + nscal]  # gets the num of scaling to
    # turns scalvarpre into a list of values
    scalvar = [[int(scalvarpre.iloc[0][0])]]
    i = 0
    x = scalvarpre[0].str.split(',', expand=True).astype('float')
    x = x.values

    while i != scalvar[0][0]:
        scalvar.append(x[i + 1, :])
        i += 1

    # Sets up to read functional parameters
    if Nfunc == 0:
        funcvar = [[0]]
    else:
        funcvarpre = settings[nvar + 3 + nscal:nvar + 4 + nscal + Nfunc]  # gets the num of scaling to
        # turns scalvarpre into a list of values
        funcvar = [[int(funcvarpre.iloc[0][0])]]
        i = 0
        x = funcvarpre[0].str.split(',', expand=True).astype('float')
        x = x.values

        while i != Nfunc:
            funcvar.append(x[i + 1, :])
            i += 1

    Nlist = 7
    nband = int((n - Nlist - nvar - 4 - nscal - Nfunc) / 2)  # counts the number that isn't there
    Nwinfunc = nvar + 1 + 2 + nscal + Nfunc + 1  # counter cause im lazy

    # extracts the harmonic windowing information
    Harmonicwindowing = settings[Nwinfunc:Nwinfunc+1]
    Harmonicwindowing = Harmonicwindowing.iloc[0][0].split(',')
    Harmonicwindowing = [float(array) for array in Harmonicwindowing] # gets rid of \t and converts  all values to floats

    filler = []
    if Harmonicwindowing[0] == 0: # This is for the case of square windowing
        filler.append("squarewindow")
    elif Harmonicwindowing[0] == 1: # This is for the convolutional guassin function from paper ~ 2014 mash?
        filler.append("Conv_guassian")
        filler.append(Harmonicwindowing[1]) # guassian std for function
    else:
        print("WARNING INCORRECT WINDOWING FUNCTION USED, HOPEFULLY DOESN'T KILL ANYTHING")
    Harmonicwindowing = filler

    # Seperates the CMA input data
    Nbarndersnatch = nvar + 1 + 2 + nscal + Nfunc + 2  # counter cause im lazy
    bandwidth = settings[Nbarndersnatch:Nbarndersnatch + nband]
    harm_weights = settings[Nbarndersnatch + nband:Nbarndersnatch + 2 * nband]

    # Changes the seperated input values into
    bandwidth = bandwidth[0].str.split(',', expand=True).astype('float')
    bandwidth = bandwidth.values  # changes from a df to np.array
    harm_weights = harm_weights[0].str.split(',', expand=True).astype('float')
    harm_weights = harm_weights.values

    # Optimization settings
    Np = 2 * nband + Nbarndersnatch - 1  # starting point

    # truncation settings for time series analysis
    truntime = settings.iloc[Np + 1][0]
    truntime = truntime.split(',')

    x = [float(truntime[0]),truntime[1].strip()]        # strip is to remove tabs
    if x[1] == 'MAX': # allows for name input
        pass
    else:
        x[1] = float(x[1])
    truntime = x

    op_settings = [0, 0, 0, 0, 0]  # prealocation
    op_settings[0] = int(settings.iloc[Np + 2][0])  # This one setts the optimization method
    op_settings[1] = int(settings.iloc[Np + 4][0])  # this one does the num cores
    datatype = int(settings.iloc[Np + 3][0])  # Gets the type of data were compairing against
    op_settings[2] = 2 ** float(settings.iloc[Np + 5][0])  # sets the number of datapoints to compair in the current
    op_settings[3] = float(settings.iloc[Np + 6][0])  # sets the tolx value
    op_settings[4] = float(settings.iloc[Np + 7][0])  # sets the initial sigma

    return var, Harmonicwindowing, bandwidth, harm_weights, op_settings, datatype, scalvar, funcvar, truntime

# STAND CMA-ES TOTCURR Runner
def STAND_CMAES_TOTCURR(var, op_settings, Method,MEC_set, Excurr):
    #CMA_MEC(scalvar, **kwargs):

    x = 1
    N = np.random.randint(1000, 9999)

    Dim = var.shape[0]  # Extracts the number of varibles being used
    popsize = int(3*np.log(Dim))  # calculates number of function evals per iteration (popsize)
    options = {'seed': N, 'tolx': op_settings[3], 'bounds': [-x, x],
               'tolstagnation': int(100 + 100 * Dim ** 1.5 / (popsize * 2)),
               'maxiter': 1000}  # Keyword assignment for CMA, 'tolfun':100 'tolx':0.05

    # Allocates wrapper  around iterative MECSim function
    model = standard_CMAES_Model(*[Excurr, Method], **MEC_set)
    run = model.simulate

    # need the paralleliztion
    DCAC_method = MEC_set.get('DCAC_method')
    # Will need to put an if loop to not do a parallel
    Ncore = DCAC_method[2]  # gets the number of usable cores

    # initial std (Needs to be more personalized) output would be a list so just needs modification from there
    sigma = op_settings[4]  # 33 probs better but sebs

    # logger = cma.CMAESDataLogger()      # records whats going on iterativly
    es = cma.CMAEvolutionStrategy(Dim * [0], sigma, options)  # optimizes to a minimal point
    # logger = cma.CMADataLogger().register(es)
    "Might need an if loop for single core systems"
    # func = CMA_Iterative_MECSim(*num,**MEC_set)
    # with cma.fitness_transformations.EvalParallel(DCAC_method[2]*2) as eval_all:
    logtot_X = []
    logtot_Y = []
    from operator import methodcaller
    while not es.stop():
        # the multiprocessing needs to go in here
        X = es.ask()

        # scales functions of bounds -1,1 to dimensional inputs
        #deepcopy is to get around immutability of the arrays
        listin = tpseries.CMA_varible_scaler(deepcopy(X),  **MEC_set)

        # input output list to do parallel in wrapper
        with Pool(processes=Ncore) as p:

            FIT = p.map(run, listin)

        es.tell(X, FIT)  # Does the CMA-ES fitting

        # makes a output file
        logtot_X += X
        logtot_Y += FIT

    # combines logtot X and Y and dimensionalizes the values
    logtot = []
    Ndim = len(logtot_Y)
    i = 0
    while i != Ndim:
        x = tpseries.varible_scaler(logtot_X[i], var)
        x = np.append(x,logtot_Y[i])
        logtot.append(x)
        i += 1

    results = es.result

    return results, logtot


# Final translation
def CMA_output(res,**kwargs):

    var = kwargs.get('var')
    DCAC_method = kwargs.get('DCAC_method')
    ExP = kwargs.get('Exp_data')
    Nsimdeci = kwargs.get('Nsimdeci')

    N = var.shape[0]
    var_out = list(np.zeros(N))

    res1 = res[0]  # reassigns to best solution [5] better with noise

    j = 0
    while j != N:
        if var.iloc[j][3] == 0:  # Non-logaritmic scalable constant
            if res1[j] > 0:
                var_out[j] = (var.iloc[j][2] - var.iloc[j][0]) * res1[j] + var.iloc[j][0]  # Parameter scaler
            else:
                var_out[j] = (var.iloc[j][0] - var.iloc[j][1]) * res1[j] + var.iloc[j][0]  # Parameter scaler

            j += 1
            # function added in to adjust logarithmic values for cma-es
        elif var.iloc[j][3] == 1:  # logaritmic scalable constant

            if res1[j] < 0:
                var_out[j] = 10 ** ((np.log10(var.iloc[j][0] / var.iloc[j][1])) * res1[j] + np.log10(
                    var.iloc[j][0]))  # Parameter scaler
            else:
                var_out[j] = 10 ** ((np.log10(var.iloc[j][2] / var.iloc[j][0])) * res1[j] + np.log10(
                    var.iloc[j][0]))  # Parameter scaler

            j += 1
        else:
            print('please put a 0 or 1 in the right input')

    # here so we can allocate mean out to run on a processor

    EXcurr = kwargs.get('Exp_data')
    if DCAC_method[0] == 1:  # AC method

        harm_weights = kwargs.get('harm_weights')  # need that flexability function
        bandwidth = kwargs.get('bandwidth')
        Scurr = Iterative_MECSim_Curr(var_out, **kwargs)
        # compares experimental harmonic EXHarmonics
        Loss = tpseries.HARM_percentagefit(ExP, Scurr,**kwargs)

    elif len(Nsimdeci) != 2:    # Fourier Transform method

        Scurr = Iterative_MECSim_Curr(var_out, **kwargs)
        # compares experimental Total current
        Loss = tpseries.Log10FTC_diffsquare(ExP, Scurr,**kwargs)

    else:  # DC method

        Scurr = Iterative_MECSim_Curr(var_out, **kwargs)
        # compares experimental Total current
        Loss = tpseries.TotCurr_diffsquare(ExP, Scurr,**kwargs)

    return var_out, Loss


# output writter for this method
def PINT_CMAES_TOTCURR_output(filename,t2, space_holder, var_out, res, mean_var_out, DCAC_method,loss,header,MEC_set):

    # t2,space_holder, var_out,Perr,res,c
    # For AC the main difference is that Perr is a matrix

    f = open(filename, "w")  # creates a blank .inp to be written

    # writes the input
    f.write('Input file name: %s\n' % (filename))
    f.write('Logic and fitting method used: %s, %s\n'% (header[1], header[2]))
    # inserts time
    f.write('Completetion time (min): %f\n\n' % (t2))

    # var_out
    f.write('Optimized paramters;{val,code,repeat}\n')
    N = len(var_out)
    x = [0, 0, 0]
    i = 0
    while i != N:
        x[0] = format_e(var_out[i])  # var out is a string
        x[1] = str(int(space_holder[i, 0]))
        x[2] = str(int(space_holder[i, 1]))
        y = ",".join(x)
        f.write('%s\n' % (y))
        i += 1
    f.write('\n')

    f.write('CMA-ES unscaled best fit\n')
    N = len(var_out)
    i = 0
    while i != N:

        y = ",".join(x)
        f.write('%f\n' % (res[0][i]))
        i += 1
    f.write('\n')

    # res
    f.write('Misc CMA_ES output\n')
    x = res[1]
    x = format_e(x)
    f.write('CMA best Fit: %s\n\n' % (x))

    # prints the loss of mean
    if DCAC_method[0] == 1: # prints the AC fit of the harmonics
        f.write('Harmonic percentage best fit:\n')
        for i in range(len(loss)-1):
            f.write(format_e(loss[i])+', ')
        f.write(format_e(loss[-1]))
        f.write('\n\n')

    # mean_var_out = str(mean_var_out)
    f.write('CMA-ES Mean parametervalues (mean,std):\n')
    i = 0
    while i != N:
        # formats values into something we can use
        x1 = mean_var_out[0][i]
        x2 = mean_var_out[1][i]
        m = format_e(x1)
        s = format_e(x2)

        f.write('Var %d: %s,%s\n' % (int(i + 1), m, s))
        i += 1

    f.write('\n')

    # theres the other output values just need to know what they are
    f.write('CMA-ES number of function evaluations: %s\n\n' % (res[3]))
    # CMAstd = res[6] STD cannot be done due to the assymetrical propities of the transform
    f.write('CMA-ES stopping criteria:\n')
    print(res)
    # this has ben removed from
    for i in res[7]:
        f.write('{'+i+' : ')
        f.write('%f}' %res[7][i])

    f.write("\n\nLogic specific outputs\n")
    if header[1] == 'Baye_HarmPerFit':
        sigma = MEC_set.get('Exsigma')
        f.write('harmonic Sigmas\n')
        f.write(np.array2string(sigma))
    f.close()


    return

