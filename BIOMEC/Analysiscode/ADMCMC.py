# sub module stuff for ADMCMC using a Pint interface
# Aurthor: Luke Gundry
# Date: 2/9/19


import pandas as pd
import numpy as np
import pints
import timeseries_modules as tpseries
from multiprocessing import Pool
from ML_signal_processing import Iterative_MECSim, Iterative_MECSim_Curr
from Script_generator import iter_logger
from MCMC_modules import *


# simulation class for MECSim
class ADMCMC_Model(pints.ForwardModel):

    def __int__(self, **kwargs):
        self.kwargs = kwargs
        self.Nsimdeci = kwargs.get('Nsimdeci')
        self.op_settings = kwargs.get('op_settings')

    def simulate(self, val_in, time):
        # Run a simulation with the given parameters for the

        # and return the simulated value
        Scurr = Iterative_MECSim_Curr(val_in, **self.kwargs)

        deci = int(len(Scurr) / self.op_settings[2])

        # Truncates Simulation harmonics at time interval
        Scurr = Scurr[::deci]
        Scurr = Scurr[self.Nsimdeci[0]:self.Nsimdeci[1]]

        return Scurr

    def n_parameters(self):
        var = self.kwargs.get('var')

        return len(var[:][0])

    def parallel(self):
        op_settings = self.kwargs.get('op_settings')

        return op_settings[1]



    def n_outputs(self):
        # Return the dimension of the parameter vector
        return 1

#Stadarnd class for case where non pints modal
class standard_ADMCMC_Model():
    # preallocation of self
    kwargs = {}
    EXcurr = None
    Method = None
    EXSigma = None

    def __init__(self,*args, **kwargs):

        # args = [Excurr, Method]

        #kwargs is the MEC_SET
        self.kwargs = kwargs

        #need someshit to specify which METHOD TO USE

        # defines the EXperimental current
        self.EXcurr = args[0]
        self.Method = args[1]

        # defines if method is a harmonic based method (needed to trac percentage error)
        if self.Method == 'HEPWsigma' or self.Method == 'Baye_HarmPerFit':
            self.EXSigma = kwargs.pop('Exsigma')
            self.Harmonicmethod = True
            self.harmfitstore = []
        elif self.Method == 'Bayes_ExpHarmPerFit':
            self.EXPperrErr = kwargs.pop("EXPperrErr")
            self.EXSigma = kwargs.pop('Exsigma')
            self.Harmonicmethod = True
            self.harmfitstore = []
        else:
            self.Harmonicmethod = False


    # runs the statistical chains used in the inference sampling
    def ADMCMC_chain(self, *args):

        # gets the input arguments from kwargs
        phi0 = self.kwargs.get('phi0')
        covar = self.kwargs.get('covar')
        prior_range = self.kwargs.get('prior_range')
        op_settings = self.kwargs.get('op_settings')

        # sets chain specific settings
        np.random.seed(args[0][1])  # sets the random number seed

        # reset the lengths of the chain with respect to the parallelization
        NMax = int(op_settings[4]/op_settings[7])  # length of overall chain
        Nsamp = int(op_settings[3]/op_settings[7])  # length of initial sampling

        chain, Naccept, Err = self.MCMC_burnin_sample(phi0, covar, prior_range, Nsamp)

        chain, rate, Err = self.adaptive_sampling(chain, Err,covar, prior_range, Nsamp, NMax, Naccept)

        if self.Harmonicmethod:
            return chain, rate, Err, self.harmfitstore
        else:
            return chain, rate, Err

    def simulate(self, *args):
        # Run a simulation with the given parameters for the
        # given times
        # and return the simulated value

        Scurr = Iterative_MECSim_Curr(args[0], **self.kwargs)

        if self.Method == 'HEPWsigma':

            Loss = tpseries.HEPWsig_Bdiffsquare(self.EXcurr, Scurr,self.EXSigma,**self.kwargs)

        elif self.Method == 'Baye_HarmPerFit':

            Loss, harmfitstore = tpseries.Baye_HarmPerFit(self.EXcurr, Scurr, self.EXSigma, **self.kwargs)

            return Loss, harmfitstore

        elif self.Method == 'Bayes_ExpHarmPerFit':
            # calculates harmonic percentage error based on simulation
            Loss, harmfitstore = tpseries.Baye_ExpHarmPerFit(self.EXcurr, Scurr, self.EXSigma, self.EXPperrErr,**self.kwargs)

            return Loss, harmfitstore

        else:
            print("incorrect fitting function header used")

        return Loss

    # the calculation step between the priori and final sampler
    def MCMC_burnin_sample(self, phi0, covar, prior_range, Nsamp):

        # allocates error
        Err = []

        ao = 1
        t = 1  # iteration counter
        N = Nsamp * len(phi0)  # max counter
        phit = [np.array(phi0)]  # sets the initial mean
        Naccept = 0  # acceptance counter

        if self.Harmonicmethod: # goes into harmonic percentage fit
            curr_fit, harmcurrentfit = self.simulate(phi0)  # experimental data gets passed in kwargs
            self.harmfitstore.append(harmcurrentfit)
        else:
            curr_fit = self.simulate(phi0)

        Err.append(curr_fit)

        while t != N + 1:

            phi_trail = np.random.multivariate_normal(phit[-1], ao*covar)
            prob = prior_test(phi_trail, prior_range)  # use this to fitler experimental results with respect to prior

            if prob != 0:  # checks to see if position is possible

                # runs simulation
                if self.Harmonicmethod:  # goes into harmonic percentage fit
                    logtrial, harmtestfit = self.simulate(phi_trail)
                else:
                    logtrial = self.simulate(phi_trail)

                # Probability comparison code
                comp = np.exp(logtrial - curr_fit)
                r = min(1, comp)

                u = np.random.uniform(0, 1)  # generates random uniform distrabution
                if u < r:

                    phit.append(phi_trail)
                    curr_fit = logtrial  # resets the sum
                    Err.append(curr_fit)
                    if self.Harmonicmethod:
                        harmcurrentfit = harmtestfit
                        self.harmfitstore.append(harmcurrentfit)

                    Naccept += 1  # accept rate

                else:
                    phit.append(phit[-1])
                    Err.append(curr_fit)
                    if self.Harmonicmethod:
                        self.harmfitstore.append(harmcurrentfit)

            else:
                phit.append(phit[-1])  # accept the current state
                Err.append(curr_fit)
                if self.Harmonicmethod:
                    self.harmfitstore.append(harmcurrentfit)

            t += 1

        return phit, Naccept, Err

    # MCMC metropolis hastings algoritm with adaptive covarience
    def adaptive_sampling(self, phit, Err, covar, prior_range, Nsamp, NMax, Naccept):

        # set up for the loop
        N = Nsamp * len(phit[0])
        t = N + 1  # iteration start point

        # define important parameters
        mean = np.array(phit[0])  # watch THIS
        a = 1

        if self.Harmonicmethod: # goes into harmonic percentage fit
            curr_fit, harmcurrentfit = self.simulate(phit[-1])  # experimental data gets passed in kwargs
        else:
            curr_fit = self.simulate(phit[-1])

        while t != NMax:

            s = t - N
            gamma = (s + 1) ** (-0.6)

            phi_trail = np.random.multivariate_normal(phit[-1], a * covar)  # gets the trail solution

            prob = prior_test(phi_trail, prior_range)  # Tests if solution is possible assuming uniform distrabution

            if prob != 0:  # checks to see if position is possible

                # runs simulation
                if self.Harmonicmethod:  # goes into harmonic percentage fit
                    logtrial, harmtestfit = self.simulate(phi_trail)  # experimental data gets passed in kwargs
                else:
                    logtrial = self.simulate(phi_trail)

                # Probability comparison code
                comp = np.exp(logtrial - curr_fit)
                r = min(1, comp)

                u = np.random.uniform(0, 1)  # generates random uniform distrabution
                if u < r:

                    phit.append(phi_trail)
                    curr_fit = logtrial  # resets the sum
                    Err.append(curr_fit)
                    if self.Harmonicmethod:
                        harmcurrentfit = harmtestfit
                        self.harmfitstore.append(harmcurrentfit)
                    accepted = 1
                    Naccept += 1  # counter for num accepted

                else:
                    phit.append(phit[-1])
                    Err.append(curr_fit)
                    if self.Harmonicmethod:
                        self.harmfitstore.append(harmcurrentfit)
                    accepted = 0

            else:
                phit.append(phit[-1])  # accept the current state
                Err.append(curr_fit)
                if self.Harmonicmethod:
                    self.harmfitstore.append(harmcurrentfit)
                accepted = 0

            # iteration probability modifiers
            covar = ((1 - gamma) * covar + gamma * np.outer((phit[-1] - mean), (phit[-1] - mean)))
            mean = (1 - gamma) * mean + gamma * phit[-1]
            a = np.exp(np.log(a) + gamma * (accepted - 0.25))  # 0.25 will need to be a varible

            t += 1

        rate = Naccept / t  # Calculates chain acceptance rate

        return phit, rate, Err

    # test for the values fit with the
    def prior_test(self, val_in, prior_range):
        N = len(val_in)
        prob = 1
        i = 0
        while i != N:
            if val_in[i] >= prior_range[i, 0] and val_in[i] <= prior_range[i, 1]:
                prob *= 1
            else:
                prob = 0
                break  # cancels while loop with prob = 0
            i += 1

        return prob



################################################################################################
######################################   ADMCMC_TOTCURR_METHOD  ################################
################################################################################################

def PINT_ADMCMC_TOTCURR_settings(settings):
    nvar = int(settings.iloc[1][0])
    n = settings[0].shape
    n = int(n[0])
    nscal = int(settings.iloc[nvar + 2][0])
    Nfunc = int(settings.iloc[nvar + nscal + 3][0])

    # Gets the scaling
    scalvarpre = settings[nvar + 2:nvar + 3 + nscal]  # gets the num of scaling to

    nband = int((n - 10 - nvar - 3 - Nfunc - nscal) / 2)

    scalvar = [[int(scalvarpre.iloc[0][0])]]
    i = 0
    x = scalvarpre[0].str.split(',', expand=True).astype('float')
    x = x.values

    # sets up scalar varibles
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

            # Seperates the CMA input data
    var = settings[2:nvar + 2]
    Nbarndersnatch = nvar + 4 + nscal + Nfunc  # counter cause im lazy
    bandwidth = settings[Nbarndersnatch:Nbarndersnatch + nband]
    harm_weights = settings[Nbarndersnatch + nband:Nbarndersnatch + 2 * nband]

    # Changes the seperated input values into
    var = var[0].str.split(',', expand=True).astype('float')
    bandwidth = bandwidth[0].str.split(',', expand=True).astype('float')
    bandwidth = bandwidth.values  # changes from a df to np.array
    harm_weights = harm_weights[0].str.split(',', expand=True).astype('float')
    harm_weights = harm_weights.values

    # starts count for settings
    Np = nvar + 2 * nband + nscal + 3 + Nfunc

    # truncation settings for time series analysis
    truntime = settings.iloc[Np + 1][0]
    truntime = truntime.split(',')

    x = [float(truntime[0]), truntime[1].strip()]  # strip is to remove tabs
    if x[1] != 'MAX':  # allows for name input
        x[1] = float(x[1])
    truntime = x


    # Optimization settings
    op_settings = [0, 0, 0, 0, 0, 0, 0, 0]  # prealocation
    op_settings[0] = int(settings.iloc[Np + 2][0])  # This one setts the optimization method
    datatype = int(settings.iloc[Np + 3][0])  # Gets the type of data were compairing against
    op_settings[1] = int(settings.iloc[Np + 4][0])  # this one does the num cores
    op_settings[2] = 2 ** float(settings.iloc[Np + 5][0])  # sets the number of datapoints to compair in the current

    # MCMC settings
    #bay_settings = [0, 0, 0, 0, 0]
    op_settings[3] = int(settings.iloc[Np + 6][0])  # MCMC initial algoritm trail's per varible
    op_settings[4] = int(settings.iloc[Np + 7][0])  # number of trail for overall chain
    op_settings[5] = float(settings.iloc[Np + 8][0])  # burnin period as ratio of chain legth
    # treatment of the multi noise

    op_settings[6] = float(settings.iloc[Np + 9])
    #op_settings[6] = op_settings[3].str.split(',', expand=True).astype('float')
    #op_settings[6] = op_settings[3].values[0, :]  # noise of fit

    op_settings[7] = int(settings.iloc[Np + 10][0])  # number of chain to run at once

    return var, bandwidth, harm_weights, op_settings, datatype, scalvar, funcvar, truntime

# standard class type for ADMCMC sampling of moduled MECSim
def STAND_ADMCMC_TOTCURR(var, op_settings, Method, MEC_set, Excurr):

    # sets up for the
    phi0, covar, prior_range = prior_man(var)
    MEC_set['phi0'] = phi0
    MEC_set['covar'] = covar
    MEC_set['prior_range'] = prior_range

    # Allocates wrapper  around iterative MECSim function
    model = standard_ADMCMC_Model(*[Excurr, Method], **MEC_set)
    run = model.ADMCMC_chain
    X0 = []

    Ncore = op_settings[1]

    for i in range(Ncore):
        v0 = []
        for j in range(len(phi0)):
            v0 = []
            v0.append(phi0[j])
        v0.append(np.random.randint(10**6))
        X0.append(v0)

    #X0.append([val_in"""initial values + some random seed"""

    # input output list to do parallel in wrapper
    with Pool(processes=Ncore) as p:
        output = p.map(run, X0)

    # Something to unpack output into (chain, rate, Err) for each respective chain
    Mchains = []
    Mrate = []
    MErr = []
    Mharmperfit = []
    for i in range(Ncore):
        Mchains.append(output[i][0])
        Mrate.append(output[i][1])
        MErr.append(output[i][2])
        if Method == 'HEPWsigma' or Method == 'Baye_HarmPerFit'or Method == 'Bayes_ExpHarmPerFit':
            Mharmperfit.append(output[i][3])

    return Mchains, Mrate, MErr, Mharmperfit

# PINTS ADMCMC TOTCURR Runner
def PINT_ADMCMC_TOTCURR(op_settings, DCAC_method, MEC_set,  Excurr):

    # will need to get values out of the MEC_set

    # sum function so that multiple optimizers can be called
    # Fix
    
    """mAY NEED TO PRELIM DECIMATE eXCURR"""
    dummy = np.linspace(0, len(Excurr), len(Excurr))  # empty time varible FIX
    # Create an object with links to the model and time series
    var, MEC_set, noise, noisebool = tpseries.Noise_handler(MEC_set)

    # fix MEC SETTINGS HERE AND EXTRACT NOISE PARAMETER
    model = ADMCMC_Model()
    model.__int__(**MEC_set)

    x0 = np.array(var.iloc[:][0])
    if noisebool:
        x0 = np.append(x0, noise[0])
    else:
        pass
    x0 = [x0]   # to allow use of one chain
    # model = partial(model1.simulate(**kwargs), **MEC_set)     # attaches the kwargs to the function
    problem = pints.SingleOutputProblem(model, dummy, Excurr)

    # Set up the uniform priors
    """FIX SO THAT LOG10 SCALE WORKS"""
    priors = prior_set(var, noise, noisebool)
    # sets up related covarience matrix
    covar = tpseries.sigma_set(priors)

    if noisebool:
        # Create a log-likelihood function (adds an extra parameter!)
        log_likelihood = pints.GaussianLogLikelihood(problem)
    else:  # noise is a parameter
        log_likelihood = pints.GaussianKnownSigmaLogLikelihood(problem, op_settings[6])

    # Create a posterior log-likelihood (log(likelihood * prior))
    log_prior = pints.UniformLogPrior(priors[0], priors[1])

    # Create a posterior log-likelihood (log(likelihood * prior))
    log_posterior = pints.LogPosterior(log_likelihood, log_prior)

    # variation for Multichain method
    if DCAC_method[2] == 1: # single chain
        pass
    else:
        i = 0
        x = []
        while i != DCAC_method[2]:
            x.append(x0[0]) # add something here to add x0 random
            i += 1
        x0 = x

    # Create mcmc routine
    # mcmc = pints.AdaptiveCovarianceMCMC(x0, sigma0=None)
    mcmc = pints.MCMCController(log_posterior, DCAC_method[2], x0,sigma0=covar, method=pints.AdaptiveCovarianceMCMC)

    # Add stopping criterion
    mcmc.set_max_iterations(int(op_settings[4]/op_settings[7]))
    mcmc.set_initial_phase_iterations(int(op_settings[3]*len(var.iloc[:][0])/op_settings[7]))

    for sampler in mcmc.samplers():
        sampler.in_initial_phase()
        sampler.set_target_acceptance_rate(rate=0.234)
        sampler.set_initial_phase(op_settings[3]*log_likelihood.n_parameters())


    # optionalization paralization (non optimal but it'll do
    Para_chain = False
    if DCAC_method[2] != 1:
        Para_chain = True
    mcmc.set_parallel(parallel=Para_chain)

    #mcmc.set_log_to_screen(False)

    #mcmc.parallel(op_settings[1])
    chains = mcmc.run()
    #accept_rate = mcmc.acceptance_rate()

    return chains


# output writter for this method
def PINT_ADMCMC_TOTCURR_output(filename,MCMC_time,Niter,accept_rate,MCMC_mean,colapsed_std,phi0,covar_std,cal_covar,corr_matrix,R_hat,header,MEC_set):

    N = len(phi0)

    f = open(filename, "w+")  # creates a blank .inp to be written

    # writes the input
    f.write('Input file name: %s\n' % (filename))
    f.write('Logic and fitting method used: %s, %s\n' % (header[1], header[2]))
    # date of the colour
    f.write('Date: %s\n' % (datetime.datetime.today().strftime('%d-%m-%Y')))
    # inserts time
    f.write('Completetion time: %f\n\n' % (MCMC_time))

    f.write('number of function Iterations: %i\n' % (Niter))

    f.write('Acceptance rate of MCMC: %f\n\n' % (accept_rate))

    # Statistics
    f.write('Statistical values\n')

    # mean_var_out = str(mean_var_out)
    f.write('MCMC 1D Values Mean parametervalues (mean,std):\n')
    i = 0
    while i != N:
        m = format_e(MCMC_mean[i])
        s = format_e(colapsed_std[i])

        f.write('Var %d: %s,%s\n' % (int(i + 1), m, s))
        i += 1

    f.write('\n')

    f.write('Comparison between input mean and MCMC mean:\n')
    i = 0
    while i != N:
        m = MCMC_mean[i] / phi0[i]
        f.write('Var %d: %.5f\n' % (int(i + 1), m))
        i += 1

    f.write('\n')

    # Perr different to AC
    # x = format_e(Perr)
    # f.write('Percentage Error: %s\n\n' % (x))

    f.write('Multivariate normal distribution statistics\n\nCovarience standard deviation\n')
    i = 0
    while i != N:
        f.write('Var %d: %s\n' % (int(i + 1), format_e(covar_std[i])))
        i += 1

    f.write('\n')  # seperator

    f.write('Covarience Matrix\n')

    i = 0
    while i != N:
        j = 0
        while j != N:
            f.write('%s\t' % (format_e(cal_covar[i, j])))
            j += 1
        f.write('\n')
        i += 1
    f.write('\n')  # seperator

    f.write('Correlation Matrix\n')
    i = 0
    while i != N:
        j = 0
        while j != N:
            f.write('%s\t' % (format_e(corr_matrix[i, j])))
            j += 1
        f.write('\n')
        i += 1

    if type(R_hat) != []:
        f.write('\n')
        f.write('Gelman and Rubin (R-hat)\n')
        for i in range(len(R_hat)):
            f.write('%.5f\t' %R_hat[i])

    f.write("\n " + header[1] + " specific outputs\n")
    if header[1] == 'Baye_HarmPerFit':
        sigma = MEC_set.get('Exsigma')
        f.write('harmonic Sigmas\n')
        f.write(np.array2string(sigma))
    f.close()

    return


########################################################################################################################
#######################################    CUSTOM ADMCMC ALGORITHIM    #################################################
########################################################################################################################

# sets up the ranges for the uniform priors
def prior_man(var):
    N = var.shape[0]  # gets the # of varibles

    # preallocation
    phi0 = np.zeros(N)
    covar = np.zeros((N, N))
    prior_range = np.zeros((N, 2))

    i = 0
    while i != N:
        phi0[i] = var.values[i][0]
        prior_range[i, 0] = var.iloc[i][1]
        prior_range[i, 1] = var.iloc[i][2]
        covar[i, i] = (prior_range[i, 1] - prior_range[i, 0]) ** 2 / 12  # varience of a uniform distrabution

        i += 1

    return phi0, covar, prior_range

# Extracts Covarience statistics from MCMC output
def covarinence_stat(phit, burnin):
    N = len(phit[0, :])  # extracts the cma best fit values

    Nmax = phit.shape[0]  # Gets the overall number of iterations

    MCMC_mean = []
    i = 0
    while i != N:
        MCMC_mean.append(np.mean(phit[int(Nmax * burnin):, i]))  # calculates and adds the mean value
        i += 1

    # single value standard deviations
    colapsed_std = []
    i = 0
    while i != N:
        colapsed_std.append(np.std(phit[int(Nmax * burnin):, i],ddof=1))  # calculates and adds the mean value
        i += 1

    # multivarible dist stats
    cal_covar = np.cov(np.transpose(phit[int(Nmax * burnin):, :]),ddof=1)  # gets the covarence matrix of the

    covar_std = np.sqrt(np.diag(cal_covar))  # extracts the sigmas from covar matrix

    sig_sq = np.outer(covar_std, covar_std)  # gets a sigma squared matrix

    corr_matrix = cal_covar / sig_sq  # Gets the coleation matrix

    return MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix
