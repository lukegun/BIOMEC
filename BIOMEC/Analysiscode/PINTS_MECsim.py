"""
PINTS MECSim interface





"""
# general house keeping modules
import time
import datetime
import sys
import os
import shutil
# Custom function imports
import ML_signal_processing as MLsp   #General ML_processing funciont
import CMA_modules as CMAmod
import MCMC_modules as MCMCmod
import timeseries_modules as tpseries
import Script_generator as Scripgen
import window_func as Wint
import plotting_scripts as plotter
import numpy as np
import pints
import matplotlib.pyplot as plt
# module type functions CMAES and ADMCMC
import CMAES as standCMA
import ADMCMC as standADMCMC

t1 = time.time()

# Set up input and output
input_name = sys.argv[1]
outputfname = input_name.split('.')[0] +'_output_'+str(datetime.date.today())

# genertes the output file checking if a previous one exists
outputfname = Scripgen.outputfilegenertor(outputfname)


shutil.copy(input_name,outputfname)

# Collects input data
settings, data, Exp_data = Scripgen.globalinreader(input_name)  # Takes the input data and seperates to its sections

#number of experimental inputs
Expnumber = len(Exp_data)

# changes CMA_settings to usable data
header, var, bandwidth, harm_weights, op_settings, datatype, scalvar, funcvar, truntime, Harmonicwindowing = tpseries.PINTS_Varibles(settings)

spaces = Scripgen.data_finder(data)

AC_freq, AC_amp = MLsp.AC_info(data, spaces)

#bandwidth adjustments to standardize the bandwidth notation
if Harmonicwindowing[0] == "squarewindow":
    # this is to meet standard notation of badwidth = whole range around frequency harmonic therefore square window uses half
    # This error has been corrected for other window functions
    for i in range(1,bandwidth.shape[0]):
        bandwidth[i][:] = bandwidth[i][:]/2


# note that for NORMHARM Excurr is the set of harmonics Decimates the experimental result
Excurr, Exp_t,Extime,sigma = tpseries.Exp_data_spliter(Exp_data, datatype, AC_freq, bandwidth,spaces,header,op_settings,Harmonicwindowing)  # splits the experimental input depending on input tyharm_weights = Scripgen.harm_weights_trans(harm_weights)  # reshapes harm_inputs into row

harm_weights = Scripgen.harm_weights_trans(harm_weights)

DCAC_method = MLsp.ACDC_method(AC_amp, op_settings)
print("Initial Set up succesful, starting " + header[2] + " using " + header[1] + " objective function.")
# General truncation of experiment files
if DCAC_method[0] == 1 and (header[1] == 'HarmPerFit' or header[1] == 'HEPWsigma' or header[1] == 'Baye_HarmPerFit' or header[1] == 'Bayes_ExpHarmPerFit'):  # AC method {DC = 0, AC = 1} FIX
    # Excurr truncates and sets up time truncation and time

    Excurr, sigma, Nsimdeci, Nex = MLsp.EXPharmtreatment(Excurr,sigma,Extime, truntime, op_settings)

    if Harmonicwindowing[0] == "Conv_guassian": # This is to generate the simulation windowing for ease
        filters = Wint.RGguassconv_filters(2**int(data.iloc[6]),bandwidth,Extime/2**int(data.iloc[6]),AC_freq,Harmonicwindowing[1])
    else:
        filters = None

    # Excurr is carried over harmstore set up previous function
elif header[1] == 'TCDS':

    # truncates and passes around experimental harmonics and find truncation points
    TotalExcurr = Excurr # here so that we can plot the output harmonics
    Excurr, Nsimdeci, Nex = MLsp.EXPtotcurrtreatment(Excurr, Extime, truntime, op_settings)
elif header[1] == 'FTC' or header[1] == 'Log10FTC':
    TotalExcurr = Excurr
    # Translates the current into the frequency domain
    Excurr = np.fft.rfft(Excurr)
    Excurr = Excurr/len(Excurr)  # adjustment for amplitude spectrum as fft does not come out that way
    frequens = np.fft.rfftfreq(len(Excurr), d=Extime/len(Excurr))
    Excurr, Nsimdeci, Nex = MLsp.EXPFTtreatment(Excurr, frequens, truntime)

    # settings to remove background issue in log10fit
    if header[1] == 'Log10FTC':
        Excurr, nsimdecionoff = MLsp.log10harmtunc(Excurr, frequens, bandwidth,truntime,AC_freq,Nsimdeci)
        Nsimdeci = nsimdecionoff  # carries the inputs into the Excurr

# bunch of set up stuff
space_holder = var[[4, 5]].values  # retains values for use in printed values in np.array
var, scalvar = Scripgen.line_assigner(var, spaces, scalvar)
cap_series = tpseries.capacatance0(data,spaces,var)  # extracts the input capacitance
funcvar_holder = MLsp.functionholder(var)  # issolates the functional parameters
phi0 = var.iloc[:][0].values # extracts starting points

if datatype == 2: #adjusment for time series comparison to
    dt_exp = Exp_t[1] # time interval
    op_settings[2] = len(Exp_t)


if header[1] == 'HarmPerFit':  # AC method

    if header[2] == 'CMAES':  # CMA-ES for current total comparison

        freq_set = {'harm_weights': harm_weights, 'bandwidth': bandwidth}  # here to allow functional ACDC functionality

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data': Excurr,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'harm_weights': harm_weights,'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq':AC_freq,
                   'filters':filters}

        # Results are the best fit and logtot are all trails tried
        results, logtot = standCMA.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 1  # sets up to use the Perr method
        MEC_set['DCAC_method'] = DCAC_method
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out  # sets up for pinted inp file

        # sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                            results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    elif header[2] == 'ADMCMC':  # CMA-ES for current total comparison
        print('This method and logic is currently not developed yet')

elif header[1] == 'Baye_HarmPerFit':  # AC method
    print("Method depreciated due to incorrect noise model")
    exit()
    """# correction to the truncated harmonics of sigma
    sigma = tpseries.Harm_sigma_normal(Excurr,sigma,Expnumber)

    if header[2] == 'CMAES':  # CMA-ES for current total comparison

        freq_set = {'harm_weights': harm_weights, 'bandwidth': bandwidth}  # here to allow functional ACDC functionality

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data':Excurr,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'harm_weights': harm_weights,'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq':AC_freq}

        # Results are the best fit and logtot are all trails tried
        results, logtot = standCMA.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 1  # sets up to use the Perr method
        MEC_set['DCAC_method'] = DCAC_method
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out  # sets up for printed inp file

        # sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        #nned to fix this

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                           results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    elif header[2] == 'ADMCMC':  # CMA-ES for current total comparison

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq, 'harm_weights': harm_weights}

        # Results are the best fit and logtot are all trails tried
        Mchains, Mrate, MErr, Mharmperfit = standADMCMC.STAND_ADMCMC_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # creates the parameters output files
        tpseries.MCMC_para_logger2(outputfname+'/'+'MCMC_parallel_log.txt', Mchains)  # prints the multi chain output

        # combines and plots the multplie chains to one
        dist = []
        accept_rate = sum(Mrate) / len(Mrate)
        Err = []

        harmperfit = [] # this is aonly for this method
        for i in range(int(op_settings[4] / op_settings[7])):
            for j in range(op_settings[7]):
                dist.append(Mchains[j][i])
                Err.append(MErr[j][i])
                harmperfit.append(Mharmperfit[j][i])

        # this is only for the harmonic fit
        MCMCmod.MCMC_harmfit_logger(outputfname + '/' + 'MCMC_harmper_log.txt', harmperfit)

        dist = np.array(dist)
        Err = np.array(Err)

        MCMCmod.MCMC_dist_logger(outputfname+'/'+'Iter_log.txt', dist, Err)  # saves the distrabution to a text file

        # something to calculate the covarence and statistic values
        MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = MCMCmod.covarinence_stat(dist, op_settings[5])

        fit = MCMC_mean  # sets up for printed inp file

        # prints and calculates the series of R-gelman statistic
        if op_settings[7] != 1:
            R_hat = MCMCmod.Gel_rub_test(outputfname+'/'+'R_hat',Mchains, Err, op_settings[5])
        else:
            R_hat = []

        # something to calculate the percentage best fit of the mean calculated harmonic
        Perr = tpseries.MCMC_harmfitpercentage(MCMC_mean, Excurr, True, **MEC_set)

        # something to print output text
        standADMCMC.PINT_ADMCMC_TOTCURR_output(outputfname+'/'+"result_file.txt", completion_time, len(dist), accept_rate, MCMC_mean,
                                               colapsed_std, phi0, covar_std, cal_covar, corr_matrix,R_hat, header, MEC_set)"""

elif header[1] == 'Bayes_ExpHarmPerFit':  # AC method

    EXPperrErr, sigma, EXPperrmean = tpseries.Harm_sigma_percentge(Exp_data, Excurr, truntime, AC_freq, bandwidth, spaces,  op_settings,Harmonicwindowing)

    # correction to the truncated harmonics of sigma
    #sigma = tpseries.Harm_sigma_normal(Excurr,sigma,Expnumber)

    if header[2] == 'CMAES':  # CMA-ES for current total comparison

        freq_set = {'harm_weights': harm_weights, 'bandwidth': bandwidth}  # here to allow functional ACDC functionality

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data':Excurr,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'harm_weights': harm_weights,'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq':AC_freq,
                   "EXPperrErr":EXPperrErr, 'filters':filters}

        # Results are the best fit and logtot are all trails tried
        results, logtot = standCMA.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 1  # sets up to use the Perr method
        MEC_set['DCAC_method'] = DCAC_method
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out  # sets up for printed inp file

        # sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        #nned to fix this

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                           results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    elif header[2] == 'ADMCMC':  # CMA-ES for current total comparison

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq, 'harm_weights': harm_weights,
                   "EXPperrErr":EXPperrErr, 'filters':filters}

        # Results are the best fit and logtot are all trails tried
        Mchains, Mrate, MErr, Mharmperfit = standADMCMC.STAND_ADMCMC_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # creates the parameters output files
        tpseries.MCMC_para_logger2(outputfname+'/'+'MCMC_parallel_log.txt', Mchains)  # prints the multi chain output

        # combines and plots the multplie chains to one
        dist = []
        accept_rate = sum(Mrate) / len(Mrate)
        Err = []

        harmperfit = [] # this is aonly for this method
        for i in range(int(op_settings[4] / op_settings[7])):
            for j in range(op_settings[7]):
                dist.append(Mchains[j][i])
                Err.append(MErr[j][i])
                harmperfit.append(Mharmperfit[j][i])

        # this is only for the harmonic fit
        MCMCmod.MCMC_harmfit_logger(outputfname + '/' + 'MCMC_harmper_log.txt', harmperfit)

        dist = np.array(dist)
        Err = np.array(Err)

        MCMCmod.MCMC_dist_logger(outputfname+'/'+'Iter_log.txt', dist, Err)  # saves the distrabution to a text file

        # something to calculate the covarence and statistic values
        MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = MCMCmod.covarinence_stat(dist, op_settings[5])

        fit = MCMC_mean  # sets up for printed inp file

        # prints and calculates the series of R-gelman statistic
        if op_settings[7] != 1:
            R_hat = MCMCmod.Gel_rub_test(outputfname+'/'+'R_hat',Mchains, Err, op_settings[5])
        else:
            R_hat = []

        # something to calculate the percentage best fit of the mean calculated harmonic
        Perr = tpseries.MCMC_harmfitpercentage(MCMC_mean, Excurr, True, **MEC_set)

        # something to print output text
        standADMCMC.PINT_ADMCMC_TOTCURR_output(outputfname+'/'+"result_file.txt", completion_time, len(dist), accept_rate, MCMC_mean,
                                               colapsed_std, phi0, covar_std, cal_covar, corr_matrix,R_hat, header, MEC_set)

elif header[1] == 'TCDS':  # PINTS (yt-Isim)**2 for total current

    Excurr = MLsp.deci_exp_sim(Excurr, DCAC_method[3])  # decimates the experimental to length of simulation to deci number

    if header[2] == 'CMAES' or header[1] == 'any other optimiser':   # CMA-ES for current total comparison

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq,'harm_weights':harm_weights}

        # Results are the best fit and logtot are all trails tried
        results, logtot = tpseries.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1)/60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 0  # sets up to use tdiffsquare
        MEC_set['DCAC_method'] =  DCAC_method
        MEC_set.update({'Exp_data':Excurr})
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out  # sets up for pinted inp file

        #sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                           results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    elif header[2] == 'ADMCMC'or header[1] == 'any other MCMC':   # MCMC for current total

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data': Excurr,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq,'harm_weights':harm_weights}

        # chains are just the log tot Chains is an array of arrays

        chains = standADMCMC.PINT_ADMCMC_TOTCURR(op_settings, DCAC_method, MEC_set, Excurr)
        completion_time = (time.time() - t1)/60
        # treat output here
        accept_rate = 0.234

        """fix with relative to accep_list"""
        # checks if multiple chains was run
        if op_settings[1] == 1:
            dist = chains[0]
        else:
            dist = tpseries.Para_combiner_pints(chains)
            Mchains = chains
        Niter = len(dist[:,0])

        MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = tpseries.covarinence_stat(dist, op_settings[5]) # op_settings[5] = burnin

        fit = MCMC_mean  # sets up for printed inp file

        # prints and calculates the series of R-gelman statistic
        if op_settings[7] != 1:
            R_hat = pints.rhat_all_params(chains)
            #R_hat = MCMCmod.Gel_rub_test(outputfname + '/' + 'R_Hat.txt', chains, Err, op_settings[5])
        else:
            R_hat = []

        standADMCMC.PINT_ADMCMC_TOTCURR_output(outputfname+'/'+"result_file.txt", completion_time, len(dist), accept_rate, MCMC_mean,
                                               colapsed_std, phi0, covar_std, cal_covar, corr_matrix,R_hat, header, MEC_set)

        # prints output as txt file
        tpseries.MCMC_para_logger2(outputfname+'/'+'MCMC_parallel_log.txt', chains)  # prints the multi chain output
        tpseries.MCMC_dist_logger2(outputfname+'/'+'Iter_log.txt', dist)  # saves the distrabution to a text file

elif header[1] == 'FTC':  # PINTS (yt-Isim)**2 for fourier transform of the current

    if header[2] == 'CMAES' or header[1] == 'any other optimiser':   # CMA-ES for current total comparison

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq,'harm_weights':harm_weights
                   ,'Extime':Extime,'bandwidth':bandwidth}

        # Results are the best fit and logtot are all trails tried
        results,logtot = tpseries.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1)/60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 0  # sets up to use tdiffsquare
        MEC_set['DCAC_method'] =  DCAC_method
        MEC_set.update({'Exp_data':Excurr})
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out  # sets up for printed inp file

        #sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                           results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    elif header[2] == 'ADMCMC'or header[1] == 'any other MCMC':   # MCMC for current total

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data': EX_hil_store,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq,'harm_weights':harm_weights
                    ,'Extime':Extime}

        # chains are just the log tot Chains is an array of arrays

        chains = standADMCMC.PINT_ADMCMC_TOTCURR(op_settings, DCAC_method, MEC_set, Excurr)
        completion_time = time.time() - t1
        # treat output here
        accept_rate = 0.2

        """fix with relative to accep_list"""
        # checks if multiple chains was run
        if op_settings[1] == 1:
            dist = chains[0]
        else:
            dist = tpseries.Para_combiner_pints(chains)
        Niter = len(dist[:,0])

        MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = tpseries.covarinence_stat(dist, op_settings[5]) # op_settings[5] = burnin

        fit = MCMC_mean  # sets up for printed inp file

        standADMCMC.PINT_ADMCMC_TOTCURR_output(outputfname+'/'+"global_out.txt", completion_time, len(dist), accept_rate, MCMC_mean,
                                               colapsed_std, phi0, covar_std, cal_covar, corr_matrix,R_hat, header, MEC_set)

        # prints output as txt file
        tpseries.MCMC_para_logger2(outputfname+'/'+'MCMC_parallel_logger.txt', chains)  # prints the multi chain output
        tpseries.MCMC_dist_logger2(outputfname+'/'+'global_logger.txt', dist)  # saves the distrabution to a text file

elif header[1] == 'Log10FTC':  # PINTS (yt-Isim)**2 for fourier transform of the current

    Excurr = np.abs(Excurr)     # gets teh power amplitude
    Excurr = np.log10(Excurr)   # sets up the experimental

    if header[2] == 'CMAES' or header[1] == 'any other optimiser':   # CMA-ES for current total comparison

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,'Extime':Extime,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq,'harm_weights':harm_weights}

        # Results are the best fit and logtot are all trails tried
        results,logtot = tpseries.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1)/60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 0  # sets up to use tdiffsquare
        MEC_set['DCAC_method'] =  DCAC_method
        MEC_set.update({'Exp_data':Excurr})
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out  # sets up for pinted inp file

        #sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                           results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    elif header[2] == 'ADMCMC'or header[1] == 'any other MCMC':   # MCMC for current total
        print('ADMCMC has not currently been developed for Log10FT logic')
        exit()
        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data': EX_hil_store,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,'Extime':Extime,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq,'harm_weights':harm_weights}

        # chains are just the log tot Chains is an array of arrays

        chains = standADMCMC.PINT_ADMCMC_TOTCURR(op_settings, DCAC_method, MEC_set, Excurr)
        MCMC_time = (time.time() - t1)/60
        # treat output here
        accept_rate = 0.234

        """fix with relative to accep_list"""
        # checks if multiple chains was run
        if op_settings[1] == 1:
            dist = chains[0]
        else:
            dist = tpseries.Para_combiner_pints(chains)
        Niter = len(dist[:,0])

        MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = tpseries.covarinence_stat(dist, op_settings[5]) # op_settings[5] = burnin

        fit = MCMC_mean  # sets up for printed inp file

        standADMCMC.PINT_ADMCMC_TOTCURR_output(outputfname+'/'+"global_out.txt", MCMC_time, len(dist), accept_rate, MCMC_mean,
                                               colapsed_std, phi0, covar_std, cal_covar, corr_matrix,R_hat, header, MEC_set)

        # prints output as txt file
        tpseries.MCMC_para_logger2(outputfname+'/'+'MCMC_parallel_logger.txt', chains)  # prints the multi chain output
        tpseries.MCMC_dist_logger2(outputfname+'/'+'global_logger.txt', dist)  # saves the distrabution to a text file

elif header[1] == 'HEPWsigma':
    # CMA_ES method for multivariate method
    if header[2] == 'CMAES':

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method, 'Exp_data': Excurr,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth,'Exsigma':sigma,'AC_freq':AC_freq,'harm_weights':harm_weights, 'filters':filters}

        # Results are the best fit and logtot are all trails tried
        results, logtot = standCMA.STAND_CMAES_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # plot logger file
        Scripgen.iter_logger(outputfname+'/Iter_log.txt',logtot)

        DCAC_method[0] = 1  # sets up to use the Perr method
        MEC_set['DCAC_method'] = DCAC_method
        var_out, Perr = standCMA.CMA_output(results, **MEC_set)

        fit = var_out # sets up for pinted inp file

        # sets up mean CMA_es values and STD
        mean_var_out = CMAmod.CMA_meanval(results, var)

        # treat output here
        standCMA.PINT_CMAES_TOTCURR_output(outputfname+'/result_file.txt',completion_time, space_holder, var_out,
                                           results, mean_var_out, DCAC_method,Perr,header,MEC_set)

    # multivariate ADMCMC
    elif header[2] == 'ADMCMC':

        MEC_set = {'data': data, 'var': var, 'spaces': spaces, 'DCAC_method': DCAC_method,
                   'scalvar': scalvar, 'funcvar': funcvar, 'funcvar_holder': funcvar_holder,
                   'cap_series': cap_series, 'op_settings': op_settings, 'Nsimdeci': Nsimdeci,
                   'bandwidth': bandwidth, 'Exsigma': sigma, 'AC_freq': AC_freq, 'harm_weights': harm_weights, 'filters':filters}

        # Results are the best fit and logtot are all trails tried
        Mchains, Mrate, MErr,Mharmperfit = standADMCMC.STAND_ADMCMC_TOTCURR(var, op_settings, header[1], MEC_set, Excurr)

        completion_time = (time.time() - t1) / 60

        # creates the parameters output files
        tpseries.MCMC_para_logger2(outputfname+'/'+'MCMC_parallel_log.txt', Mchains)  # prints the multi chain output

        #combines and plots the multplie chains to one
        dist = []
        accept_rate = sum(Mrate)/len(Mrate)
        Err = []
        for i in range(int(op_settings[4]/op_settings[7])):
            for j in range(op_settings[7]):
                dist.append(Mchains[j][i])
                Err.append(MErr[j][i])

        dist = np.array(dist)
        Err = np.array(Err)

        MCMCmod.MCMC_dist_logger(outputfname+'/'+'Iter_log.txt', dist,Err)  # saves the distrabution to a text file

        # something to calculate the covarence and statistic values
        MCMC_mean, colapsed_std, cal_covar, covar_std, corr_matrix = MCMCmod.covarinence_stat(dist, op_settings[5])

        fit = MCMC_mean  # sets up for printed inp file

        # prints and calculates the series of R-gelman statistic
        if op_settings[7] != 1:
            R_hat = MCMCmod.Gel_rub_test(outputfname+'/'+'R_Hat.txt',Mchains, Err, op_settings[5])
        else:
            R_hat = []

        # something to calculate the percentage best fit of the mean calculated harmonic
        Perr = tpseries.MCMC_harmfitpercentage(MCMC_mean,Excurr,True, **MEC_set)

        # something to print output text
        standADMCMC.PINT_ADMCMC_TOTCURR_output(outputfname+'/'+"result_file.txt", completion_time, len(dist), accept_rate, MCMC_mean,
                                               colapsed_std, phi0, covar_std, cal_covar, corr_matrix,R_hat, header, MEC_set)

    else:
        print('Incorrect optimisation method passed for HEPWsigma Fitting method')

else:
    print('Incorrect header files where used')
    exit()


 # rextracts the experimental data
x = header[1] # holds this so that the method can go through and only get the 1st current
# extracts total current
header[1] = 'TCDS'
Extotcurr, Exp_t, Extime, sigma = tpseries.Exp_data_spliter(Exp_data, datatype, AC_freq, bandwidth, spaces, header,
                                                             op_settings,Harmonicwindowing)  # splits the experimental input depending on input

# extravt harmonics
if AC_amp != 0:
    header[1] = 'HarmPerFit'
    Exharmcurr, Exp_t, Extime, sigma = tpseries.Exp_data_spliter(Exp_data, datatype, AC_freq, bandwidth, spaces,
                                                                    header, op_settings,Harmonicwindowing)  # splits the experimental input depending on input

header[1] = x

#something here to print the MECSim input best fit file
with MLsp.cd(outputfname):
    # writes the mecsim output to file
    fit = MLsp.MECsim_noisepararemover(fit,var)
    MLsp.MECSiminpwriter(fit, **MEC_set)

    # simulate the best fit
    Scurr = MLsp.Iterative_MECSim_Curr(*[fit], **MEC_set)

    # gets a linear time array
    Exptimearray = np.linspace(0,Extime,len(Extotcurr))

    Simtimearray = np.linspace(0,Extime,len(Scurr))
    Simvoltage = Scripgen.outputsetter(data,AC_freq, [AC_amp])

    # plots and saves output from mecsim
    MECheader = Scripgen.MECoutsetter(data, AC_freq, [AC_amp])
    if datatype != 2:
        Scripgen.outputwriter('MEC_Sim_output_bfit', MECheader, Simvoltage, Scurr, Simtimearray)

    # plot the total current
    #Scurr, Nsimdeci, Nex = MLsp.EXPtotcurrtreatment(Scurr, Extime, truntime, op_settings)
    plt.figure()
    plt.plot(Exptimearray, Extotcurr,label='Experimental')
    Ndeci = int(len(Scurr) / op_settings[2])
    plt.plot(Simtimearray, Scurr,label='Simulated')
    plt.legend()
    plt.xlabel('Time (sec)')
    plt.ylabel('Current (A)')
    plt.savefig('totcurrent.png')

    plt.close()

    #calculates  the current output
    if header[1] == 'FTC' or header[1] == 'Log10FTC':
        log10ex,freqex = plotter.logten_extractor(Extotcurr,Exptimearray[-1], truntime,bandwidth, AC_freq)
        log10sim, freqsim = plotter.logten_extractor(Scurr, Simtimearray[-1], truntime, bandwidth, AC_freq)

        plt.figure()
        plt.plot(freqex, log10ex)
        plt.plot(freqsim,log10sim)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Log10(I) (A)')
        plt.savefig('powerspectrumcomparison.png')
        plt.close()

    # plot the harmonics
    if AC_amp != 0:
        plotter.harmplot('harmonicplots', Scurr, Simtimearray, Exharmcurr, Exptimearray, bandwidth, AC_freq, spaces)

    #plot the probability distrabutions
    if header[2] == 'ADMCMC':
        # plot postieriour probability
        variblenames = plotter.density_plotter('Probabilityplots', dist, space_holder, op_settings[5])

        os.makedirs('Convergenceplots')

        # plot convergence of parameters
        if op_settings[7] != 1: # plot for multiple convergences chains
            plotter.multiconverg('Convergenceplots',Mchains,variblenames)
        else:           # plot for a singular convergence
            plotter.sinularconverg('Convergenceplots',dist,variblenames)

    elif header[2] == 'CMAES':
        # gets the variable names
        os.makedirs('Convergenceplots')
        variblenames = plotter.name_Allo(space_holder)	
        plotter.sinularconverg('Convergenceplots',logtot,variblenames)

    print('Calculation complete')
