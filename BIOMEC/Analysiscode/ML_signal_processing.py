# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:36:24 2018

@author: Luke Gundry
"""

from scipy.fftpack import rfft, irfft, rfftfreq
from scipy.signal import hilbert as anal_hil_tran # NEED signal to make a smooth envolope
import numpy as np
import os
from subprocess import run
from Script_generator import *
import time
from shutil import copyfile 
import mecsim  # imports modulated mecsim
# extract how mony harmonics  NEED TO INVOLVE THIs IN GLOBAl

#truncate harmonics at window regions

# returns AC frequency Line Number in MECSim input
def AC_info(data,spaces):
    
    ii = 1 # counter
    x = np.zeros((spaces[4], 2))   # preallocates a AC bin
    while ii != spaces[4] + 1:   # extracts all # of AC signals
        
        x[ii-1, :] = data.iloc[spaces[0] + 1 + ii][0].split(',')     #takes each AC input and splits frequency and amp       
        ii += 1
    
    AC_freq = np.sort(x[:,1]) # sorts frequency into lowest to highest
    AC_amp = sum(x[:,0])      # Sum of all the AC signals
    
    return AC_freq, AC_amp

# defines frequency windows to be cut and passed through inverse-fft
def ifft_windowing(AC_freq, freq, bandwidth,spaces):

    nharm = bandwidth.shape[1]
    g = freq[1] # gets gradent for translation to n space 
    f_window = np.zeros((spaces[4]*nharm, 2))
    
    j = 0   #  AC counter
    while j != spaces[4]:   # collects the AC input freqency
        f = AC_freq[j]
        i = 0   # harmonic counter
        while i != nharm:   # collects the harmonics for each freqency
            width = bandwidth[j + 1,i]
            f_window[(j*nharm + i),0] = (i + 1)*f - width
            f_window[(j*nharm + i),1] = (i + 1)*f + width
            i += 1
        j += 1
    
    # move frequency to n space
    n_window = (2/g)*f_window         #ROUND COULD BE THERE
    
    return n_window, f_window, g

# generates harmonics envolopes from MECSim current output    
def Harmonic_gen(fft_res, n_window, spaces, bandwidth, g): # cuts the harmonic fourier space data out {NEEDS TO EXCLUDE FUNDEMENTAL}

    Mwin = 2*np.max(bandwidth) # finds max bandwidth   {COULD BE SEPRATE}
    nharm = bandwidth.shape[1]     # counts the number of harmonics
    Np = len(fft_res)    # gets the number of datapoints
    # harmstore = np.zeros(nharm*spaces[4] ,int(Mwin*(2/g)))
    hil_store = np.zeros((nharm*spaces[4] + 1 ,Np))
    N_ifft = np.zeros((nharm*spaces[4] + 1))
    
    
    #extracts the fundimental haronic
    N_fund = int(round(bandwidth[0,0]*2/g))
    x  = fft_res[0:N_fund]  # Cuts out windowed harmonics 
    y = np.zeros(Np)
    jj=0
    while jj != N_fund:  
        y[jj] = x[jj]
        jj += 1
    hil_store[0,:] = irfft(y)  # generates fundimental harmonics from cut window
       
    
    j = 0
    while j != spaces[4]:
    # something to take fft and cut windows into the four space thing
        
        i = 0
        while i!= nharm:
            wl = int(round(n_window[j*nharm + i,0]))
            wh = int(round(n_window[j*nharm + i,1]))
            x = fft_res[wl:wh]   #This need to be truncated in the future
            y = np.zeros(Np)
            jj=wl
            while jj != wh:
                y[jj] = x[jj-wl]
                jj += 1
            harmonic = irfft(y)   # generates harmonics
            N_ifft[(j*nharm + i) + 1] = int(len(x))      # gets N length of harmonic
            N = int(N_ifft[(j*nharm + i) + 1])       # turns number to interger
            hil_store[(j*nharm + i) + 1, :] = abs(anal_hil_tran(harmonic)) #uses HILBERT TRANSFORM to generate the envolope
            # using the abs fixed an issue with the complexes disapearing in the return
            i += 1
        j += 1
 
    return hil_store, N_ifft

# DC fitter compares experimental DC current to simulated
def DC_Fitt(Curr, ExpI):  #Compairs current, wieghting functions can be put in here
    
    Nex = len(ExpI)
    Nsim = len(Curr)
    
    #plt.plot(np.arange(Nex),ExpI,np.arange(Nsim),Curr)
    
    # this is done to allow different input size comparison
    if Nex > Nsim: 
        C = Nex/Nsim    # C needs to e larger then one ....
        N = Nsim
    else:
        C = Nsim/Nex    # far nicer better
        N = Nex
    
    
    #   preallocation of error
    Ier = np.zeros(Nsim)

    i = 0   # Counter for while loop that goes through every i value
    while i != N - 1:   # Non-linear fitting for the current difference section
        j = int(round(np.round(i*C)))  # sets point on larger array to take
        if Nex > Nsim:  # If more experimental data then sim
            Ier[i] = abs(ExpI[j] - Curr[i])**2
        else:    # If more experimental data then sim
            Ier[i] = (ExpI[i] - Curr[j])**2
            
        i += 1
    
    MAPE = sum(Ier)  #calculates and shows Mean difference error
    
    return(MAPE)
    

# AC fitter compares experimental Harmonics to simulated
def AC_Fitt(hil_store, N_ifft, Exp_harm,harm_weights): # compairs harmonic evelopes of the reaction 
    
    nc = hil_store.shape[1]  
    nr = hil_store.shape[0] # gets the number of harmonics in the system 

    # general prellocation
    Ier = np.zeros((nr,nc))  # fitting maths
    Peak_err = np.zeros(nr)  # somthing to test the peaks   
    error = np.zeros(nr) 
    t1 = time.time()
    # compair hil_store to experimental harmonics
    j = 0
    while j != nr:
        
        Exp_env = np.array(np.trim_zeros(Exp_harm[j,:]))
        sim_env = np.array(np.trim_zeros(hil_store[j,:]))
        
        Nex = len(Exp_env)
        Nsim = len(sim_env)
            
        if Nex > Nsim: 
            C = Nex/Nsim # C needs to e larger then one ....
            N = Nsim
            iex = 1 # boolen data to speed up if loop        
        else:
            C = Nsim/Nex    # far better
            N = Nex
            iex = 0
        
        i = 0
        while i != N:
            
            jj = int(np.rint(i*C))  # might need better
            if iex == 1:  # If more experimental data then sim
                Ier[j,i] = (Exp_env[jj] - sim_env[i])**2
            else:    # If more experimental data then sim
                Ier[j,i] = (Exp_env[i] - sim_env[jj])**2
            i += 1

        error[j] = (sum(Ier[j,:]) + Peak_err[j])    
        
        j += 1
        
    S = np.dot(harm_weights, error)
    MAPE = np.sum(S)
    
    return MAPE

#Automaticly identifies if AC or DC methods should be used then allocates it
def ACDC_method(AC_amp,op_settings):     
        
    S = np.sum(AC_amp)
    DCAC_method = [0,0,0,0]     # here to assign the fitting method and ACDC method
    
    #Optimization settings
    #op_settings = [0,0]     #prealocation
    DCAC_method[1] = op_settings[0]  # sets the fitting function with {0 = absfit, 1 = %fit}
    DCAC_method[2] = op_settings[1]  # Sets the number of cores
    DCAC_method[3] = op_settings[2]  # Sets the decimation number
    
    if  S != 0:
        DCAC_method[0] = 1
        
    else:
        
        DCAC_method[0] = 0
                
    return DCAC_method

      
# this is the draft of the loop to go inside the machine learning
"It might be better to use harm_weightss as a setting varible for ML"
def Iterative_MECSim(*args, **kwargs): # (data, var, val_in, spaces, DCAC_method, Exp_data,harm_weights):
    
    val_in = args[0] # extracts args to the value being modified as a list
    
    # Extracts the varibles 
    data = kwargs.get('data')
    var = kwargs.get('var')
    spaces = kwargs.get('spaces')
    DCAC_method = kwargs.get('DCAC_method')
    Exp_data = kwargs.get('Exp_data')
    scalvar = kwargs.get('scalvar')
    funcvar = kwargs.get('funcvar')
    funcvar_holder = kwargs.get('funcvar_holder')
    cap_series = kwargs.get('cap_series')
    
    # extracts and removes the counter
    """Counter needs to be removed"""
    #val_in = val_in[0:-1]
    Ncore = DCAC_method[2]
    
    t1 = time.time()
    data = Var2data(data,var,val_in,spaces)  # turns the varibles into data 
    
    # sets up for scalar dependant varibles
    if scalvar[0][0] != 0:
        data = Scalvar2data(scalvar, val_in, data, spaces)
    else:
        pass

    # functional varibles
    if funcvar[0][0] != 0:
        funclist = listin_func(val_in, funcvar, funcvar_holder)
        data = funcvar2data(funclist, spaces, data)
    else:
        pass

    # modifies the modelled double layer capacitance
    if cap_series[0]:

        data = cap2data(val_in, var, spaces, data,cap_series[1])

    else:
        pass
    
    """start module reader"""
    MECSettings,Numberin,Reacpara,reacMech,Capaci,Diffvalues,ACSignal,EModown,\
        Currenttot, Error_Code,Ecustom,Tcustom = ModMec_form(data,spaces)
        
    mecsim.mecsim_main(MECSettings, Numberin, Reacpara, reacMech, Capaci, Diffvalues, ACSignal, EModown,Ecustom,Tcustom, Currenttot,
                       Error_Code)    # runs mecsim
    
    Currenttot = Currenttot[0:int(Numberin[0]),]    # truncates the static varibles of Fortran77    
    """end module reader"""
        
        
    deci = int(len(Currenttot)/DCAC_method[3]) # sets the decimation number of the harmonics
    
    if DCAC_method[0] == 1: #experiment uses AC Method
        
        harm_weights = kwargs.get('harm_weights') # need that flexability function
        bandwidth = kwargs.get('bandwidth')
        
        #Exp_harm = Exp_data
        
        fft_res = rfft(Currenttot)
        # this can be done seperatly So lok into it
        # calculate time step (dt = (2*(Efin-Est)/v)/Ndata)
        dt = (2*(MECSettings[3]-MECSettings[2])/MECSettings[5])/Numberin[0]
        
        freq = rfftfreq(len(fft_res),d = dt)   # for ML could be done at the start to stop iterations ALSO DOUBLE VECTOR

        # start of AC signal processing {TEST FROM HERE}
        AC_freq, AC_amp = AC_info(data,spaces)      # extracts the AC signal
        
        n_window, f_window, g =  ifft_windowing(AC_freq, freq, bandwidth, spaces)     # Extracts the frequency and n windowing {CAN BE SEPERATE}
        
        hil_store, N_ifft = Harmonic_gen(fft_res, n_window, spaces, bandwidth, g)  # extract envolope data
        
        hil_store = hil_store[:,0::deci] # decimates the simulation data for effecency        
        
        if DCAC_method[1] == 0: # use abs squared fitter
            c = AC_Fitt(hil_store, N_ifft, Exp_data,harm_weights)
        elif DCAC_method[1] == 1:    # use percentage fitter
            c = percentage_errorAC(hil_store, N_ifft, Exp_data,harm_weights)
            
        else:
            pass
        
    else:   # experiment use DC methd
        
        Scurr = Currenttot[0::deci]
        
        if DCAC_method[1] == 0: # use abs squared fitter
            c = DC_Fitt(Scurr, Exp_data)  # May not work
        elif DCAC_method[1] == 1:    # use percentage fitter
            c = percentage_errorDC(Scurr, Exp_data)
        else:
            pass
    print(time.time() - t1)

    return c


def functionholder(var):

    funcvar_holder = []
    # issolates up the functional parameters
    i = 0
    while i != len(var.iloc[:][0]):
        if var.iloc[i][4] == 41 or var.iloc[i][4] == 42:
            x = var.iloc[i].values
            x = np.append(x, i + 1)  # adds the varible column data to the area we need
            funcvar_holder.append(x)

        else:
            pass
        i += 1

    return funcvar_holder


# this is the draft of the loop to go inside the machine learning
"It might be better to use harm_weightss as a setting varible for ML"


def Iterative_MECSim_Curr(*args, **kwargs):  # (data, var, val_in, spaces, DCAC_method, Exp_data,harm_weights):

    val_in = args[0]  # extracts args to the value being modified as a list
    #print(val_in)
    # Extracts the varibles
    data = kwargs.get('data')
    var = kwargs.get('var')
    spaces = kwargs.get('spaces')
    #Exp_data = kwargs.get('Exp_data')
    scalvar = kwargs.get('scalvar')
    #DCAC_method = kwargs.get('DCAC_method')
    funcvar = kwargs.get('funcvar')
    funcvar_holder = kwargs.get('funcvar_holder')
    cap_series = kwargs.get('cap_series')
    # extracts and removes the counter
    #counter = val_in[-1]
    #val_in = val_in[0:-1]
    #Ncore = DCAC_method[2]

    data = Var2data(data, var, val_in, spaces)  # turns the varibles into data

    # sets up for scalar dependant varibles
    if scalvar[0][0] != 0:
        data = Scalvar2data(scalvar, val_in, data, spaces)
    else:
        pass

    if funcvar[0][0] != 0:
        funclist = listin_func(val_in, funcvar, funcvar_holder)
        data = funcvar2data(funclist, spaces, data)
    else:
        pass

    # checks to see if capacitance is present


    if cap_series[0]:

        data = cap2data(val_in, var, spaces, data,cap_series[1])

    else:
        pass

    """start module reader"""
    MECSettings, Numberin, Reacpara, reacMech, Capaci, Diffvalues, ACSignal, EModown, \
    Currenttot, Error_Code,Ecustom,Tcustom = ModMec_form(data, spaces)
    t1 = time.time()

    try:
        mecsim.mecsim_main(MECSettings, Numberin, Reacpara, reacMech, Capaci, Diffvalues, ACSignal, EModown,Ecustom,Tcustom, Currenttot,
                       Error_Code)  # runs mecsim
    except:
        print("MECSim iteration obvioulsy failed in some way")

    Scurr = Currenttot[0:int(Numberin[0]), ]  # truncates the static varibles of Fortran77

    #print("completed")
    #plt.figure()
    #plt.plot(Scurr)
    #plt.savefig("Current"+str(np.random.randint(20000000))+".png")

    print(time.time() - t1)

    #del globals()[Diffvalues]

    return Scurr


def MECSiminpwriter(val_in, **kwargs):  # (data, var, val_in, spaces, DCAC_method, Exp_data,harm_weights):

    # Extracts the varibles
    data = kwargs.get('data')
    var = kwargs.get('var')
    spaces = kwargs.get('spaces')
    #Exp_data = kwargs.get('Exp_data')
    scalvar = kwargs.get('scalvar')
    DCAC_method = kwargs.get('DCAC_method')
    funcvar = kwargs.get('funcvar')
    funcvar_holder = kwargs.get('funcvar_holder')
    cap_series = kwargs.get('cap_series')
    # extracts and removes the counter
    #counter = val_in[-1]
    #val_in = val_in[0:-1]
    #Ncore = DCAC_method[2]

    data = Var2data(data, var, val_in, spaces)  # turns the varibles into data

    # sets up for scalar dependant varibles
    if scalvar[0][0] != 0:
        data = Scalvar2data(scalvar, val_in, data, spaces)
    else:
        pass

    if funcvar[0][0] != 0:
        funclist = listin_func(val_in, funcvar, funcvar_holder)
        data = funcvar2data(funclist, spaces, data)
    else:
        pass

    # checks to see if capacitance is present


    if cap_series[0]:


        data = cap2data(val_in, var, spaces, data,cap_series[1])

    MECSimwriter(data)

    return

#removes the noise parameter for the writing function
def MECsim_noisepararemover(fit,var):
    x = []
    for i in range(len(fit)):
        if var.iloc[i][4] == 13:
            pass
        else:
            x.append(fit[i])
    return x

#extracts data from Mecsim file  (data is a pandas dataframe)
def input_sep(data, spaces,nstoremax = 10000000):
    
    # predefine input arrays
    #Numberin = [2**datapoints, NEVramp + 2, ACsource, Nspecies, NCap+1, Nreactions]

    # Set up parameters
    PReacpara = []
    PreacMech = []
    PCapaci = []
    PDiffvalues = []    # diff [conc,diff,surf]
    Numberin = [0,0,0,0,0,0]
    PACSignal = []     # should be size (Nac, 2) [Amp,freq]
    PEModown = [] #[AdEst, AdEend, E_rev1, E_rev2] # should be size 4 + inf 
    MECSettings = [] # should be size 32
    
    # set up parameters
    NVramp = int(data.iloc[19][0]) # gets number of NVramp

    # Numberin
    Numberin[0] = 2**float(data.iloc[6][0])     # Ndatapoint
    Numberin[1] = int(data.iloc[19][0])   # Nev
    Numberin[2] = int(data.iloc[33 + NVramp][0]) # NAC
    Numberin[3] = int(data.iloc[34 + NVramp + Numberin[2]][0]) #Nspecies
    Numberin[4] = int(data.iloc[35 + NVramp + Numberin[2] + Numberin[3]][0]) + 1 # Ncap+1
    Numberin[5] = len(data[:][0]) - (37 + NVramp + Numberin[2] + Numberin[3] + Numberin[4]) #Nmech

    #Extract Data into MECSettings 
    """MECSettings = [Temp,Resistance, E_start, E_rev, Ncyc, scanrate, datapoints, 
               Digipotcompad, outputtype, ECtype, pre_equalswitch, fixtimestep,
               Nfixedsteps, beta, Dstar_min, maxvoltstep, timeres, debugout, 
               Advolt, NEVramp,Geotype, planararea, Nsphere, sphererad, Ncyilinders,
               rad_cylinder, cylLen, spacialres, RDErad, RDErots, RDEkinvis, Epzc]"""
    i = 0
    while i != len(data.iloc[:][0]):
        # thing to skip the boring shit
        
        if i == 19:  #constant value of NEVramp
            MECSettings.append(float(data.iloc[i][0]))
            PEModown.append(float(data.iloc[i + 1][0]))   # gets mod V ramp start
            PEModown.append(float(data.iloc[i + 2][0]))   # gets mod V ramp end
            j = 0
            while j != NVramp:
                PEModown.append(float(data.iloc[i + 2 + j + 1][0]))
                j += 1
                
            i = int(i + 2 + NVramp)
        
        elif  i == 19 + 2 + NVramp + 12: #kinematics point ready for AC
            
            j = 0
            while j != Numberin[2]:
                x = []  # holder for AMP,freq
                y = data.iloc[i + 1+j][0].split(',')
                x.append(float(y[0])) # get amp
                x.append(float(y[1])) # get freq
                PACSignal.append(x)
                j += 1
            i = i + 1 + Numberin[2]
            
            
        elif i == spaces[1]: # Diffusion values to be put into the crap
            j = 0
            while j != Numberin[3]:
                x = []  # holder for AMP,freq
                y = data.iloc[spaces[1]+j][0].split(',')
                x.append(float(y[0])) # get Conc
                x.append(float(y[1])) # get Diff
                x.append(float(y[2])) # get surf
                PDiffvalues.append(x)
                j += 1
            i = i - 1 + Numberin[3]

        elif i == spaces[1]+Numberin[3]: # Cap values to be put into the crap
            PCapaci.append(float(data.iloc[spaces[1]+Numberin[3]][0]))
            MECSettings.append(float(data.iloc[i+1][0]))
            j = 0
            while j != Numberin[4]:
                PCapaci.append(float(data.iloc[spaces[1]+2+j+Numberin[3]][0]))
                j += 1 
            i = spaces[2] - 1
            
        elif  i == spaces[2]: # Kinetic mechanism
            
            j = 0
            while j != Numberin[5]:
                x1 = []  # Mechanism
                x2 = []  # mech para
                y = data.iloc[spaces[2]+j][0].split(',')

                k = 0
                while k !=  Numberin[3]+1:
                    
                    x1.append(float(y[k]))
                    k += 1
                PreacMech.append(x1)
                while k != Numberin[3]+6:
                    x2.append(float(y[k]))
                    k += 1
                PReacpara.append(x2)
                j += 1
            i = len(data.iloc[:][0]) -1 

        else:

            MECSettings.append(float(data.iloc[i][0]))
            
        i += 1

    """example    
    reacMech1 = [[0,-1,1,0],[2,0,-1,1]]
    Reacpara1 = [[2,10 ,0,1000,0.5],[2,2 ,0.1,2000,0.5]]
    Capacipre = [NCap, cap0, cap1, cap2, cap3, cap4]
    Diffvaluespre = [[Conc1,Diff1,surf1],[Conc2,Diff2,surf2],[Conc3,Diff3,surf3]]
    ACSignalpre = [[A1, freq1],[A2, freq2]]
    EModownpre = [AdEst, AdEend, E_rev1, E_rev2]"""

    # set up some placeholder customm values
    # Empty BUT ADJUST TO FIT THE NUMBER OF CUSTOM values
    #MECSettings[18] = 0 # sets the output type to take the below values in  2 for custom

    # preallocate the functions for custom time and potential values
    # likely this will be needed to loaded in from a settings files in future
    Ecustom = np.zeros(nstoremax)
    Tcustom = np.zeros(nstoremax)
   
    return MECSettings,Numberin,PReacpara,PreacMech,PCapaci,PDiffvalues,PACSignal,PEModown,Ecustom,Tcustom
 
#something to load MASTER.file (data is a pandas dataframe)
def ModMec_form(data,spaces):

    # this is for preallocating some lengths
    nstoremax=10000000
    
    # extracts values from data 
    MECSettings,Numberin,PReacpara,PreacMech,PCapaci,PDiffvalues,PACSignal,PEModown,Ecustom,Tcustom \
    = input_sep(data, spaces,nstoremax)
    
    # turns values into 
    Reacpara,reacMech,Capaci,Diffvalues,ACSignal,EModown,Currenttot,\
     Error_Code = pre2post(PReacpara,PreacMech,\
                                      PCapaci,PDiffvalues,PACSignal,PEModown,nstoremax)


    return MECSettings,Numberin,Reacpara,reacMech,Capaci,Diffvalues,ACSignal, EModown,Currenttot, Error_Code, Ecustom,Tcustom
    


def pre2post(PReacpara,PreacMech,PCapaci,PDiffvalues,PACSignal,PEModown,nstoremax=10000000):
    
    # Outputs
    Error_Code = 0
    
    # Inputs
    # header file
    
    nrmax=30
    nsigmax=30
    nrmaxtot=100
    nmaxErev=10000
    nsp=30
    nmaxcapcoeff = 4

    # static files
    Currenttot =  np.zeros(nstoremax)   # the plus ones are due to fortran starting at zero
    EModown = np.zeros(nmaxErev+1)
    ACSignal = np.zeros((nsigmax+1,2))
    Diffvalues = np.zeros((nrmax+1,3))
    reacMech = np.zeros((nsp+1,nrmax+1))
    Reacpara = np.zeros((nsp+1,4+1))
    Capaci = np.zeros(nmaxcapcoeff+2)

    """"Formater for Numberin"""


    """Inputs values loaded into a static varibles for fortran uses lists"""
    #EModown
    i = 0
    while i != nmaxErev+1:
        if i< len(PEModown):
            EModown[i] = PEModown[i]
            i += 1
        else:
            i = nmaxErev+1
        
    #Capaci
    i = 0
    while i != nmaxcapcoeff+2:
        if i < len(PCapaci):
            Capaci[i] = PCapaci[i]
            i += 1
        else:
            i = nmaxcapcoeff+2


    #ACSignal
    i = 0
    while i != nsigmax+1:
        if i < len(PACSignal):
            ACSignal[i,:] = PACSignal[i]
            i += 1
        else:
            i = nsigmax+1
        
    #Diffvalues
    i = 0
    while i != nrmax+1:
        if i< len(PDiffvalues):
            Diffvalues[i] = PDiffvalues[i]        # seems to work better with a deep copy for some reason
            i += 1
        else:
            i = nrmax+1
        
        
    # Reaction Parameters
    i = 0
    while i != nsp+1:
        if i< len(PReacpara):
            Reacpara[i] = PReacpara[i]
            i += 1
        else:
            i = nsp+1

    #reactionMech
    i = 0
    while i != nsp+1:
        if i< len(PreacMech):
            j = 0
            while j < len(PreacMech[0]):
                reacMech[i,j] = PreacMech[i][j]
                j += 1
            i += 1
        else:
            i = nsp+1
    
    
    return Reacpara,reacMech,Capaci,Diffvalues,ACSignal,EModown,Currenttot, Error_Code


def percentage_errorDC(Curr, ExpI):
    
    # need a thing to filter small zeros
    Nex = len(ExpI)
    Nsim = len(Curr)
    
    #plt.plot(np.arange(Nex),ExpI,np.arange(Nsim),Curr)
    
    # this is done to allow different input size comparison
    if Nex > Nsim: 
        C = Nex/Nsim    # C needs to e larger then one ....
        N = Nsim
    else:
        C = Nsim/Nex    # far nicer better
        N = Nex
    
    Iersq = np.dot(ExpI,ExpI)
    
    """NEED TO filter VALUES ABOVE NOISE"""
    #   preallocation of error
    Ier = np.zeros(N)
    Exer = np.zeros(N)

    i = 0   # Counter for while loop that goes through every i value
    while i != N - 1:   # Non-linear fitting for the current difference section
        
        j = int(round(np.round(i*C)))  # sets point on larger array to take
        if Nex > Nsim:  # If more experimental data then sim
            Ier[i] = ((ExpI[j] - Curr[i])**2)
            Exer[i] = (ExpI[j]**2)
        else:    # If more experimental data then sim
            Ier[i] = ((ExpI[i] - Curr[j])**2)
            Exer[i] = (ExpI[i]**2)
        i += 1
    
    Perr = np.sqrt(abs(np.sum(Ier)/np.sum(Exer)))  # calculates and shows Mean difference error
    
    return Perr

"""Sum percentage error for all haronics present"""
def percentage_errorAC(hil_store, N_ifft, Exp_harm,harm_weights):
    
    # need a thing to filter small zeros
    nc = hil_store.shape[1]  
    nr = hil_store.shape[0] # gets the number of harmonics in the system 

    # general prellocation
    Ier = np.zeros((nr,nc))  # fitting maths
    Peak_err = np.zeros(nr)  # somthing to test the peaks    
    Perr = np.zeros(nr)
    count = np.zeros(nr)
    t1 = time.time()
    # compair hil_store to experimental harmonics
    j = 0
    
    while j != nr:
        
        Exp_env = np.array(np.trim_zeros(Exp_harm[j,:]))
        sim_env = np.array(np.trim_zeros(hil_store[j,:]))
        
        Nex = len(Exp_env)
        Nsim = len(sim_env)
        
        if Nex > Nsim:
            
            C = Nex/Nsim # C needs to e larger then one ....
            N = Nsim
            iex = 1 # boolen data to speed up if loop 
            Iersq = np.dot(Exp_env[0::round(C)],Exp_env[0::round(C)])
            
        else:
            
            C = Nsim/Nex    # far better
            N = Nex
            iex = 0
            Iersq = np.dot(Exp_env,Exp_env)
            
        i = 0
        while i != N:
            
            jj = int(round(np.rint(i*C)))  # might need better
            if iex == 1:  # If more experimental data then sim
                Ier[j,i] = (Exp_env[jj] - sim_env[i])**2
            else:    # If more experimental data then sim
                Ier[j,i] = (Exp_env[i] - sim_env[jj])**2
            i += 1
        
        Perr[j] = np.sqrt(np.sum(Ier[j,:])/(Iersq))    # assigns the Perr as a row
        j += 1

    return Perr

# translates Perr into a matrix from a row
def Perr_ACtran(Perr, spaces):
    
    N = len(Perr)
    Nac = spaces[4]
    Nharm = int((N-1)/Nac)
    x = np.zeros((Nac+1, Nharm))
    
    x[0,0] = Perr[0]
    
    j = 1
    jj = 1
    Nac = spaces[4]
    while j != Nac + 1:
        i = 0
        while i != Nharm:
            x[j,i] = Perr[jj]
            jj += 1
            i += 1
        j += 1
    
    return x    # spits out the modified Perr

# Generates the multiple files for mecsim to run in
def parallelgen(op_settings):
    
    Ncore = op_settings[1]  # Extracts the nuber of cores avalible
    dirname = 'workingdir'  # sets the name of the working directory
    cwd = os.getcwd() # gets the current working directory
    
    # Creates a folder for each core process then copys mecsim into the folder
    i = 0
    while i != Ncore:
        """Might need something here for allowing functioning"""
        s = '%s/%s%i' %(cwd,dirname, i+1)  # name for working directory
        os.makedirs(s)
        s = '%s/MECSim' %(s)
        copyfile('MECSim',s)
        
        #CHMOD Stuff for permission changing (st_mode=33279)
        st = os.stat('MECSim')
        old_mode = st.st_mode   # this is what were changing for use
        new_mode = old_mode | 33279 # 33279 allows for read and write as programs and all users
        os.chmod(s, new_mode)
        
        i += 1
    
    return

# This is the CD directory sections
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


# does the functional code to an input list
def listin_func(val_in, funcvar, funcvar_holder):
    funclist = []
    i = 0
    while i != funcvar[0][0]:
        Nfunc = int(funcvar_holder[i][6]) - 1  # gets varible row
        if funcvar[i + 1][2] == 0:  # if functional parameter is unpaired
            if funcvar_holder[i][4] == 41:  # R functional parameter
                # sets Kf
                Kf = val_in[Nfunc] * np.sin(funcvar[i + 1][3] * (np.pi / 180))
                # sets kb
                Kb = val_in[Nfunc] * np.cos(funcvar[i + 1][3] * (np.pi / 180))

                funclist.append(np.array([Kf, Kb, 1, funcvar_holder[i][5]]))

            elif funcvar_holder[i][4] == 42:  # theta functional parameter

                Kf = funcvar[i + 1][3] * np.sin(val_in[Nfunc] * (np.pi / 180))
                # sets kb
                Kb = funcvar[i + 1][3] * np.cos(val_in[Nfunc] * (np.pi / 180))

                funclist.append(np.array([Kf, Kb, 1, funcvar_holder[i][5]]))

            else:
                pass

        else:  # if functional parameter is paired
            Ncoup = int(funcvar[i + 1][3]) - 1  # gets paired parameter listin input
            if funcvar_holder[i][4] == 41:  # R functional parameter
                # sets Kf
                Kf = val_in[Nfunc] * np.sin(val_in[Ncoup] * (np.pi / 180))
                # sets kb
                Kb = val_in[Nfunc] * np.cos(val_in[Ncoup] * (np.pi / 180))

                funclist.append(np.array([Kf, Kb, 1, funcvar_holder[i][5]]))

            elif funcvar_holder[i][4] == 42:  # theta functional parameter

                Kf = val_in[Ncoup] * np.sin(val_in[Nfunc] * (np.pi / 180))
                # sets kb
                Kb = val_in[Ncoup] * np.cos(val_in[Nfunc] * (np.pi / 180))

                funclist.append(np.array([Kf, Kb, 1, funcvar_holder[i][5]]))

            else:
                pass

        # hunt for dublicates
        if funcvar[i + 1][2] == 1:  # checks if current varible is coupled
            k = i + 1  # makes sure it doesn't delete itself
            while k != funcvar[0][0] + 1:
                if funcvar[i + 1][3] == funcvar[k][1]:
                    del funcvar[k]  # deletes dublicate
                    funcvar[0][0] -= 1  # fills deleted input

                else:
                    k += 1

        else:
            pass

        i += 1

    return funclist


#decimates the experimental signal to be the same length as input
# JUST WORKS FOR WHEN n(EX) > N(Numb).
def deci_exp_sim(Excurr,Nsim):

    Nex = len(Excurr)

    C = Nex/Nsim

    rc = np.round(C)

    #check to see if multiples of 2
    if np.log2(Nex)%1 == 0 and  np.log2(Nsim)%1 == 0:
        pass
    else:
        print('data is not a multiple of 2')

    curr = Excurr[0::int(round(C))]

    return curr

#stores harmonic generator functions for ease of use
def harm_gen(Curr, time, AC_freq, bandwidth, spaces):
    fft_res = rfft(Curr)

    # this can be done seperatly So lok into it
    freq = rfftfreq(len(fft_res), d=time)  # for ML could be done at the start to stop iterations ALSO DOUBLE VECTOR

    n_window, f_window, g = ifft_windowing(AC_freq, freq, bandwidth,
                                           spaces)  # Extracts the frequency and n windowing {CAN BE SEPERATE}

    hil_store, N_ifft = Harmonic_gen(fft_res, n_window, spaces, bandwidth, g)  # extract envolope data

    return hil_store

# function for finding the nearest value in a series of array
def find_nearest(array, value):
    ''' Find nearest value is an array '''
    idx = (np.abs(array-value)).argmin()
    return idx


# truncates and passes around experimental harmonics and find truncation points
def EXPharmtreatment(EX_hil_store,sigma,Extime, truntime, op_settings):

    nEX = EX_hil_store.shape[1]

    E1 = find_nearest(np.linspace(0, Extime, nEX), truntime[0])
    Ndeci = op_settings[2]

    E1 = int(E1 - E1 % (nEX / Ndeci))  # sets a point for E1 that has Int in simulation file
    Esim1 = E1 * Ndeci / nEX

    if truntime[1] == 'MAX':
        E2 = nEX
        Esim2 = Ndeci
    else:
        E2 = find_nearest(np.linspace(0, Extime, nEX), truntime[1])
        E2 = int(E2 + ((nEX / Ndeci) - E2 % (nEX / Ndeci)))

        Esim2 = E2 * Ndeci / nEX

    Nsimdeci = [int(Esim1), int(Esim2)]
    Nex = [E1, E2]

    # truncates experimental harmonics
    EX_hil_store = EX_hil_store[:, E1:E2]
    # need to identify correct points for simulation and experimental so that its the same time spot then pass identified decimated simulation values to MEC_set

    # Truncates the sigma file if applicable to method
    if isinstance(sigma, bool):  # sigma array exists  """ or method == 'Baye_HarmPerFit'"""
        pass
    else:
        sigma = sigma[:, E1:E2]
        # Somefunction that truncates the sigma values at set time points

    return EX_hil_store, sigma, Nsimdeci, Nex

def log10harmtunc(Excurr, frequens, bandwidth,truntime,AC_freq,Nex):

    print(Nex)
    x = [] # new EXCurr array
    freqx = []
    print(AC_freq)
    print(truntime)
    N = len(bandwidth[0])
    nsimdecionoff = []
    plt.plot(frequens[0:len(Excurr)],np.log10(Excurr))
    plt.savefig('log10truncated12.png')
    plt.close()
    # DC section
    if bandwidth[0][0] != 0:
        if truntime[0] < bandwidth[0][0]:
            #Nlow = find_nearest(frequens, truntime[0])
            Nhigh = -Nex[0] + find_nearest(frequens, bandwidth[0][0])
            x.append(Excurr[0:Nhigh])
            freqx.append(frequens[0:Nhigh])
            nsimdecionoff.append([0,Nhigh])

    j = 0

    while j != len(AC_freq):
        i = 0
        while i != N:
            if (i+1)*AC_freq[j] + bandwidth[j + 1][i] > truntime[1] or (i+1)*AC_freq[j] - bandwidth[j + 1][i]  > truntime[1]:
                if (i + 1) * AC_freq[j] - bandwidth[j + 1][i] < truntime[1]:
                    Nwindlow = -Nex[0] + find_nearest(frequens, (i + 1) * AC_freq[j] - bandwidth[j + 1][i])
                    Nwindhigh = len(Excurr)
                    x.append(Excurr[Nwindlow:Nwindhigh])
                    freqx.append(frequens[Nwindlow:Nwindhigh])
                    nsimdecionoff.append([Nwindlow, Nwindhigh])
                i = N
            else:
                Nwindlow = -Nex[0] + find_nearest(frequens, (i+1)*AC_freq[j] - bandwidth[j + 1][i])
                Nwindhigh = -Nex[0] + find_nearest(frequens, (i+1)*AC_freq[j] + bandwidth[j + 1][i])
                x.append(Excurr[Nwindlow:Nwindhigh])
                freqx.append(frequens[Nwindlow:Nwindhigh])
                nsimdecionoff.append([Nwindlow,Nwindhigh])
                i += 1
        j += 1

    y = np.asarray(x[0])
    yfreq = np.asarray(freqx[0])
    i = 1
    while i != len(x):
        y = np.concatenate((y,x[i]))
        yfreq = np.concatenate((yfreq,freqx[i]))
        i += 1
    Excurr = y

    return Excurr, nsimdecionoff

# truncates and passes around experimental Fourier transformed current and find truncation points
def EXPFTtreatment(FTCurr,Exfreq, truntime,):

    E1 = int(find_nearest(Exfreq, truntime[0]))
    E2 = int(find_nearest(Exfreq, truntime[1]))

    Esim1 = E1
    Esim2 = E2

    Nsimdeci = [int(Esim1), int(Esim2)]
    Nex = [E1, E2]

    # truncates experimental harmonics
    FTCurr = FTCurr[E1:E2]

    return FTCurr, Nsimdeci, Nex

# truncates and passes around experimental harmonics and find truncation points
def EXPtotcurrtreatment(Excurr,Extime, truntime, op_settings):

    nEX = len(Excurr)
    E1 = find_nearest(np.linspace(0, Extime, nEX), truntime[0])
    Ndeci = op_settings[2]

    E1 = int(E1 - E1 % (nEX / Ndeci))  # sets a point for E1 that has Int in simulation file
    Esim1 = E1 * Ndeci / nEX

    if truntime[1] == 'MAX':
        E2 = nEX
        Esim2 = Ndeci
    else:
        E2 = find_nearest(np.linspace(0, Extime, nEX), truntime[1])
        E2 = int(E2 + ((nEX / Ndeci) - E2 % (nEX / Ndeci)))

        Esim2 = E2 * Ndeci / nEX

    Nsimdeci = [int(Esim1), int(Esim2)]
    Nex = [E1, E2]

    Excurr = Excurr[E1:E2]

    return Excurr, Nsimdeci, Nex