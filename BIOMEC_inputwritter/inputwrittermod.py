# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 10:11:50 2020

@author: luke
"""

def optsettings(header):
    
    Npara = int(input("Please input number of parameters to optimise for: "))
    
    paralist = paraopt(header,Npara)
    
    scalfunclist = scalarfuncvaribles(header,Npara,paralist)
    
    if header[1] == "TCDS":

        harmfilter = ["0 ! windowing function (Zero for square, 1 for CG windowing;std = 0.1 Hz recomended)"]
        
        bDCside = " ! Bandwidth Fundimental (BANDWIDTH FIRST)"
        wDCside = " ! Fundimental weights"
        bharmside = " ! Harmonic Bandwidth REPEAT (lowest freq to highest)"
        wharmside = " ! Harmonic Weights REPEAT (lowest freq to highest)"
        
        harmsetting = ["1,0" + bDCside,"0,0" + bharmside,"1,0" + wDCside,"0,0" + wharmside]
    else:

        derp = " ! windowing function (Zero for square, 1 for CG windowing;std = 0.1 Hz recomended)"
        correct = True
        print("Please input windowing function for harmonics\n0 = square window (Simple but more ringing)\n1 = Guassian Convolutional (Simple but more ringing)\n")
        while correct:
            logicm = int(input("Please input windowing function number : "))
            if logicm == 0 or logicm == 1:
                correct = False
                if logicm == 1:
                    fl = float(input("Please input guassian std (0.1 strongly recommended) : "))
                    x = '%i, %.4f' %(logicm, fl)
                else:
                    x = str(logicm)
            else:
                print("incorrect input try again")

        harmfilter = [x + derp]

        harmsetting = harmonichandler()
    
    # Takes in logic for frequency domain shit 
    trunsettings = timetrun(header)

    s = "%s %s %s ! Header declaration (tot Logic Method)" %(header[0],header[1], header[2])
    headertot = [s]
    
    generaloutputs = headertot + paralist + scalfunclist + harmfilter + harmsetting + trunsettings
    
    return generaloutputs

def paraopt(header,Npara):
    
    print("\nList of parameters and relavent code:\n")
    print("Resistance: 11\nKinematic Viscosity: 12\nExperimental noise (TCDS only): 13\n"\
          + "Concentration: 21\n Diffussion coefficent: 22\nForward Reaction rate: 31\n"\
          + "Backward Reaction rate: 32\nFormal Potential: 33\nElectron Transfer Rate: 34\n"\
          + "Alpha or Lambda: 35\nEquilibrium Constant magnitude (func): 41\nEquilibrium Constant Theta (func): 42\n"\
          + "Capacitance constants C0-C4: 51-55\nCapacitance Scalar Multiple: 56\n\n")
    
    s = "%s\t\t! Number of varibles" %int(Npara)
    paramterlist = [s]
    for i in range(Npara):
        s = "\nParameter #%i" %(i+1)
        print(s)
        
        code = int(input("Please give code of parameter you want to optimise: "))
        if code>20 and code<45:
            repeat = int(input("Please input which repeat line the parameter \n"\
                               "is on in MECSim file (Start at 1): "))
        else:
            repeat = 1
            
        if header[2] == "CMAES":
            logs = int(input("is this parameter optimised in a log scale (yes = 1 or no = 0): "))
        logs = 0
        correctinput = False
        while not correctinput:
            xmed = float(input("Please input median parameter value: "))
            xsmall = float(input("Please input smallest parameter range: "))
            xlarge = float(input("Please input largest parameter range: "))
            if xsmall < xmed and xmed < xlarge:
                correctinput = True
            else:
                print("\n Parameters range is not correct, please redo order of small, med large is correct")
        
        xlist = "%s, %s, %s, %i, %i, %i" %(format_e(xmed),format_e(xsmall),format_e(xlarge),logs,code,repeat)
        
        paramterlist.append(xlist)
        
    return paramterlist

def scalarfuncvaribles(header,Npara,paramterlist):
    
    Ns = int(input("Please input number of scaled varibles: "))
    s = "%i\t\t! number of scaled varibles" %Ns
    Nscallist = [s]
    
    side = " \t! scaling reletive varibles (scaling factor, varible number, parameter code, repeat line)"
    for i in range(Ns):
        Npara = int(input("Please input which parameter is being scaled in the input parameter (Starts at one): "))
        scal = float(input("Please input the scalar value: "))
        Scode = int(input("Please input output parameter varible code: "))
        Srepeat = int(input("Please input output parameter reapeat line in MECSim: "))
     
        s = "%.4f, %i, %i, %i" %(scal, Npara, Scode,Srepeat)
        Nscallist.append(s + side)
        
    # functional varibles settings
    Nf = int(input("Please input number of fuctional varibles: "))
    s = "%i\t\t! Number of functional scaling" %Nf
    Nscallist.append(s)
    
    side = " \t! functional parameter (function, varible number, paired value scalar, {paired value, varible row})"
    for i in range(Nf):
        Npara = int(input("Please input function: "))
        scal = int(input("Please input which parameter is being used in function (Starts at one): "))
        Scode = int(input("Please input paired value scalar (use 0 if nothing): "))
        Srepeat = float(input("Plese input paired value or paired varible input row: "))
        
        if Scode == 1:
            s = "%i, %i, %i, %i" %(Npara, scal, Scode,int(Srepeat))
        else:
            s = "%i, %i, %i, %s" %(Npara, scal, Scode,format_e(Srepeat))
        Nscallist.append(s + side)
        
        
    return Nscallist

def harmonichandler():    
    
    Nac = int(input("Please input number of AC frequencies used in experimental data: "))
    Nharm = int(input("Please input number of harmonics present (excluding DC): "))
    
    DCband = float(input("Please input bandwidth of DC component: "))
    DCwieght = float(input("Please input scalar weight of DC component: "))
    
    
    # sets up for the DC somponent
    DCbandhold = "%.2f" %DCband
    DCwieghthold = "%.2f" %DCwieght
    
    for i in range(Nharm - 1):
        DCbandhold += ", 0"
        DCwieghthold += ", 0"
        
    # 
    if Nac > 1:
        print("As you are using more then one AC frequency note that\n" +\
              "you will be asked to input multple sets of data with the\n" +\
              "these will be analyised by the code from LOWEST FREQUENCY TO "+\
              "HIGHEST FREQUENCY.")
        
    tothbandhold = []
    tothwieghthold = [] 
    
    for i in range(Nac):
        harmband = float(input("Please input bandwidth of fundimental harmonic: "))
        harmwieght = float(input("Please input scalar weight of fundimental harmonic: "))
        
        harmbandhold = "%.2f" %harmband
        harmwieghthold = "%.2f" %harmwieght
        
        correct_inp = False
        while not correct_inp:
        
            same = int(input("Are the harmonic weights and bandwidths the same for all harmonics? (0 = No; 1 = Yes): "))
            if same == 1:
                for j in range(Nharm-1):
                    harmbandhold += ", %.2f" %harmband
                    harmwieghthold += ", %.2f" %harmwieght
                    
                    correct_inp = True
                        
            elif same == 0:
                for j in range(Nharm - 1):
                    hb = float(input("Please input bandwidth of the #%i harmonic: " %(j+1) ))
                    hw = float(input("Please input scalar weight of the #%i harmonic: " %(j+1)))
                    
                    harmbandhold += ", %.2f" %hb
                    harmwieghthold += ", %.2f" %hw
                    
                    correct_inp = True    
            else:
                print("Incorrect input on same")
                
        tothbandhold.append(harmbandhold) 
        tothwieghthold.append(harmwieghthold) 
        
    #sets the side stings
    bDCside = " ! Bandwidth Fundimental (BANDWIDTH FIRST)"
    wDCside = " ! Fundimental weights"
    bharmside = " ! Harmonic Bandwidth REPEAT (lowest freq to highest)"
    wharmside = " ! Harmonic Weights REPEAT (lowest freq to highest)"
        
    # something to put this all in as an output
    # band first
    harmsetting = [DCbandhold + bDCside]
    for i in range(Nac):
        harmsetting.append(tothbandhold[i] + bharmside)
    
    # harmoni
    harmsetting.append(DCwieghthold + wDCside)
    for i in range(Nac):
        harmsetting.append(tothwieghthold[i] + wharmside)
        
        
    return harmsetting

# does truncation time with exception to frequency domain plus junk parameters
def timetrun(logic):
    
    timeend = "\t\t\t! Truncation points (sec) (0,MAX; for not applicable)"
    
    correct_inp = False
    while not correct_inp:
        trun = int(input("Are you truncating the time series for analysis? (0 = No, 1 = Yes): "))
        # truncation time
        if trun == 1:
            if logic == "Log10FTC" or logic == "FTC":
                mintime = float(input("Please input the minimum truncation Frequency (Hz): "))
                maxtime = float(input("Please input the maximum truncation Frequency (Hz): "))
            else:
                mintime = float(input("Please input the minimum truncation time (sec): "))
                maxtime = float(input("Please input the maximum truncation time (sec): "))
            
            s = "%.3f,%.3f" %(mintime,maxtime)
            correct_inp = True
            
        elif trun == 0:
            s = "0,MAX"               
            correct_inp = True
        
        else:
            print("ERROR: Incorrect previous input, retry")
        
    trunsettings = [s + timeend]
    
    # Fitting method used
    s = "1\t\t\t! Fitting method used  (0 = absdiff, 1 = %diff)"
    trunsettings.append(s)
    
    # Experimental input type
    s = "1\t\t\t! Experimental input type (0 = MECSim, 1 = FTACV, 2=CHI)"
    trunsettings.append(s)
    
    return trunsettings
    

def CMAsettings(header):
    
    settingslist = []
    
    # Number of multiprocesses
    print("\nThe calculation is generally run accross multiple CPU processors to speed up calculations")
    print("Warning for CMAES do not exceed int(3*log(#parameters)")
    Ncore = int(input("Please input number of CPU processors: "))
    s = "%i\t\t\t! Number of cores to be used" %Ncore
    settingslist.append(s)
    
    # 2^N comparison points
    print("\nTo speed up calculations only certain number of experimental data points are optimised of order of 2^N")
    Ndata = int(input("Please input N in above statement: "))
    s = "%i\t\t\t! 2^N comparison data points per current" %Ndata
    settingslist.append(s)
    
    # tolx IMPORTANT TO EXPLAIN
    print("\nTolx is the termination criterion for CMAES, it gives the level of of accuracy till it stops (low = 0.05, mid = 0.025, high = 0.01)")
    tolx = float(input("Please input Tolx in above statement: "))
    s = "%.4f\t\t\t! tolx, value of x as %%*range/2 needed before fit is meet (0.05,0.025,0.01 recomended)"  %tolx
    settingslist.append(s)
    
    # initial sigma value as %*range (0.33 recomendanded)
    print("\nIntial sigma for CMA-ES is the starting range as a ratio of starting scan range")
    sigma = float(input("Please input sigma in above statement (0.33 is strongly recomended): "))
    s = "%.4f\t\t\t! initial sigma value as %%*range (0.33 recomendanded)"  %sigma
    settingslist.append(s)
    
    return settingslist

def ADMCMCsettings(header):
    
    settingslist = []
    
    # 2^N comparison points
    print("\nTo speed up calculations only certain number of experimental data points are optimised of order of 2^N")
    Ndata = int(input("Please input N in above statement: "))
    
    # number of overall chains
    Nchain =str(input("\nPlease input number of MCMC chains desired (soft limit 4): "))
    
    Ncore = Nchain
    
    s = "%s\t\t\t! Number of cores to be used" %Ncore
    settingslist.append(s)
    s = "%s\t\t\t! 2^N comparison data points per current" %Ndata
    settingslist.append(s)  
    
    # prior sampling
    print("\nIn ADMCMC prior sampling is inital run as part of the algorithm.")
    Nprior = int(input("Please put initial prior sampling per varible (100-250 recomended): "))
    s = "%i\t\t\t! MCMC initial prior sampling per varible" %Nprior
    settingslist.append(s)
    
    # Number of trails in overall chain
    print("\nThis is number of simulations for the entirly, remember that each chain = totchain/#")
    chaintot = int(input("Please give the number of total simulations to run overall: "))
    s = "%i\t\t\t! number of trail for overall chain" %chaintot
    settingslist.append(s)
    
    # burnin period
    print("\nAt the end of the MCMC calculation a percentage of the trails are thrown away, this is refered to as burnin")
    burnin = float(input("Please put in the burnin % (%/100 so between 0-1): "))
    s = "%.4f\t\t\t! burnin period as ratio of chain legth (%%/100)" %burnin
    settingslist.append(s)
    
    # noise (mostly obsolte)
    sigma = "1.0e-6\t\t\t! noise (%%/100) for each frequency (lowest freq to highest)"
    settingslist.append(sigma)
    
    s = Nchain + "\t\t\t! number of chains to run" 
    settingslist.append(s)
    
    return settingslist

def filewritter(filename, filetot):
    
     f = open(filename,"w")
     for lines in filetot:
         if lines != "\n":
             f.write(lines + "\n")
         else:
             f.write(lines)
             
def format_e(n):
    
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'e' + a.split('E')[1]


