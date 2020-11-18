# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 10:11:41 2020

@author: luke

"""

from inputwrittermod import *
import os
import sys

# inputs seperators
breaker1 = "!Settings   (value,min,max,log sensitivity, parameter code, repeat line)"
breaker2 = "MECSim settings   ! DO NOT REMOVE THESE SEPERATORS"
breaker3 = "Experimental settings ! {NameEX.txt} Repet for triplicates effects method Important do not remove this line" 


#################### PART THAT TELLS THE CODE WHAT TO DO #####################

# header settings
tot = "STANDARD"
# sets up optimisation method
correctinput = False
while not correctinput:
    optimisation = str(input("Please input optimisation method (CMAES or ADMCMC): "))
    if optimisation == "CMAES" or  optimisation == "ADMCMC":
        correctinput = True
    else:
        print("Incorrect optimisation input")
        
correctinput = False
while not correctinput:
    print("List of logic methods:\nHarmPerFit\nBaye_HarmPerFit\nTCDS\nFTC\nLog10FTC ")
    logicm = str(input("Please input logic method : "))
    if logicm == "HarmPerFit" or  logicm == "Baye_HarmPerFit" or logicm == "TCDS" or logicm == "Log10FTC" or logicm == "FTC":
        correctinput = True
    else:
        print("Incorrect optimisation input, please input exactly from above list")

header = [tot, logicm, optimisation]


####################### Start OF OPT SETTINGS ################################

gensetting = optsettings(header)

# custom optimisation settings
if optimisation == "CMAES":
    outputlist = CMAsettings(header)
elif optimisation == "ADMCMC":
    outputlist = ADMCMCsettings(header) 
else:
    print("ERROR in assigning header settings, TERMINATING PROGRAM")
    sys.exit()
    
file_startpart = gensetting + outputlist

####################### PART THAT TAKES IN MECFILE AND EXP ###################

# Takes the input file for the simulation
mecreal = False
while not mecreal:
    MECfilename = str(input("Please input MECSim input file name: "))
    # check to see if the MECsim file is real
    if os.path.exists(MECfilename):
        mecreal = True
        with open(MECfilename) as f:
            MECfile = f.read().splitlines() 
    else:
        print("MECSim file was incorrect and could not be found, try again.")
        

Nexp = int(input("Number of experimental files you are comparing the simulation to: "))
print("put something here so it works with logic FIX")
experimentalfile = []
print("Please put the file input location relative to \nwhere you're running MECSim analytics package\neg. NameEX.txt or file/NameEX.txt ")
for i in range(Nexp):
    x = str(input("Experimental filename " + str(i +1) + ": "))
    experimentalfile.append(x)

# sets up for the bottom part
file_endpart =  ["\n"] + [breaker2] + MECfile + ["\n"] +[breaker3] + experimentalfile

filetot = [breaker1] + file_startpart + file_endpart


######################## PART THAT WRITES THE FILE ###########################
print("settings have been input succsessfully")
filename = str(input("Please input name of text file these settings will be written to? (include .txt at the end): "))

filewritter(filename, filetot)

print("process has been completed with settings saved to text file " + filename)

