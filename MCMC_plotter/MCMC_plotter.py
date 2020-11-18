# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 06:19:43 2018

-Fixed 41 an 42 exception

@author: luke
"""

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from MCMC_plotter_modules import *
import os
import numpy as np
import pandas as pd
# load data
# load MCMC density

#plt.figure(2)
#plot density plots
#plt.ticklabel_format(style='plain', axis='y'))
#sns.jointplot(df[:,4], df[:,2],kind="kde",norm_hist=False,cbar=True)

# define file inputs
filename = 'Iter_log.txt'
outputfname = 'Final_paper_post'
# some other stuff if needed

histplot = True
burnin = 0.25
font = 14

df = pd.read_csv(filename,sep='\t',index_col=False,header = None)

df = np.array(df.values)
df = df[int(burnin*len(df[:,0]))::,:-1]   # removes the nan value


# input varible setttings (code,repeat)
var_all = [[11,1],[33,0],[35,0],[56,0]]

# functional varible constant ( put in in sequentual order) [is it pair (0 = no, 0 = yes), constant]
[[0,],[]]

# create files
os.makedirs(outputfname)

# scaling to standard units
n = len(var_all)
i = 0
while i != n:
    if var_all[i][0] == 21: # concentration
        df[:,i] = 1000*df[:,i]  # scales to mols per litre
    
    # exemption for polar coordintes
    elif var_all[i][0] == 41:   #R
        pass # can't scale but need something to show magnitude
    elif var_all[i][0] == 42:   #theta
        df[:,i] = np.tan(df[:,i]*np.pi/180) 
        
    else:
        pass
    i += 1
    
# something to calculate the means
meanvalues = []
i = 0
for x in var_all:
    me = df[:,i].mean()
    meanvalues.append(me)
    i += 1
    
# something to calc scaled props    
scalarall, suffall = unitsuffix(meanvalues,var_all)

# something to calculate the scalled values
i = 0
for x in var_all:
    df[:,i] = df[:,i]*scalarall[i]
    i += 1


#Allocates name
name = name_Allo(var_all,suffall)

if histplot:
    i = 0
    plt.rc('xtick',labelsize=font)
    plt.rc('ytick',labelsize=font)
    for varibles in var_all:
    
        #plot stuff 1D
        me = df[:,i].mean()
        
        plt.figure()

        ax1 = sns.distplot(df[:,i], kde=False, hist_kws={'weights': np.full(len(df[:,0]), 1/len(df[:,0])), 'color':"blue"})
        plt.axvline(me, color='k', linestyle='dashed', linewidth=1)
        plt.axvline(me + 2*np.std(df[:,i]), color='r', linestyle=':', linewidth=1)
        plt.axvline(me - 2*np.std(df[:,i]), color='r', linestyle=':', linewidth=1)
        ax1.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
        ax1.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
        plt.locator_params(axis='x', nbins=7)   # sets number of ticks
        plt.ylabel('Density',fontdict={'size':font})
        plt.xlabel(name[i],fontdict={'size':font})
        
        # Save as png
        s = '%s/histVar_%i' %(outputfname,i)
        plt.savefig( s, bbox_inches='tight')
        plt.close()
        
        #plot 2D
        j = i + 1
        while j != len(var_all):
            
            # sets up the print varibles
            cov = np.cov(df[:,i],df[:,j])
            sig1 = np.sqrt(cov[0,0])
            sig2 = np.sqrt(cov[1,1])
            rho = cov[0,1]/(sig1*sig2)
        
            textstr = '$\sigma_x$ = %s\n$\sigma_y$ = %s\n$\\rho_{xy}$ = %.3f' %(format_e(sig1),format_e(sig2),rho)
            
            # plot 2d hist
            plt.figure()
            ax2 = sns.kdeplot(df[:,i], df[:,j],cmap="Blues",n_levels=10,shade=True)
            plt.scatter(me,df[:,j].mean(),color='k',marker = 'x')
            ax2.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
            ax2.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
            plt.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=font,
                     verticalalignment='top',size=font,bbox=dict(facecolor='white', alpha=0.5,boxstyle="square"))
            plt.locator_params(axis='x', nbins=7)   # sets number of ticks
            plt.locator_params(axis='y', nbins=7)   # sets number of ticks
            plt.ylabel(name[j],fontdict={'size':font})
            plt.xlabel(name[i],fontdict={'size':font})
            
            s = '%s/KDE_%i-%i' %(outputfname,i,j)
            plt.savefig(s, bbox_inches='tight') # Saves the image
            plt.close()
            # save 
            j += 1
        
        i += 1
# move to a seperate section
#elif compairplot:
    #pass
    
