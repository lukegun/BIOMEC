# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 06:20:06 2018

@author: luke
"""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



def format_e(n):
    
    a = '%E' % n
    
    return a.split('E')[0].rstrip('0').rstrip('.')[0:4] + 'e' + a.split('E')[1]

def round_sf(number, significant):
    return round(number, significant - len(str(number)))


def name_Allo(var_all,suffall):
    
    name = []
    #n = len(var_all)
    
    i = 0   
    for var in var_all:  # this whole section could be made more efficent but I don't know.
        
        if var_all[i][0] == 11: # resistance
        
            s = "$R_u$ / $%s \Omega$" %suffall[i]
            name.append(s)
        
        elif var_all[i][0] == 12: # kinomatic viscosity
            name.append('Kinematic_Viscosity ($cm^2$/s)')
            
        elif var_all[i][0] == 13: # kinomatic viscosity
            s = '$\sigma$ / $%sA$' %suffall[i]
            name.append(s)
        
        elif var_all[i][0] == 21: # concentration
            s = '$C%s / $%sM$' %(naminguderscore(var),suffall[i])  # FIx
            name.append(s)
            
        elif var_all[i][0] == 22: # diffusion
            s = '$D%s / $cm^2\;s^{-1}$' %(naminguderscore(var))
            name.append(s)
            
        elif var_all[i][0] == 31: # reaction rate forward
            s = '$k_{f}%s' %(naminguderscore(var))
            name.append(s)
                        
        elif var_all[i][0] == 32: # reaction rate backward
            s = '$k_{b}%s$' %(naminguderscore(var))
            
        elif var_all[i][0] == 33: # Reaction potential Eo
            s = '$E^{0}%s / $%sV$' %(naminguderscore(var),suffall[i])
            name.append(s)
            
        elif var_all[i][0] == 34: # electon kinetic reaction rate
            s = '$k^{0}%s / $cm\;s^{-1}$' % (naminguderscore(var))
            name.append(s)
            
        elif var_all[i][0] == 35: # alp or lambda
            s = '$\\alpha%s' % (naminguderscore(var))
            name.append(s)
            
        elif var_all[i][0] == 41 or var_all[i][0] == 42: # alp or lambda
            s = '$K%s' % (naminguderscore(var))
            name.append(s)
            
        elif var_all[i][0] == 51: # alp or lambda
            s = '$C^0_{dl}$ / $%sF\;cm^{-2}$' % (suffall[i])
            name.append(s)
            
        elif var_all[i][0] == 56: # alp or lambda
            name.append('S')
            
        else:
            print('error in allocation module')
        
        i += 1
    
    
    return name


def unitsuffix(meanvaluesin,var_all):
    
    scalarall = []
    suffall = []
    
    # for negitive values don't screw up the scaling
    meanvalues = [abs(y) for y in meanvaluesin]
    
    i = 0
    for var in var_all:
        if var[0] == 22 or var[0] == 56 or var[0] == 34 or var[0] == 35 or var[0] == 41 or var[0] == 42:    #exception for D and k0
            unit = ""
            scalar = 1
            
        elif meanvalues[i] > 1 and meanvalues[i] < 10**(3):
            unit = ""
            scalar = 1

        elif meanvalues[i] > 10**(-3) and meanvalues[i] < 1:
            unit = "m"
            scalar = 10**(3)

        elif meanvalues[i] > 10**(-6) and meanvalues[i] < 10**(-3):
            unit = "\mu "
            scalar = 10 ** (6)
            
        elif meanvalues[i] > 10**(-9) and meanvalues[i] < 10**(-6):
            unit = "n"
            scalar = 10 ** (9)

        elif meanvalues[i] > 10 ** (-12) and meanvalues[i] < 10 ** (-9):
            unit = "p"
            scalar = 10 ** (12)

        else:
            print("to small a current for the plot axises")
            exit()
        
        scalarall.append(scalar)
        suffall.append(unit)
        i += 1
    
    
    return scalarall, suffall

def naminguderscore(var):
    
    if var[1] == 0:
        s = "$"
    
    elif var[1] == 0 and (var[1] == 31 or var[1] == 32):
        s = "^{%i}$" %var[1]
    else:
        s = "_{%i}$" %var[1]
    
    return s

def KDP_plotter(i,df,name,outputfname):
    
    font = 14
    
    #plot 2D
    j = i + 1
    while j != len(df):
            
        # sets up the print varibles
        cov = np.cov(df[i],df[j])
        sig1 = np.sqrt(cov[0,0])
        sig2 = np.sqrt(cov[1,1])
        rho = cov[0,1]/(sig1*sig2)
        
        textstr = '$\sigma_x$ = %s\n$\sigma_y$ = %s\n$\\rho_{xy}$ = %.3f' %(format_e(sig1),format_e(sig2),rho)
            
        # plot 2d hist
        plt.figure()
        ax2 = sns.kdeplot(df[i], df[j],cmap="Blues",n_levels=10,shade=True)
        plt.scatter(np.mean(df[i]),np.mean(df[j]),color='k',marker = 'x')
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
    
    
    return

def hist_plotter(i,df,name,outputfname):
    
    font = 14
    
    plt.figure()

    ax1 = sns.distplot(df[i], kde=False, hist_kws={'weights': np.full(len(df[0]), 1/len(df[0])), 'color':"blue"})
    plt.axvline(np.mean(df[i]), color='k', linestyle='dashed', linewidth=1)
    plt.axvline(np.mean(df[i]) + 2*np.std(df[i]), color='r', linestyle=':', linewidth=1)
    plt.axvline(np.mean(df[i]) - 2*np.std(df[i]), color='r', linestyle=':', linewidth=1)
    ax1.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    plt.locator_params(axis='x', nbins=7)   # sets number of ticks
    plt.ylabel('Density',fontdict={'size':font})
    plt.xlabel(name[i],fontdict={'size':font})
        
    # Save as png
    s = '%s/histVar_%i' %(outputfname,i)
    plt.savefig( s, bbox_inches='tight')
    plt.close()
    
    return