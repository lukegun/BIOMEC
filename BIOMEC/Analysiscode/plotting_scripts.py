"""
PINTS MECSim plotting functions





"""


import numpy as np
import ML_signal_processing as MLsp
import Script_generator as Scripgen
import matplotlib.pyplot as plt
import os
import seaborn as sns


#extracts the fourier transform for plotting


def logten_extractor(curr,Extime, truntime,bandwidth,AC_freq):
    # Translates the current into the frequency domain
    curr = np.fft.rfft(curr)
    curr = curr / len(curr)  # adjustment for amplitude spectrum as fft does not come out that way
    frequens = np.fft.rfftfreq(len(curr), d=Extime / len(curr))
    curr, Nsimdeci, Nex = MLsp.EXPFTtreatment(curr, frequens, truntime)

    # settings to remove background issue in log10fit

    low = MLsp.find_nearest(frequens, truntime[0])
    high = MLsp.find_nearest(frequens, truntime[1])

    curr = np.log10(curr[low:high]).real
    freq = np.linspace(truntime[0],truntime[1],len(curr))

    return curr, freq

def harmplot(outputfname,Simcurr,Simtime,EX_hil_store,Exp_t,bandwidth,AC_freq,spaces):

    #EX_hil_store = MLsp.harm_gen(EXcurr, Exp_t[1], AC_freq, bandwidth, spaces) Not needed as passed in straight

    Sim_hil_store = MLsp.harm_gen(Simcurr, Simtime[1], AC_freq, bandwidth, spaces)

    os.makedirs(outputfname)

    N = len(bandwidth[0])

    if EX_hil_store.shape[1] > Sim_hil_store.shape[1]:  # dessimate experimental current

        deci = int(EX_hil_store.shape[1] / Sim_hil_store.shape[1])
        EX_hil_store = EX_hil_store[:, ::deci]
        Exp_t = Exp_t[::deci]

    elif EX_hil_store.shape[1] < Sim_hil_store.shape[1]: # dissimate simulated current

        deci = int(Sim_hil_store.shape[1]/EX_hil_store.shape[1])
        Sim_hil_store = Sim_hil_store[:, ::deci]
        Simtime = Simtime[::deci]

    else:
        pass

    pererror = []

    # plot harmonics
    i = 0
    while i != N + 1:


        # calculate percentage error
        x1 = sum((EX_hil_store[i, :] - Sim_hil_store[i, :])**2)/sum(EX_hil_store[i, :]**2)
        x1 = np.sqrt(x1)

        pererror.append(x1)

        if i == 0:
            textstr = 'DC Signal\n%%Err = %.3f' % (x1 * 100)
        elif i == 1:
            textstr = '1st Harmonic\n%%Err = %.3f ' % (x1 * 100)
        elif i == 2:
            textstr = '2nd Harmonic\n%%Err = %.3f  ' % (x1 * 100)
        elif i == 3:
            textstr = '3rd Harmonic\n%%Err = %.3f  ' % (x1 * 100)
        else:
            textstr = '%ith Harmonic\n%%Err = %.3f  ' % (i,x1 * 100)
        plt.figure()
        fig, ax = plt.subplots()
        plt.plot(Exp_t,EX_hil_store[i,:],color='k',label='Experimental')
        plt.plot(Simtime, Sim_hil_store[i, :], color='r',label='Simulated',linestyle='-.')
        plt.ylabel('Current (Amps)')
        plt.xlabel('Time (sec)')
        plt.text(0.40, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                 verticalalignment='top', size=10, bbox=dict(facecolor='white', alpha=0.5, boxstyle="square"))
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        ax.legend(loc = 'upper right')
        s = '%s/Harmonic%i' % (outputfname, i)
        plt.savefig(s, bbox_inches='tight')
        plt.close()

        # error as function of time
        i += 1

    return

# allocates name as varibles
def name_Allo(var_all):
    name = []
    n = len(var_all)

    i = 0
    while i != n:  # this whole section could be made more efficent but I don't know.

        if var_all[i][0] == 11:  # resistance
            name.append('Resistance ($\Omega$)')

        elif var_all[i][0] == 12:  # kinomatic viscosity
            name.append('Kinematic_Viscosity ($cm^2$ $s^{-1}$)')

        elif var_all[i][0] == 13:
            name.append('Experimental Noise (A)')

        elif var_all[i][0] == 21:  # concentration
            name.append('Concentration %i (mol $L^{-1}$)' % (var_all[i][1]))  # FIx

        elif var_all[i][0] == 22:  # diffusion
            name.append('Diffusion %i ($cm^2$ $s^{-1}$)' % (var_all[i][1]))

        elif var_all[i][0] == 31:  # reaction rate forward
            name.append('Reaction Rate Forward %i' % (var_all[i][1]))

        elif var_all[i][0] == 32:  # reaction rate backward
            name.append('Reaction Rate Backward %i' % (var_all[i][1]))

        elif var_all[i][0] == 33:  # Reaction potential Eo
            name.append('Formal Potential %i (V)' % (var_all[i][1]))

        elif var_all[i][0] == 34:  # electon kinetic reaction rate
            name.append('ET rate constant %i (cm $s^{-1}$)' % (var_all[i][1]))

        elif var_all[i][0] == 35:  # alp or lambda
            name.append('Alpha %i' % (var_all[i][1]))

        elif var_all[i][0] == 41 or var_all[i][0] == 42:  # alp or lambda
            name.append('K%iapp' % (var_all[i][1]))

        elif var_all[i][0] >= 51 and var_all[i][0] <= 55:  # alp or lambda
            name.append('Cap%i (F $cm^$-2$)' % (var_all[i][1] - 1))

        elif var_all[i][0] == 56:  # alp or lambda
            name.append('Capacity Scalar')

        else:
            print('error in allocation module')

        i += 1

    return name


# plot the intergration
def density_plotter(filename, df,var_all,burnin):
    font = 14

    df = df[int(burnin * len(df[:, 0]))::, :]  # removes the nan value

    name = name_Allo(var_all) # gets names for the distrabution

    # makes the file to save the pics in
    os.makedirs(filename)

    # scaling to standard units
    n = len(var_all)
    i = 0
    while i != n:
        if var_all[i][0] == 21:  # concentration
            df[:, i] = 1000 * df[:, i]  # scales to mols per litre

        # exemption for polar coordintes
        elif var_all[i][0] == 41:  # R
            pass  # can't scale but need something to show magnitude
        elif var_all[i][0] == 42:  # theta
            df[:, i] = np.tan(df[:, i] * np.pi / 180)

        else:
            pass
        i += 1

    if True:
        i = 0
        plt.rc('xtick', labelsize=font)
        plt.rc('ytick', labelsize=font)
        for varibles in var_all:

            # plot stuff 1D
            me = df[:, i].mean()

            plt.figure()

            ax1 = sns.distplot(df[:, i], kde=False,
                               hist_kws={'weights': np.full(len(df[:, 0]), 1 / len(df[:, 0])), 'color': "blue"})
            plt.axvline(me, color='k', linestyle='dashed', linewidth=1)
            plt.axvline(me + 2 * np.std(df[:, i], ddof=1), color='r', linestyle=':', linewidth=1)
            plt.axvline(me - 2 * np.std(df[:, i], ddof=1), color='r', linestyle=':', linewidth=1)
            ax1.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
            plt.locator_params(axis='x', nbins=7)  # sets number of ticks
            plt.ylabel('Probability Density', fontdict={'size': font})
            plt.xlabel(name[i], fontdict={'size': font})

            # Save as png
            s = '%s/histVar_%i' % (filename, i)
            plt.savefig(s, bbox_inches='tight')
            plt.close()

            # plot 2D
            j = i + 1
            while j != len(var_all):
                # sets up the print varibles
                cov = np.cov(df[:, i], df[:, j], ddof=1)
                sig1 = np.sqrt(cov[0, 0])
                sig2 = np.sqrt(cov[1, 1])
                rho = cov[0, 1] / (sig1 * sig2)

                textstr = '$\sigma_x$ = %s\n$\sigma_y$ = %s\n$\\rho_{xy}$ = %.3f' % (
                Scripgen.format_e(sig1), Scripgen.format_e(sig2), rho)

                # plot 2d hist
                plt.figure()
                ax2 = sns.kdeplot(df[:, i], df[:, j], cmap="Blues", n_levels=10, shade=True)
                plt.scatter(me, df[:, j].mean(), color='k', marker='x')
                ax2.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
                ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
                plt.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=font,
                         verticalalignment='top', size=font, bbox=dict(facecolor='white', alpha=0.5, boxstyle="square"))
                plt.locator_params(axis='x', nbins=7)  # sets number of ticks
                plt.locator_params(axis='y', nbins=7)  # sets number of ticks
                plt.ylabel(name[j], fontdict={'size': font})
                plt.xlabel(name[i], fontdict={'size': font})

                s = '%s/KDE_%i-%i' % (filename, i, j)
                plt.savefig(s, bbox_inches='tight')  # Saves the image
                plt.close()
                # save
                j += 1

            i += 1

    return name

# singular convergence plotter
def sinularconverg(filename,dist,variblenames):

    i = 0
    for names in variblenames:
        x = []
        for j in range(len(dist[:])):
            x.append(dist[j][i])
        plt.figure()
        plt.plot(x)
        plt.xlabel(names)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        s = filename + '/Var'+str(i)
        plt.savefig(s,bbox_inches='tight')
        plt.close()
        i += 1

    return


#multiple convegence plotter
def multiconverg(filename, Mchains, variblenames):

    i = 0
    for names in variblenames:
        plt.figure()
        j = 0
        while j != len(Mchains):
            x = []
            for l in range(len(Mchains[1][:])):
                x.append(Mchains[j][l][i])
            plt.plot(x)
            j += 1
        plt.xlabel(names)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        s = filename + '/Var'+str(i)
        plt.savefig(s,bbox_inches='tight')
        plt.close()

        i += 1

    return