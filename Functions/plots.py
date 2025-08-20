from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import numpy as np

CM = 1/2.54  # centimeters in inches

clColorStr = 'Blues'
kColorStr = 'Oranges'
expColorStr = 'Purples'

kColors = getattr(cm, kColorStr)(np.linspace(0.3, 0.8, 4))
clColors = getattr(cm, clColorStr)(np.linspace(0.3, 0.8, 4))
expColors = getattr(cm, expColorStr)(np.linspace(0.3, 0.8, 4))

kColor = kColors[-2]
clColor = clColors[-1]
expColor = expColors[-1]

def plotResults(input, results, plotlocs = ['CA1PYR_TK21[0].soma[0]'], savename = None):
    plt.rcParams['lines.linewidth'] = 0.75
    plt.rc('font', size=8)

    #plots an overview of the results of a simulation
    for recloc in plotlocs:
        loc = str(recloc)
        fig, axs = plt.subplots(4,1 ,figsize = (14.5*CM,8.5*CM)) 
        fig.patch.set_facecolor('white')
        #-------------------------------------------voltages + stimulation timings------------------------------------------------------------
        ax1 = axs[0]
        ax1.set_ylabel('v [mV]', color='black')
        ax1.plot(np.array(results.t)/1e3, np.array(results.recordings.v[loc]), color='black', zorder = 10)
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.set_ylim([-100,60])
        ax1.set_xticks([])
        ax1.set_zorder(3)
        ax1.set_frame_on(False)
        ax1.set_xlim([0,input.duration/1e3])

        ax2 = ax1.twinx()
        ax2.set_ylabel('Optical stim', color='tab:red')

        ax2.plot(np.array(results.t)/1e3, np.array(results.stimulations.Iopto)/np.max(np.abs(np.array(results.stimulations.Iopto))),  color='tab:red', label = 'optical stimulation (' + str(input.stimopt.Ostimparams.amp)+ ' W/m2)')
        ax2.set_xlim([0,input.duration/1e3])
        ax2.set_ylim([-0.05,1.05])
        ax2.set_zorder(1)

        ax3 = ax1.twinx()
        ax3.spines.right.set_position(("axes", 1.2))
        ax3.set_ylabel('Electrical stim', color='grey')
        if np.max(np.abs(np.array(results.stimulations.iStim))) != 0: ax3.plot(np.array(results.t)/1e3, results.stimulations.iStim/np.max(np.abs(np.array(results.stimulations.iStim))), '--', color='grey', label = 'current injection (' + str(round(input.stimopt.Estimparams.amp,2)) + ' nA)', zorder = 1)
        else: ax3.plot(np.array(results.t)/1e3, np.array(results.stimulations.iStim), '--', color='grey', label = 'current injection (' + str(round(input.stimopt.Estimparams.amp,2)) + ' nA)')

        #ax3.legend(loc = 'upper right')
        ax3.set_xlim([0,input.duration/1e3])
        ax3.set_ylim([-0.05,1.05])
        ax3.set_zorder(2)
        ax3.set_xlim([0,input.duration/1e3])

        #-------------------------------------------opsin dynamics------------------------------------------------------------ 
        ax1 = axs[1]
        ax1.set_ylabel('i_' + input.cellsopt.opsin_options.opsinmech + ' \n [mA/cm2]', color= 'black')  
        ax1.plot(np.array(results.t)/1e3, np.array(results.recordings.iopto[loc]), color= 'black')
        ax1.tick_params(axis='y', labelcolor='black')
        if input.axlims: ax1.set_ylim([-.05,0.25])
        ax1.set_xticks([])
        ax1.set_xlim([0,input.duration/1e3])

        ax3 = ax1.twinx()
        ax3.set_ylabel('g [S/cm2]', color='tab:green')
        ax3.plot(np.array(results.t)/1e3, np.array(results.recordings.gopto[loc]),  color='tab:green', label = 'g_opsin')
        ax3.set_xlim([0,input.duration/1e3])

        #------------------------------------------- chloride dynamics ------------------------------------------------------------
        ax1 = axs[2]
        ax1.set_ylabel('E_rev_Cl \n [mV]', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.plot(np.array(results.t[2:])/1e3, np.array(results.recordings.ecl[loc][2:]),  color='black')
        c = max([np.abs(np.max(results.recordings.ecl[loc][2:])-results.recordings.ecl[loc][2]), np.abs(results.recordings.ecl[loc][2] - np.min(results.recordings.ecl[loc][2:]))])
        if not np.isnan(c): ax1.set_ylim([results.recordings.ecl[loc][2]- 1.1*c,results.recordings.ecl[loc][2]+1.2*c])
        if input.axlims: ax1.set_ylim([-90,-30])
        ax1.set_xlim([0,input.duration/1e3])

        ax2 = ax1.twinx()
        ax2.set_ylabel('[Cl]_e [mM]', color=clColors[-1])
        ax2.plot(np.array(results.t[1:])/1e3, np.array(results.recordings.clo[loc][1:]),  color=clColors[-1])
        ax2.tick_params(axis='y')
        c = 1.1*max([np.max(results.recordings.clo[loc][2:])-input.cellsopt.ion_options.cextra0['cl'], input.cellsopt.ion_options.cextra0['cl']-np.min(results.recordings.clo[loc][2:])])
        if not np.isnan(c): ax2.set_ylim([input.cellsopt.ion_options.cextra0['cl']-c,input.cellsopt.ion_options.cextra0['cl']+c])
        if input.axlims: ax2.set_ylim([123,153])
        ax2.set_xlim([0,input.duration/1e3])

        ax3 = ax1.twinx()
        ax3.spines.right.set_position(("axes", 1.2))
        ax3.set_ylabel('[Cl]_i [mM]', color=clColors[0])
        ax3.plot(np.array(results.t[1:])/1e3, np.array(results.recordings.cli[loc][1:]),  color=clColors[0])
        ax3.tick_params(axis='y')
        c = 1.1* max([np.max(results.recordings.cli[loc][2:])-input.cellsopt.ion_options.cintra0['cl'], input.cellsopt.ion_options.cintra0['cl']-np.min(results.recordings.cli[loc][2:])])
        if not np.isnan(c): ax3.set_ylim([input.cellsopt.ion_options.cintra0['cl']-c,input.cellsopt.ion_options.cintra0['cl']+c])

        if input.axlims: ax3.set_ylim([3,33])
        ax1.set_xticks([])
        ax3.set_xlim([0,input.duration/1e3])
        #------------------------------------------- potassium dynamics ------------------------------------------------------------
        ax1 = axs[3]
        ax1.set_ylabel('E_rev_K \n [mV]', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.plot(np.array(results.t[2:])/1e3, results.recordings.ek[loc][2:],  color='black')
        c = max([np.abs(np.max(results.recordings.ek[loc][2:])-results.recordings.ek[loc][2]), np.abs(results.recordings.ek[loc][2] - np.min(results.recordings.ek[loc][2:]))])
        if not np.isnan(c): ax1.set_ylim([results.recordings.ek[loc][2]- 1.1*c,results.recordings.ek[loc][2]+1.2*c])
        if input.axlims: ax1.set_ylim([-90,-30])
        ax1.set_xlim([0,input.duration/1e3])

        ax2 = ax1.twinx()
        ax2.set_ylabel('[K]_e[mM]', color=kColors[-1])
        ax2.plot(np.array(results.t[1:])/1e3, results.recordings.ko[loc][1:],  color=kColors[-1])
        c = 1.1*max([np.max(results.recordings.ko[loc][2:])-input.cellsopt.ion_options.cextra0['k'], input.cellsopt.ion_options.cextra0['k']-np.min(results.recordings.ko[loc][2:])])
        if not np.isnan(c): ax2.set_ylim([input.cellsopt.ion_options.cextra0['k']-c,input.cellsopt.ion_options.cextra0['k']+c])
        if input.axlims: ax2.set_ylim([0,15])
        ax2.set_xlim([0,input.duration/1e3])

        ax3 = ax1.twinx()
        ax3.spines.right.set_position(("axes", 1.2))
        ax3.set_ylabel('[K]_i [mM]', color=kColors[0])
        ax3.plot(np.array(results.t[1:])/1e3, results.recordings.ki[loc][1:],  color=kColors[0])
        ax1.set_xlabel('time [s]')
        c = 1.1* max([np.max(results.recordings.ki[loc][2:])-input.cellsopt.ion_options.cintra0['k'], input.cellsopt.ion_options.cintra0['k']-np.min(results.recordings.ki[loc][2:])])
        if not np.isnan(c): ax3.set_ylim([input.cellsopt.ion_options.cintra0['k']-c,input.cellsopt.ion_options.cintra0['k']+c])
        if input.axlims: ax3.set_ylim([85,100])
        ax3.set_xlim([0,input.duration/1e3])

        plt.tight_layout()
        if savename is not None: plt.savefig(savename)
        plt.show()