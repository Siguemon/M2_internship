import sys ; this_filename = sys.argv[0].split('.')[0]
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def get_fidcrop(key):
    '''
    get the fidcrop from the format of the key in the data (.npz file output from the make_total_spectra package)
    <key> must be in fmt "XX_<fidcrop>mm"
    returns the fidcrop in mm
    '''
    return float(key.split('_')[1][:-2])

def get_fiducial_volume(fidcrop):
    '''
    calculates the fiducial volume value from the fidcrop
    <fidcrop> in mm
    '''
    height = 10.64 #mm
    radius = 15.368 #mm

    return (height - fidcrop)/height * (radius - fidcrop)**2/(radius)**2 #%

def plot_spectra(ax, data, err_data, err_kw={}, color='k', ls='-'):
    '''
    Plots the spectra from given dataset, over a specified axis of the
    matplotlib figure and for given bins. The conversion counts->DRU is made
    within the function.

    * Inputs :
    <ax> the matplotlib's figure axis
    <data> the energies vector you want to display the spectra of
    <err_data> the vertical errors associated to the data
    <err_kw> dictionnary for the errorbars format
    <color> the color of the spectra and errorbars
    <ls> the linestyle (matplotlib's parameter)
    '''
    ax.step((bins[1:]+bins[:-1])/2, data, color=color, ls=ls, where='mid')
    ax.errorbar(.5*(bins[1:]+bins[:-1]), data, yerr=err_data, **err_kw, color=color)

f = lambda x,a : a #fit model function for scipy.optimize.curve_fit

def fit(f,xdata,ydata,yerr):
    '''
    make the fit on the data with its error and get the estimated value with its error on the constant model.
    <f> the model (here the uniform function defined above)
    <xdata> the bins data
    <ydata> the spectra data
    <yerr> the spectra error data
    '''
    a_esti, pcov = curve_fit(f, xdata, ydata, sigma=yerr)
    return a_esti[0], np.sqrt(pcov[0,0]) #return estimated value + error

####

filename = sys.argv[1] #the location of the output file from make_total_spectra.sh you want to plot
save     = sys.argv[2] #the directory the plot will be

tmp = np.load(fr'{filename}.npz')

bins = tmp['bins'] #keV
mask = (bins[:-1]>=2) & (bins[:-1]<=7)

fidcrops = list()

for i in tmp.files:
    if not i=='bins':
        category = i.split('_')[0]
        fidcrop  = get_fidcrop(i) #mm

        if not fidcrop in fidcrops : fidcrops.append(fidcrop) #mm


for fidcrop in fidcrops:
    print(f'{fidcrop:.1f}mm\t{get_fiducial_volume(fidcrop):.1f}%')

colors = [f'C{i}' for i in range(len(fidcrops))]
####

plt.style.use(r'./style.mplstyle')
supfig = plt.figure(figsize=(16,8))
grid = plt.GridSpec(1,2,hspace=.3)
fig = supfig.add_subfigure(grid[0])
totax = supfig.add_subplot(grid[1], title=f'Spectre Total (monosite + multisite)', xlabel='Énergie [keV]', ylabel='DRU [1/keV/kg/day]', yscale='log', xticks=np.arange(2,7.1,1), xlim=[2,7])#, ylim=[1E2,1E3])
del grid
grid = plt.GridSpec(2,2, top=.85, hspace=.4, left=.2)
ax = fig.add_subplot(grid[0,0], title='Fiduciel Monosite', ylabel='DRU [1/keV/kg/day]', yscale='log', xticks=np.arange(2,7.1,1), xlim=[2,7])#, ylim=[3E1,1E3])
bx = fig.add_subplot(grid[0,1], title='Fiduciel Multisite', yscale='log', sharey=ax, sharex=ax)
cx = fig.add_subplot(grid[1,0], xlabel='Énergie [keV]', ylabel='DRU [1/keV/kg/day]', title='Surface Monosite', yscale='log', sharex=ax, sharey=ax)
dx = fig.add_subplot(grid[1,1], title='Surface Multisite', xlabel='Énergie [keV]', yscale='log', sharex=bx, sharey=ax)
del grid


err_kw = {'fmt':'o', 'ms':.5, 'capsize':2, 'lw':1}

for color, fidcrop in zip(colors, fidcrops):
    print(f'{fidcrop:.1f}')
    #left side plots

    F1_mean, F1_err = fit(f, bins[1:], tmp[f'F1_{fidcrop:.1f}mm'][mask], tmp[f'ERRF1_{fidcrop:.1f}mm'][mask])
    F2_mean, F2_err = fit(f, bins[1:], tmp[f'F2_{fidcrop:.1f}mm'][mask], tmp[f'ERRF2_{fidcrop:.1f}mm'][mask])
    S1_mean, S1_err = fit(f, bins[1:], tmp[f'S1_{fidcrop:.1f}mm'][mask], tmp[f'ERRS1_{fidcrop:.1f}mm'][mask])
    S2_mean, S2_err = fit(f, bins[1:], tmp[f'S2_{fidcrop:.1f}mm'][mask], tmp[f'ERRS2_{fidcrop:.1f}mm'][mask])

    print(f'|\tF1 : {F1_mean:.0f} +/- {F1_err:.0f}')
    print(f'|\tF2 : {F2_mean:.0f} +/- {F2_err:.0f}')
    print(f'|\tS1 : {S1_mean:.0f} +/- {S1_err:.0f}')
    print(f'|\tS2 : {S2_mean:.0f} +/- {S2_err:.0f}')

    plot_spectra(ax,
                 tmp[f'F1_{fidcrop:.1f}mm'],
                 tmp[f'ERRF1_{fidcrop:.1f}mm'],
                 err_kw, color=color)

    plot_spectra(bx,
                 tmp[f'F2_{fidcrop:.1f}mm'],
                 tmp[f'ERRF2_{fidcrop:.1f}mm'],
                 err_kw, color=color)

    plot_spectra(cx,
                 tmp[f'S1_{fidcrop:.1f}mm'],
                 tmp[f'ERRS1_{fidcrop:.1f}mm'],
                 err_kw, color=color)

    plot_spectra(dx,
                 tmp[f'S2_{fidcrop:.1f}mm'],
                 tmp[f'ERRS2_{fidcrop:.1f}mm'],
                 err_kw, color=color)

    #right side plots
    F_mean, F_err = fit(f, bins[1:], tmp[f'F1_{fidcrop:.1f}mm'][mask] + tmp[f'F2_{fidcrop:.1f}mm'][mask],
                                     np.sqrt(tmp[f'ERRF1_{fidcrop:.1f}mm'][mask]**2 + tmp[f'ERRF2_{fidcrop:.1f}mm'][mask]**2))
    S_mean, S_err = fit(f, bins[1:], tmp[f'S1_{fidcrop:.1f}mm'][mask] + tmp[f'S2_{fidcrop:.1f}mm'][mask],
                                     np.sqrt(tmp[f'ERRS1_{fidcrop:.1f}mm'][mask]**2 + tmp[f'ERRS2_{fidcrop:.1f}mm'][mask]**2))

    print(f'|\tFiduTot : {F_mean:.0f} +/- {F_err:.0f}')
    print(f'|\tSurfTot : {S_mean:.0f} +/- {S_err:.0f}')
    print(f'|\tTotTot : {S_mean+F_mean:.0f} +/- {np.sqrt(S_err**2+F_err**2):.0f}')

    plot_spectra(totax,
                 tmp[f'F1_{fidcrop:.1f}mm'] + tmp[f'F2_{fidcrop:.1f}mm'],
                 np.sqrt(tmp[f'ERRF1_{fidcrop:.1f}mm']**2 + tmp[f'ERRF2_{fidcrop:.1f}mm']**2),
                 err_kw, color=color, ls='-')

    totax.hlines(F_mean, (bins[0]+bins[1])/2, (bins[-1]+bins[-2])/2, color=color, ls='--')
    totax.text((bins[-1]+bins[-2])*1.05/2, F_mean, fr'{F_mean:.0f}$\pm${F_err:.0f} DRU', ha='center', va='center', color=color)

    # plot_spectra(totax,
    #              tmp[f'S1_{fidcrop:.1f}mm'] + tmp[f'S2_{fidcrop:.1f}mm'],
    #              np.sqrt(tmp[f'ERRS2_{fidcrop:.1f}mm']**2 + tmp[f'ERRS2_{fidcrop:.1f}mm']**2),
    #              err_kw, color=color, ls='--')

###

totax.yaxis.set_ticks_position('both')
# totax.plot([],[],color='C0',ls='-', label='Fiduciel (Monosite + Multisite)') #; totax.plot([],[],color='C0',ls='--', label='Surface (Monosite + Multisite)')
# totax.legend()
fig.legend([f'{100*get_fiducial_volume(i):.0f}\\%' for i in fidcrops], loc=(.125,.9), title='Volume Fiduciel/Volume Détecteur', ncols=4)


###
# plt.show()
supfig.savefig(rf'{save}.pdf', dpi=300, transparent=True)
