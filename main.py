import time ; start = time.time()
import sys ; this_filename = sys.argv[0].split('.')[0]
import json
import numpy as np
import matplotlib.pyplot as plt
import uproot

from make_analysis import analysis_file, NoEvent
from conversion import counts_to_dru, NoActivityVolume, NoActivityContaminant, NoSourceMass
from end_prgm import stop

#debug (display the informations while running)
debug = True

#files constants
path        = sys.argv[1] #the path for the simulations (here only one for the example)
simus_names = sys.argv[2] #the path for the names of the simulations
out_filename = sys.argv[3] #the path for the output file that will be plotted with plot_total_spectra.py

#analysis constants
fidcrops = [1*1E-3, 1.5*1E-3, 2*1E-3, 2.5*1E-3] #m
bins     = np.arange(1, 8, .5) #keV
ROI      = bins[-1] - bins[0] #keV
factor   = bins[1] - bins[0] #keV

#check if the name of the results file you will write on is the one you want
to_print = dict()
to_print['bins'] = bins

print(f'Data file output name : {out_filename}.npz\nAccept ? Y/n ', end='')
choice = input()
if not choice.upper()=='Y': raise KeyboardInterrupt('Filename not agreed')

#iterate over the fiducial cuts
for fidcrop in fidcrops:
    #the spectra vectors
    DRU_F1, DRU_F2, DRU_S1, DRU_S2 = (np.array([0 for i in range(bins.size - 1)], dtype='float64'),
                                      np.array([0 for i in range(bins.size - 1)], dtype='float64'),
                                      np.array([0 for i in range(bins.size - 1)], dtype='float64'),
                                      np.array([0 for i in range(bins.size - 1)], dtype='float64'))

    #the error on spectra vectors
    VAR_F1, VAR_F2, VAR_S1, VAR_S2 = (np.array([0 for i in range(bins.size - 1)], dtype='float64'),
                                      np.array([0 for i in range(bins.size - 1)], dtype='float64'),
                                      np.array([0 for i in range(bins.size - 1)], dtype='float64'),
                                      np.array([0 for i in range(bins.size - 1)], dtype='float64'))

    if debug : print(f'| Fidcrop : {fidcrop*1E3} mm')

    with open(simus_names, 'r') as f:

        data = json.load(f)

        #iterate over every volume of the geometry
        for key, item in data.items():
            volume = key[len('radiogenicGpsIn'):]

            if debug : print(f'|\t| Volume : {volume}')

            #iterate over every contaminant for the volume
            for contaminant in item:
                try :
                    if debug : print(f'|\t|\t| Contaminant : {contaminant}')

                    pathtofile = f'{path}/{key}_{contaminant}.root'

                    #get the histogram in counts
                    try :
                        if not ((volume=='LeadOuter_Top' or volume=='LeadOuter_Mid') and contaminant=='Th230-U238'): #used for possible focus over a contaminant or a volume..

                            try :
                                (counts_F1, counts_F2, counts_S1, counts_S2) = analysis_file(pathtofile, fidcrop, bins, debug=False)
                            except NoEvent : print(f'|\t|\t|\t| No event for : {volume}_{contaminant}.')

                            #convert the counts in DRUs
                            try :
                                DRU_F1 += counts_to_dru(counts_F1, volume, contaminant, ROI, factor)
                                DRU_F2 += counts_to_dru(counts_F2, volume, contaminant, ROI, factor)
                                DRU_S1 += counts_to_dru(counts_S1, volume, contaminant, ROI, factor)
                                DRU_S2 += counts_to_dru(counts_S2, volume, contaminant, ROI, factor)

                                VAR_F1 += counts_to_dru(np.sqrt(counts_F1), volume, contaminant, ROI, factor)**2
                                VAR_F2 += counts_to_dru(np.sqrt(counts_F2), volume, contaminant, ROI, factor)**2
                                VAR_S1 += counts_to_dru(np.sqrt(counts_S1), volume, contaminant, ROI, factor)**2
                                VAR_S2 += counts_to_dru(np.sqrt(counts_S2), volume, contaminant, ROI, factor)**2

                            #check for the custom exceptions of the modules i call and display the warning if this happens
                            except NoActivityVolume :      print(f'|\t|\t|\t| Could not find activity for : {volume}.')
                            except NoActivityContaminant : print(f'|\t|\t|\t| Could not find activity for : {volume}_{contaminant}.')
                            except NoSourceMass :          print(f'|\t|\t|\t| Could not find source mass for : {volume}.')

                    #check for an exception i met for some simulations that i could not open and display the warning if this happpens
                    except uproot.exceptions.KeyInFileError : print(f'|\t|\t|\t| Could not open simulation file for : {volume}_{contaminant}.')


                except : #in case any problem occurs during the whole combination, the results are saved in the file before stopping. Paricularly useful if you see this is taking too long but do not want to lose the results (via Ctrl+C in the terminal typically)
                    to_print.update({f'F1_{fidcrop*1E3}mm' : DRU_F1})
                    to_print.update({f'S1_{fidcrop*1E3}mm' : DRU_S1})
                    to_print.update({f'F2_{fidcrop*1E3}mm' : DRU_F2})
                    to_print.update({f'S2_{fidcrop*1E3}mm' : DRU_S2})

                    to_print.update({f'ERRF1_{fidcrop*1E3}mm' : np.sqrt(VAR_F1)})
                    to_print.update({f'ERRS1_{fidcrop*1E3}mm' : np.sqrt(VAR_S1)})
                    to_print.update({f'ERRF2_{fidcrop*1E3}mm' : np.sqrt(VAR_F2)})
                    to_print.update({f'ERRS2_{fidcrop*1E3}mm' : np.sqrt(VAR_S2)})

                    stop(to_print, start, out_filename)
                    raise KeyboardInterrupt(f'filed saved in {out_filename}.npz')


    #back to the fidcrop looping level, filling the dictionnary you will use to write the results file
    to_print.update({f'F1_{fidcrop*1E3}mm' : DRU_F1})
    to_print.update({f'S1_{fidcrop*1E3}mm' : DRU_S1})
    to_print.update({f'F2_{fidcrop*1E3}mm' : DRU_F2})
    to_print.update({f'S2_{fidcrop*1E3}mm' : DRU_S2})

    to_print.update({f'ERRF1_{fidcrop*1E3}mm' : np.sqrt(VAR_F1)})
    to_print.update({f'ERRS1_{fidcrop*1E3}mm' : np.sqrt(VAR_S1)})
    to_print.update({f'ERRF2_{fidcrop*1E3}mm' : np.sqrt(VAR_F2)})
    to_print.update({f'ERRS2_{fidcrop*1E3}mm' : np.sqrt(VAR_S2)})


##end the program
stop(to_print, start, out_filename)
