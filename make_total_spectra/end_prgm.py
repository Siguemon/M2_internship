import numpy as np
import time

def stop(to_print, start, name):
    '''
    to_print is a dict for which the key is the name of the vector, and the associated item the vector of course. Can be used with totot.files for toto = np.load(...)
    <to_print> dictionnary so that you can save in a .npz file vectors linked to keys. This way you can use toto.files in another file to see the different keys and arrays you have in the file.
    <start> the time at which you main.py program started, so that the runtime can be shown.
    <name> the name of the file you want to write your total spectra data on.
    '''

    ### save everything to txt file (to avoid wasting time for a plot)
    np.savez(fr'{name}.npz', **to_print)
    print(fr'Saved file as {name}.npz')

    ### end the prgm
    runtime = int(time.time() - start)
    print(f'Runtime : {runtime//60}m {runtime%60}s')
