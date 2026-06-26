import numpy as np
import fid_class as C
import uproot

def analysis_file(pathtofile, fidcrop, bins, debug=False):
    '''
    give the counts for the four different categories for a given simulation.
    * Inputs :
    <pathtofile> the path to the simulation file you want to analyze
    <fidcrop> the distance between the cylinder surfaces of the crystal volume
    and the fiducial volume, in meters.
    <bins> the bounds of the bins as for standard use of numpy.histogram()
    <debug> bool to display the process if a problem occurs
    * Outputs :
    the counts for monosite or multiste events in the fiducial or surface
    volume (F1 stands for monosite fiducial, S2 for multisite surface).
    '''
    ## open the file and set the data
    tree = uproot.open(pathtofile+':ESO')
    if debug : print(f'Opening {pathtofile}')


    (X,Y,Z,E) = (tree[key].array(library='np') for key in ('HitX','HitY','HitZ','Edep'))

    if not X.size :
        raise NoEvent('no event registered for this simulation.')
    else :
        theta = -30 * np.pi/180 #rotation for the detector selection
        (X,Y) = (np.cos(theta)*X - np.sin(theta)*Y,
                 np.sin(theta)*X + np.cos(theta)*Y)
        del theta

        E = 1E3 * E #conversion MeV to keV

        ## get what we want from the data (the sorted energies)
        tmp = (C.event(x,y,z,e,fidcrop=fidcrop) for x,y,z,e in zip(X,Y,Z,E))
        output = np.array([i.get() for i in tmp])
        del tmp

        (i,j) = 1,0 #ids for the FID147
        multisite = np.logical_and(output[:,i,j,0], output[:,i,j,1])
        Mono_Fidu = np.where(~multisite, output[:,i,j,0], 0)
        Mono_Surf = np.where(~multisite, output[:,i,j,1], 0)
        Duo_Fidu  = np.where( multisite, output[:,i,j,0], 0)
        Duo_Surf  = np.where( multisite, output[:,i,j,1], 0)

        del output

        counts_F1, _ = np.histogram(Mono_Fidu, bins=bins)
        counts_F2, _ = np.histogram(Duo_Fidu,  bins=bins)
        counts_S1, _ = np.histogram(Mono_Surf, bins=bins)
        counts_S2, _ = np.histogram(Duo_Surf,  bins=bins)

        if debug : print(f'End of analysis for {pathtofile}\n')
        return (counts_F1, counts_F2, counts_S1, counts_S2)

### custom exceptions
class NoEvent(Exception): pass
