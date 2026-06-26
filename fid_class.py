import numpy as np
from numba.experimental import jitclass
from numba import float64, int32

@jitclass
class event:
    r'''
    EVENT CLASS :
    =============
    /!\ calling this class requires the use of numba (file top lines) for
    compilation of the function and gain of time (~factor 10 !).
    /!\ the offset correction depends on the geometry considered. You need to
    update them manually in localR and zoffset_correction methods if you do not
    use RUN015 geometry for your simulations.

    * Inputs:
    | - data : np.ndarray of shape (4, Nhits). Its 4 rows are X (lab-frame of
    |reference hit X); Y (lab-frame of reference hit Y); Z (lab-frame of
    |reference hit Z); E (deposited energy per hit)
    |
    | - detbounds : the coordinates for separating the detectors (yet only in
    |the XY plane since there is no stack of detectors along the Z axis)
    |
    | - fidcrop : the distance between the cylinder surfaces of the crystal
    |volume and the fiducial volume, in meters.
    |
    +---------------------------------------------------------------------------
    |
    * Methods:
    | - filter : stores the data in the detector and fiducial category. This
    |method is already automatically applied when declaring an event, so it is
    |not supposed to be called.
    |
    | - get : to obtain the data calculed with filter when initializing the
    |event variable. Returns a numpy array in the shape (3,3,2), where the two
    |first coordinates (3,3) are the detector position index (X,Y) and the third
    |coordinate is for the fiducial and surface event, respectively.
    |
    | - localR : uses X and Y event lab coordinates to obtain a local coordinate
    |in the crystal for clarifying the selection in the event selection. Returns
    |an array of the same size as the input arrays in the event initialization.
    |
    | - localZ : allow to obtain the local-z coordinates from the
    |global-z coordinates by identifying the detector with the x and y global
    |coordinates. Does not return anything, it just changes the z vector.
    |
    | - getx, gety, getz, getr, gete : returns the associated position of energy
    |vectors, if you ever want to study space distribution for example.
    |
    +---------------------------------------------------------------------------
    '''
    radius:    float64
    height:    float64
    X:      float64[:]
    Y:      float64[:]
    Z:      float64[:]
    R:      float64[:]
    E:      float64[:]
    size:      int32
    energies: float64[:,:,:]
    #under the detector convention :
    #[[11,12,13],
    # [21,22,23],
    # [31,32,33]]

    def __init__(self, X, Y, Z, E, detbounds=(.025, .04), fidcrop=.5):
        self.radius = 15.358/1000 #m
        self.height = 10.640/1000 #m

        self.X = X
        self.Y = Y
        self.Z = self.localZ(Z, detbounds) #remove the z-offset of the global hit coordinates
        self.R = self.localR(X, Y, detbounds) #remove the x-y offset and create the local-r coordinate
        self.E = E

        self.size = (self.X).size

        #declare and fill the energies table
        self.energies = np.zeros((3,3,2), dtype=np.float64) #3*3 detectors * 2 subvolumes (Fiducial and Surface)
        self.filter(detbounds, fidcrop)

    def filter(self, detbounds, fidcrop):
        '''
        This function creates the energy tab in the 3*3*2 (for 3*3 detectors * 2
        volumes, fiducial and surface).
        * Inputs :
        | - detbounds : the coordinates for separating the detectors (yet only
        |in the XY plane since there is no stack of detectors along the Z axis
        |in RUN015 geometry).
        | - fidcrop : the distance between the cylinder surfaces of the crystal
        |volume and the fiducial volume, in meters.
        |
        * Outputs :
        | - the energy table of the event : storing the energy deposits in the
        |right place considering its position. The tab contains 18 places, for
        |each of the 9 detectors of RUN015 and for (0) fiducial volume and (1)
        |surface volume.
        '''
        (xbound, ybound) = detbounds #for the detector selection

        rbound = self.radius - fidcrop
        zbound = self.height/2 - fidcrop

        conZ = (np.abs(self.Z) <= zbound)
        conR = (self.R <= rbound)

        FIDCOND = (conZ & conR)

        #checking over each detector and filling the energies table
        for (i,j),_ in np.ndenumerate(self.energies[:,:,0]):
            ##x selection
            if   i==0: conx = (self.X < -xbound) #left column
            elif i==1: conx = (np.abs(self.X) < xbound) #center column
            else:      conx = (xbound < self.X) #right column
            ##y selection
            if   j==0: cony = (ybound < self.Y) #top row
            elif j==1: cony = (np.abs(self.Y) < ybound) #middle row
            else:      cony = (self.Y < -ybound) #bottom row

            DETCOND = (conx & cony)

            #build the energies depending on the volume and store them
            E_FIDU = np.sum( self.E[DETCOND & FIDCOND] )
            E_SURF = np.sum( self.E[DETCOND & (~FIDCOND)] )

            self.energies[i,j,0] = E_FIDU
            self.energies[i,j,1] = E_SURF

            #end of the loop

    def localR(self, X, Y, detbounds):
        '''
        Get the local radius coordinate in each detectors frame of reference.
        * Inputs :
        | - X,Y : the lab-frame of reference horizontal hit coordinates.
        | - detbounds : the coordinates for separating the detectors (yet only
        |in the XY plane since there is no stack of detectors along the Z axis
        |in RUN015 geometry).
        |
        * Outputs :
        | - the Oxy radius such that r^2=x^2+y^2 in the detector's frame of
        |reference.
        '''
        (xbound, ybound) = detbounds

        #the list of offsets on RUN015 geometry
        Xoff = [.055, 0, -.055]
        Yoff = [.068, 0, -.068]

        X = np.where(        X < -xbound, X + Xoff[0], X) #left column correction
        X = np.where(np.abs(X) < xbound,  X + Xoff[1], X) #middle column correction
        X = np.where(        X > xbound,  X + Xoff[2], X) #right column correction

        Y = np.where(        Y < -ybound, Y + Yoff[0], Y) #left column correction
        Y = np.where(np.abs(Y) < ybound,  Y + Yoff[1], Y) #middle column correction
        Y = np.where(        Y > ybound,  Y + Yoff[2], Y) #right column correction

        return np.sqrt(X**2 + Y**2)

    def localZ(self, Z, detbounds):
        '''
        Get the local vertical coordinate in each detectors frame of reference
        (since the detectors are not aligned in z).
        * Inputs :
        | - Z : the lab-frame of reference vertical hit coordinate.
        | - detbounds : the coordinates for separating the detectors (yet only
        |in the XY plane since there is no stack of detectors along the Z axis
        |in RUN015 geometry).
        |
        * Outputs :
        | - the vertical coordinate in the detector's frame of reference.
        '''
        (_, ybound) = detbounds

        Zoff = [1.15585, 1.16435, 1.17285]

        Z = np.where(self.Y < -ybound, Z - Zoff[0], Z)        #for lower  stage
        Z = np.where(np.abs(self.Y) < ybound, Z - Zoff[1], Z) #for middle stage
        Z = np.where(self.Y >  ybound, Z - Zoff[2], Z)        #for higher stage

        return Z

    def get(self):
        '''
        Get the energy table created in the "filter" method.
        * Outputs :
        | - self.energies : the energies of the table.
        '''
        return self.energies

    def getx(self): return self.X
    def gety(self): return self.Y
    def getr(self): return self.R
    def getz(self): return self.Z
    def gete(self): return self.E
