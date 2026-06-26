import numpy as np
import json
import csv

#every path used and that you need to change when using it
pathtoactivities = r'./data/activities.json'
pathtomapping    = r'./data/material_mapping.json'
pathtomasses     = r'./data/masses_RUN015.txt'


def counts_to_dru(counts, volume, contaminant, ROI, factor):
    '''
    <counts> the vector that contains the histogram counts that you may convert
    <volume> the name of the volume you have contaminated for the simulation
    <contaminant> the name of the contaminant you have used for the simulation
    <ROI> the Region of Interest, in keV
    <factor> is the amount of keVs inside a bin : this allows to have proper DRU estimation instead of 1 event/day/kg/0.1 keV for 0.1 keV bins for example
    '''
    nprimaries    = get_nprimaries() #nounit
    activity      = get_activity(volume, contaminant) #Bq/kg
    source_mass   = get_source_mass(volume) #kg
    detector_mass = get_detector_mass() #kg

    return ((counts/nprimaries) *
            activity *
            (source_mass/detector_mass) *
            24*3600 *
            1/ROI *
            1/factor) #1/keV/day/kg

def get_nprimaries():
    '''
    returns the number of primaries (yet it was the same for every simulation i used)
    '''
    return 2E6 #yet considered constant for every simulation (to check !!!!!)

def get_activity(volume, contaminant):
    '''
    returns the activity of the given volume for its given contamination (if its measurement exists !)
    <volume> the name of the volume you have contaminated for the simulation
    <contaminant> the name of the contaminant you have used for the simulation
    '''


    with open(pathtomapping, 'r') as f_mapping:
        with open(pathtoactivities, 'r') as f_activities:
            data_mapping = json.load(f_mapping)
            data_activities = json.load(f_activities)

            if not volume in data_mapping.keys(): raise NoActivityVolume(f'''volume '{volume}' not found in the mapping file''')
            else :
                goodname = data_mapping[volume]
                if not contaminant in data_activities[goodname].keys(): raise NoActivityContaminant(f'''contaminant '{contaminant}' not found in the activities file''')
                else: return data_activities[goodname][contaminant] * 1E-3 #Bq/kg

def get_detector_mass():
    '''
    return the mass of the detectors. Since this is the same for every simulation, i just return its value.
    '''
    return 42E-3 #kg, constant for every simulation since detector mass is always fixed

def get_source_mass(volume):
    '''
    returns the mass of the volume contaminated. This requires the creation of a
    dictionnary of masses created by the function "make_mass_dict" given below
    and called at the end of this program.
    <volume> the name of the volume you have contaminated for the simulation
    '''
    if not volume in source_masses.keys(): raise NoSourceMass('volume not found in the masses file')
    else: return source_masses[volume] #kg, sample function before obtaining the masses of the sources

#for the masses
def make_mass_dict():
    '''
    makes the dictionnary of masses so that the "get_source_mass" can be used
    properly. This function is automatically called at the end of this program.
    '''
    with open(pathtomasses, newline='') as f:
        rep_dict = dict()
        next(f) #skip the title line
        tmp = csv.reader(f, delimiter=' ', quotechar='|')
        for i,row in enumerate(tmp):
            # if i<301 : rep_dict[row[2]] = float(row[1]) #kg
            if 302<i<355 : rep_dict[row[1]] = float(row[2]) #kg

        return rep_dict

source_masses = make_mass_dict()

#define custom exception for the different cases
class NoActivityVolume(Exception) : pass
class NoActivityContaminant(Exception) : pass
class NoSourceMass(Exception) : pass


################################################################################
# when you run this program, you have the list of the different volumes that have registered masses to see if everything is alright
if __name__ == '__main__':
    print('Running program directly. Function test displayed.')

    print('source masses dict keys')
    for key in source_masses.keys():
        print(key)
