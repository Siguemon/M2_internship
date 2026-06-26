import time; start = time.time()
import sys
import numpy as np
import matplotlib.pyplot as plt
import uproot
import fid_class as C

try : #get the parameters from the arguments ran with the python call
    pathtofile = str(sys.argv[1])
    file = pathtofile.split(r'/')[-1]
    path = pathtofile[:-1-len(file)]
    filename = file.split('.')[0]
    contaminant, isotope = tuple(filename.split('_'))
    print(filename)

except :
    raise ValueError(f'Usage : $ py {sys.argv[0]} <path/to/file> <path/to/save=current> <upper_bin=3E3>')


try :    savepath = str(sys.argv[2])
except : savepath = '.' #current directory if no savepath specified

try : ncrops = int(sys.argv[3]) #get the number of points for the plot, default as 2
except : ncrops = 2

#some debug display
print(f'ncrops={ncrops}')
print(f'Running {sys.argv[0]} for {filename} (', end='')

################################################################################
#the dimensions of the crystal to plot the proportions as a function of the fiducial volume, and not the fidcrop parameter
crystal_radius = 15.358/1000 #m
crystal_height = 10.640/1000 #m

vol_prop = lambda crop : (crystal_radius - crop)**2/crystal_radius**2 * (crystal_height - crop)/crystal_height #proportion, then percentage unit !

################################################################################
#get the simulation file and make the data corrections

tree = uproot.open(pathtofile+':ESO')
(X,Y,Z,E) = (tree[key].array(library='np') for key in ['HitX','HitY','HitZ','Edep'])

theta = -30*np.pi/180
X, Y = (np.cos(theta)*X - np.sin(theta)*Y,
        np.sin(theta)*X + np.cos(theta)*Y)
del theta

E = 1E3 * E
NN = X.size
print(f'{NN} evts)')
################################################################################

plt.style.use(r'./style.mplstyle')
fig = plt.figure(figsize=(12,5))

#grid = plt.GridSpec(1, 2, hspace=.3, left=.15, right=.9, top=.85) #vertical disposition
grid = plt.GridSpec(1, 2, wspace=.3, bottom=.17, top=.8) #horizontal disposition

ax = fig.add_subplot(grid[0], title='Monosite', ylabel='\\#Events ratio (\\%)', xlabel='Fiducial Volume/Total Volume (\\%)')
bx = fig.add_subplot(grid[1], title='Multisite', xlabel='Volume Fiduciel/Volume Cristal (\\%)')#, ylabel='\\#Events ratio (\\%)')

del grid

#the data that we will plot
crops = np.linspace(0, crystal_height, ncrops) #m
(yF1, yF2, yS1, yS2) = (np.array([]), np.array([]), np.array([]), np.array([]))

#loop over every x and get the proportion of event for fiducial/surface monosite/multisite events categories
for i,crop in enumerate(crops):
    print(f'{i/crops.size*100:.2f}%', end='\r') #some debug information

    #get the energies table (see documention of fid_class for more information)
    tmp = (C.event(x,y,z,e, fidcrop=crop) for x,y,z,e in zip(X,Y,Z,E))
    output = np.array([i.get() for i in tmp])
    del tmp

    #total number of events
    Nevts = np.sum(np.where(output, 1, 0))

    #get the arrays of energies depending on the category
    multisite = np.logical_and(output[:,:,:,0], output[:,:,:,1])
    Mono_Fidu = np.where(~multisite, output[:,:,:,0], 0)
    Mono_Surf = np.where(~multisite, output[:,:,:,1], 0)
    Duo_Fidu  = np.where( multisite, output[:,:,:,0], 0)
    Duo_Surf  = np.where( multisite, output[:,:,:,1], 0)

    #get the number of events using the arrays given above
    NMono_Fidu = np.sum(np.where(Mono_Fidu, 1, 0))
    NMono_Surf = np.sum(np.where(Mono_Surf, 1, 0))
    NDuo_Fidu = np.sum(np.where(Duo_Fidu, 1, 0))
    NDuo_Surf = np.sum(np.where(Duo_Surf, 1, 0))

    #save the proportion of events in the data to plot
    yF1 = np.append(yF1, NMono_Fidu/Nevts*100)
    yS1 = np.append(yS1, NMono_Surf/Nevts*100)
    yF2 = np.append(yF2, NDuo_Fidu/Nevts*100)
    yS2 = np.append(yS2, NDuo_Surf/Nevts*100)
    #end of loop

#plot every curve proportion
ax.plot(vol_prop(crops), yF1,color='C0')
ax.plot(vol_prop(crops), yS1,color='C3')
bx.plot(vol_prop(crops), yF2,color='C0')
bx.plot(vol_prop(crops), yS2,color='C3',ls='--')

#some figure parameters
ax.set_xticks(np.arange(0,1.1,.2), [f'{i*100:.0f}\\%' for i in np.arange(0,1.1,.2)])
bx.set_xticks(np.arange(0,1.1,.2), [f'{i*100:.0f}\\%' for i in np.arange(0,1.1,.2)])

fig.legend(['Fiduciel','Surface'], loc=(.4,.85), fontsize=12, ncols=2)
fig.suptitle(f'''Proportion d'Events vs volume fiduciel : {filename}''')

#show or save the figure
# plt.show()
fig.savefig(fr'{savepath}/{filename}.pdf', dpi=300)
print(fr'Picture saved as {savepath}/{filename}.png')
###
print(f'Prgm runtime : {time.time() - start:.2f} s')
