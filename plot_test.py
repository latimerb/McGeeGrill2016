import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb
import h5py

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


numPUDaff   = 1
numPelaff   = 1
numINmplus  = 1
numINmminus = 1
numIND      = 1
numFB       = 1
numSPN      = 1

config_file = "simulation_config.json"
v_file = "./output/cell_vars.h5"
s_file = "PUD_spikes.h5"

df = pd.read_csv("output/spikes.csv",delimiter=' ')
rast = df.values


i = h5py.File(v_file,'r')
mem_pot = i['report']['LUT']['data']


tsim = 24000

#print("Rin = {}".format((mem_pot[6000,0]-mem_pot[4000,0])/0.1))
#val = mem_pot[4000,0]+0.632*(mem_pot[6000,0]-mem_pot[4000,0])
#ind = find_nearest(mem_pot[:,0],val)
#print("tau = {}".format(np.arange(0,1000,1/10)[ind]-500))
plt.figure()
plt.subplot(6,1,1)
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,0])
plt.title('Pel')
plt.subplot(6,1,2)
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,1])
plt.title('PUD')
plt.subplot(6,1,3)
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,2])
plt.title('INm+')
plt.subplot(6,1,4)
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,3])
plt.title('INm-')
plt.subplot(6,1,5)
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,4])
plt.title('IND')
plt.subplot(6,1,6)
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,5])
plt.title('FB')

plt.figure()
plt.plot(np.arange(0,tsim,1/10),mem_pot[:,6])
plt.title('SPN')

labels_string = ['Pel','PUD','INm+','INm-','IND','FB','SPN']
rast = np.delete(rast,1,1)
fig,ax = plt.subplots()
for i in np.arange(0,len(labels_string)):
    rast[rast[:,1]==i,1] = labels_string[i]
ax.plot(rast[:,0],rast[:,1],'b.')
plt.show()

# Plot Bladder afferent, EUS afferent, PAG afferent, IND, and Hypo on one figure
#plt.figure()
#Bladspkt = rast[np.in1d(rast[:,1],Blad_gids),:]
#plt.plot(Bladspkt[:,0],Bladspkt[:,1],'b.',label='Bladder afferent')
#
#EUSspkt = rast[np.in1d(rast[:,1],EUS_gids),:]
#plt.plot(EUSspkt[:,0],EUSspkt[:,1],'r.',label='EUS afferent')
#
#PAGspkt = rast[np.in1d(rast[:,1],PAG_gids),:]
#plt.plot(PAGspkt[:,0],PAGspkt[:,1],'m.',label='PAG afferent')
#
#plt.xlabel('Time (t) [ms]')
#plt.title('Afferent Firing Times')
#plt.legend()
#
## Plot INM+, INM-, IND, FB --------------------------------
#plt.figure()
#INMpspkt = rast[np.in1d(rast[:,1],INmplus_gids),:]
#plt.plot(INMpspkt[:,0],INMpspkt[:,1]-40,'b.',label='INm+')
#
#INMmspkt = rast[np.in1d(rast[:,1],INmminus_gids),:]
#plt.plot(INMmspkt[:,0],INMmspkt[:,1]-40,'r.',label='INM-')
#
#INDspkt = rast[np.in1d(rast[:,1],IND_gids),:]
#plt.plot(INDspkt[:,0],INDspkt[:,1]-30,'g.',label='IND')
#
#FBspkt = rast[np.in1d(rast[:,1],FB_gids),:]
#plt.plot(FBspkt[:,0],FBspkt[:,1]-50,'k.',label='FB')
#
#plt.xlabel('Time (t) [ms]')
#plt.title('Interneuron Firing Times')
#plt.legend()
#
## Plot SPN(PGN), Hypo, MPG, IMG ---------------------------
#plt.figure()
#PGNspkt = rast[np.in1d(rast[:,1],PGN_gids),:]
#plt.plot(PGNspkt[:,0],PGNspkt[:,1]-50,'b.',label='PGN')
#
#IMGspkt = rast[np.in1d(rast[:,1],IMG_gids),:]
#plt.plot(IMGspkt[:,0],IMGspkt[:,1]-80,'r.',label='IMG')
#
#MPGspkt = rast[np.in1d(rast[:,1],MPG_gids),:]
#plt.plot(MPGspkt[:,0],MPGspkt[:,1]-70,'c.',label='MPG')
#
#Hypospkt = rast[np.in1d(rast[:,1],Hypo_gids),:]
#plt.plot(Hypospkt[:,0],Hypospkt[:,1]-40,'m.',label='Hypo')
#
#plt.xlabel('Time (t) [ms]')
#plt.title('Ganglion/Preganglionic Firing Times')
#plt.legend()
#
## Motor neurons -------------------------------------------
#plt.figure()
#EUSmnspkt = rast[np.in1d(rast[:,1],EUSmn_gids),:]
#plt.plot(EUSmnspkt[:,0],EUSmnspkt[:,1],'r.',label='EUS MN')
#
#Bladmnspkt = rast[np.in1d(rast[:,1],Bladmn_gids),:]
#plt.plot(Bladmnspkt[:,0],Bladmnspkt[:,1],'g.',label='Bladder MN')
#
#plt.xlabel('Time (t) [ms]')
#plt.title('Motor Neuron Firing Times')
#plt.legend()
#
#plt.show()
