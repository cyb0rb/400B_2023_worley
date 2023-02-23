# Homework 6
# C. Worley

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from Homeworks.Homework2.ReadFile import Read
from CenterOfMass2 import CenterOfMass

def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy: 'str'
            name of galaxy
        start: 'int'
            number of first snapshot to read in
        end: 'int'
            number of last snapshot to read in
        n: 'int'
            intervals over which COM is returned
    outputs: 
        x: 'text file'
            text file with COM position and velocity vectors at each snapshot
    """
    
    # compose the filename for output
    fileout = 'Homeworks/Homework6/Orbit_%s.txt' %galaxy
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    volDec = 2
    if galaxy == 'M33':
        volDec = 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids), 7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        prefix = '/home/cworley/Documents/400b/'    # folder where files exist
        ilbl = '000' + str(i)   # snap number
        ilbl = ilbl[-3:]        # snap number
        file = prefix + '%s_'%(galaxy) + 'VLowRes/' +'%s_'%(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(file, 2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        gal_COM_p = COM.COM_P(delta, volDec)
        gal_COM_v = COM.COM_V(gal_COM_p[0], gal_COM_p[1], gal_COM_p[2])

        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = COM.time.value / 1000, *tuple(gal_COM_p.value), *tuple(gal_COM_v.value)  # time in Gyr
        #print(i)
        
        # print snap_id to see the progress
        #print(snap_ids)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 

# OrbitCOM('M33', 0, 800, 5)
# OrbitCOM('MW', 0, 800, 5)
# OrbitCOM('M31', 0, 800, 5)

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
data_M33 = np.genfromtxt('Homeworks/Homework6/Orbit_M33.txt', dtype=None, names=True)
data_M31 = np.genfromtxt('Homeworks/Homework6/Orbit_M31.txt', dtype=None, names=True)
data_MW = np.genfromtxt('Homeworks/Homework6/Orbit_MW.txt', dtype=None, names=True)

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def vector_diff(data1, data2):
    ''' function computes the magnitude of the difference between the
        distance and velocity vectors
        
        PARAMETERS:
            data1:
                first data file
            data2:
                second data file
        OUTPUTS:
            p_mag:
                difference in position
            v_mag:
                difference in velocity
        '''
    
    # difference
    #print(np.asarray(data1['x']))
    #print(data2[0])
    diff = np.asarray(data1) - np.asarray(data2)
    # magnitude of difference in vectors
    p_mag = np.sqrt(diff['x']**2 + diff['y']**2 + diff['z']**2)
    v_mag = np.sqrt(diff['vx']**2 + diff['vy']**2 + diff['vz']**2)
    return p_mag, v_mag

# Determine the magnitude of the relative position and velocities 
# of MW and M31
relp_MW_M31, relv_MW_M31 = vector_diff(data_MW, data_M31)

# of M33 and M31
relp_M33_M31, relv_M33_M31 = vector_diff(data_M33, data_M31)

# Plot the Orbit of the galaxies 
#################################
fig, ax = plt.subplots()
ax.plot(data_M31['t'], relp_MW_M31)
fig.savefig('MW_M31_position.png')




# Plot the orbital velocities of the galaxies 
#################################

