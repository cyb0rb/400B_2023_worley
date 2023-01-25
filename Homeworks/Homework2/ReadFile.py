# Homework 2

import numpy as np
import astropy.units as u

#filename = '/home/cworley/Documents/400b/MW_000.txt'

# Function reads text file and returns data as arrays
# inputs: filename (name of data file)
# returns: time, total number of particles, and full data array (data)
# data columns: [0] particle type, [1] mass (1e10Msun), 
#               [2-4] xyz positions from center of MW (kpc)
#               [5-7] xyz velocities from center of MW (km/s)
def Read(filename):
    file = open(filename, 'r')

    # reads first two lines
    line1 = file.readline()
    line2 = file.readline()

    # splits the label from value in first 2 lines
    label, value = line1.split()
    label2, value2 = line2.split()

    # applies units to values (if applicable)
    time = float(value) * u.Myr
    particles = float(value2)
    file.close()

    # rest of file from 4th row on is converted to array w/ header info
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, particles, data
