# Homework 2

import numpy as np
import astropy.units as u

#filename = '/home/cworley/Documents/400b/MW_000.txt'

def Read(filename):
    """ Function reads text file and returns galaxy particle data
        Inputs:
            filename: 'string'
                path to particle data file
        Returns:
            time: 'astropy quantity'
                labeled time in Myr from file
            particles: 'integer'
                total number of particles in file
            data: 'numpy array'
                particle data from file with column labels
                    'type' particle type (1=halo, 2=disk, 3=bulge)
                    'm' particle mass (1e10 Msun)
                    'x' 'y' 'z' xyz positions from galactic center (kpc)
                    'vx' 'vy' 'vz' xyz velocities (km/s)
    """

    # opens file
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
