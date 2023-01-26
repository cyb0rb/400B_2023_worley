# Homework 2

import numpy as np
import astropy.units as u
from ReadFile import Read

filename = '/home/cworley/Documents/400b/MW_000.txt'

def ParticleInfo(filename, p_type, p_num):
    """ for a given file, extracts various properties about any particle of a given type
        Inputs:
            filename: 'string'
                path to particle data file
            p_type: 'integer'
                type of particle (1 = halo, 2 = disk, 3 = bulge)
            p_num: 'integer'
                index of nth particle of type p_type
        Returns:
            dist: 'astropy quantity'
                magnitude of particle distance from the galactic center (kpc)
            vel: 'astropy quantity'
                magnitude of particle velocity (km/s)
            m: 'astropy quantity'
                mass of particle (Msun)
    """

    #units in file
    mass = 10**10 * u.Msun
    v = u.km / u.second

    time, particles, data = Read(filename)     # reads in file
    type_i = np.where(data['type'] == p_type)  # index particles of type = p_type
    of_type = data[0:][type_i]                 # array of just particles of type = p_type

    p = of_type[p_num]                              # array for particle at index p_num
    m = p[1] * mass                                 # particle mass in 1e10 Msun
    x, y, z = p[2]*u.kpc, p[3]*u.kpc, p[4]*u.kpc    # particle position (x, y, z) in kpc
    vx, vy, vz = p[5]*v, p[6]*v, p[7]*v             # particle velocity (vx, vy, vz) in km/s

    # calculation of distance / velocity magnitudes
    dist = np.around(np.sqrt(x**2 + y**2 + z**2), 3)
    vel = np.around(np.sqrt(vx**2 + vy**2 + vz**2), 3)

    return dist, vel, m

dist, vel, m = ParticleInfo(filename, 2, 99)
print(dist)
print(np.around(dist.to(u.lyr), 3))
print(vel)
print(m)



