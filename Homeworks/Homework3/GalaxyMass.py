import numpy as np
import astropy.units as u
from Homeworks.Homework2.ReadFile import Read

# file paths
milkyway = '/home/cworley/Documents/400b/MW_000.txt'
m31 = '/home/cworley/Documents/400b/M31_000.txt'
m33 = '/home/cworley/Documents/400b/M33_000.txt'

def ComponentMass(filename, p_type):
    ''' function returns the total mass of the desired galaxy component
        Inputs:
            filename: 'string'
                path to particle data file
            p_type: 'integer'
                type of particle (1 = halo, 2 = disk, 3 = bulge)
        Outputs:
            M: 'astropy quantity'
                total mass of all the particles of given type (1e12 Msun)
    '''
    time, particles, data = Read(filename)     # read file

    type_i = np.where(data['type'] == p_type)  # index particles of type = p_type
    of_type = data[0:][type_i]                 # array of just particles of type = p_type

    M = np.sum(of_type['m']) / 100 * u.Msun  # sum of all particle masses
    return np.round(M, 3)


# to obtain answers for table
'''
dark_mw = ComponentMass(milkyway, 1)
disk_mw = ComponentMass(milkyway, 2)
bulge_mw = ComponentMass(milkyway, 3)
mass_mw = dark_mw + disk_mw + bulge_mw
f_mw = (disk_mw + bulge_mw) / mass_mw
print('milky way mass:', f'{mass_mw:.3f}')
print('milky way baryon fraction:', f'{f_mw:.3f}')

dark_m31 = ComponentMass(m31, 1)
disk_m31 = ComponentMass(m31, 2)
bulge_m31 = ComponentMass(m31, 3)
mass_m31 = dark_m31 + disk_m31 + bulge_m31
f_m31 = (disk_m31 + bulge_m31) / mass_m31
print('m31 mass:', f'{mass_m31:.3f}')
print('m31 baryon fraction:', f'{f_m31:.3f}')

dark_m33 = ComponentMass(m33, 1)
disk_m33 = ComponentMass(m33, 2)
bulge_m33 = ComponentMass(m33, 3)
mass_m33 =  dark_m33 + disk_m33 + bulge_m33
f_m33 = (disk_m33 + bulge_m33) / mass_m33
print('m33 mass:', f'{mass_m33:.3f}')
print('m33 baryon fraction:', f'{f_m33:.3f}')

dark_lg = dark_mw + dark_m31 + dark_m33
disk_lg = disk_mw + disk_m31 + disk_m33
bulge_lg = bulge_mw + bulge_m31 + bulge_m33
mass_lg = mass_mw + mass_m31 + mass_m33
f_lg = (disk_lg + bulge_lg) / mass_lg
print('local group mass:', f'{mass_lg:.3f}')
print('local group halo:', f'{dark_lg:.3f}')
print('local group disk:', f'{disk_lg:.3f}')
print('local group bulge:', f'{bulge_lg:.3f}')
print('local group baryon fraction:', f'{f_lg:.3f}')
'''

