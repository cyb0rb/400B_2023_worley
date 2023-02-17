import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.constants import G

# relevant functions / classes from previous
from Homeworks.Homework2.ReadFile import Read
from Homeworks.Homework4.CenterOfMass import CenterOfMass

class MassProfile:
    # class to calculate mass profile and circular velocites
    # of a galaxy using data snapshot and Hernquist profile
    def __init__(self, galaxy, snap):
        ''' Class to calculate mass profile for a 
            given galaxy at a certain snapshot in time
            
            PARAMETERS:
            galaxy: 'str'
                galaxy name
            snap: 'int'
                file number of snapshot
        '''
        # add string of filenumber to value '000'
        # and remove all but last 3 digits
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        self.filename = '%s_'%(galaxy) + ilbl + '.txt'

        # for testing purposes
        self.filename = '/home/cworley/Documents/400b/' + self.filename

        # saving galaxy name 
        self.gname = galaxy

        # read in data + apply units
        self.time, self.total, self.data = Read(self.filename)
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
        self.m = self.data['m']

    # mass enclosed function
    def MassEnclosed(self, p_type, radius):
        ''' Method computes the mass enclosed within a given
            radius of the center of mass of a specified galaxy
            and particle type

            PARAMETERS:
            p_type: 'int; 1, 2, or 3'
                particle type to use for mass enclosed
            radius: 'ndarray of floats'
                array of radii from COM in kpc

            OUTPUTS:
            mass_enc: 'ndarray of astropy quantity'
                array of enclosed masses for each radius (Msun)
        '''

        # determine COM position and magnitude
        gal_COM = CenterOfMass(self.filename, p_type)
        gal_COM_p = gal_COM.COM_P(0.1)
        gal_COM_pmag = np.sqrt(np.sum(np.square(gal_COM_p)))

        # store particles of type p_type
        p_index = np.where(self.data['type'] == p_type)
        x_type = self.x[p_index]
        y_type = self.y[p_index]
        z_type = self.z[p_index]
        m_type = self.m[p_index]

        # define which particles are enclosed within radius
        # at each array element by looping over array
        mass_enc = np.zeros_like(radius)
        for i in range(len(radius)):

            r_new = np.sqrt((x_type-gal_COM_p[0])**2 + 
                            (y_type-gal_COM_p[1])**2 + (z_type-gal_COM_p[2])**2)
            index = np.where(r_new < radius[i]*u.kpc)

            # grab particles within radius
            x_new = x_type[index]
            y_new = y_type[index]
            z_new = z_type[index]
            m_new = m_type[index]

            # sum over enclosed mass up to that radius
            mass_enc[i] = np.sum(m_new)
        
        return mass_enc * 1e10* u.Msun
    
    # total mass enclosed
    def MassEnclosedTotal(self, r):
        ''' method computes total mass enclosed at each radius of input array
        
            PARAMETERS:
            r: 'ndarray of floats'
                array of radii from COM in kpc
            OUTPUT:
            m_tot: 'ndarray of astropy.quantity'
                array of enclosed mass for each radius from all particles (Msun)
        '''

        # mass enclosed for halo and disk
        m_halo = self.MassEnclosed(1, r)
        m_disk = self.MassEnclosed(2, r)

        # total mass enclosed, includes bulge if not M33
        if self.gname != "M33":
            m_bulge = self.MassEnclosed(3, r)
            m_tot = m_halo + m_disk + m_bulge
        else:
            m_tot = m_halo + m_disk
        
        return m_tot
    
    # hernquist mass profile
    def HernquistMass(self, r, a, m_halo):
        ''' method computes mass enclosed within a given radius using Hernquist
            mass profile
            
            rho = (M*a / 2pi*r) / (r+a)**3
            M = (M_halo * r**2) / (r+a)**2

            PARAMETERS:
            r: 'float'
                distance from galactic center (kpc)
            a: 'float'
                scale factor of Hernquist profile (kpc)
            m_halo: 'astropy.quantity'
                total halo mass (Msun)
            OUTPUT:
            hern_mass: 'astropy.quantity'
                mass within given radius (Msun)
            '''
        hern_mass = (m_halo * r**2) / (a + r)**2
        return hern_mass
    
    # circular velocity
    def CircularVelocity(self, p_type, r):
        ''' method calculates circular speeds at each radius due to
            the enclosed mass of a given particle type, assuming
            spherical symmetry

            V(R) = sqrt(GM(<R) / R)

            PARAMETERS:
            p_type: 'int; 1, 2, 3'
                particle type to calcuate enclosed mass
            r: 'ndarray of floats'
                array of radii from COM in kpc
            OUTPUT:
            v_circ: 'ndarray of astropy.quantity'
                array of circular velocities (km/s)
        '''
        # gravitational constant in correct units
        Grav = G.to(u.kpc*u.km**2 / (u.Msun * u.s**2))

        # mass enclosed and velocity calculation
        m = self.MassEnclosed(p_type, r)
        v_circ = np.sqrt(Grav*m / (r*u.kpc))

        return np.round(v_circ, 2)


    # total circular velocity 
    def CircularVelocityTotal(self, r):
        ''' method calculates circular velocity at each radius of an
            array due to the effects of all particle types

            V(R) = sqrt(GM(<R) / R)
            
            PARAMETERS:
            r: 'ndarray of floats'
                array of radii from COM on kpc
            OUTPUT:
            v_circ: 'ndarray of astropy.quantity'
                array of circular velocities (km/s)
        '''
        # gravitational constant in correct units
        Grav = G.to(u.kpc*u.km**2 / (u.Msun * u.s**2))

        # total mass enclosed and velocity calculation
        m = self.MassEnclosedTotal(r)
        v_circ = np.sqrt(Grav*m / (r*u.kpc))

        return np.round(v_circ, 2)


    # hernquist circular speed
    def HernquistVCirc(self, r, a, m_halo):
        ''' method computes circular speed at a given radius using
            the mass enclosed given by the Hernquist mass profile

            V(R) = sqrt(GM(<R) / R)
        
            PARAMETERS:
            r: 'float'
                distance from galactic center (kpc)
            a: 'float'
                scale factor of Hernquist profile (kpc)
            m_halo: 'astropy.quantity'
                total halo mass (Msun)
            OUTPUT:
            v_circ: 'astropy.quantity'
                circular velocity at radius r (km/s)
        '''
        # gravitational constant in correct units
        Grav = G.to(u.kpc*u.km**2 / (u.Msun * u.s**2))

        # total mass enclosed and velocity calculation
        m = self.HernquistMass(r, a, m_halo)
        v_circ = np.sqrt(Grav*m / (r*u.kpc))
        
        return np.round(v_circ, 2)


# plotting
def plotmass(galaxy, r, a):
    ''' function plots mass profile and rotation curve
        of a galaxy from given data and using Hernquist
        profile

        PARAMETERS:
            galaxy: 'str'
                galaxy name
            r: 'ndarray of floats'
                array of radii from COM in kpc
            a: 'float'
                scale factor of Hernquist profile (kpc)
        '''
    gal = MassProfile(galaxy, 0)

    # galaxy mass components
    halo = gal.MassEnclosed(1, r)
    disk = gal.MassEnclosed(2, r)
    total = gal.MassEnclosedTotal(r)
    #print('%s halo:' %galaxy, f'{halo[-1]:.2e}')
    #print('%s total:' %galaxy, f'{total[-1]:.2e}')

    # hernquist mass - uses full halo instead of
    # the range provided by r like calculated above
    halo_index = np.where(gal.data['type'] == 1)[0]
    full_halo_mass = np.sum(gal.data['m'][halo_index] * 1e10 * u.M_sun)
    hern = gal.HernquistMass(r, a, full_halo_mass)
    #print('%s hernquist:' %galaxy, f'{hern[-1]:.2e}')

    # galaxy circular velocities
    v_halo = gal.CircularVelocity(1, r)
    v_disk = gal.CircularVelocity(2, r)
    v_total = gal.CircularVelocityTotal(r)
    v_hern = gal.HernquistVCirc(r, a, full_halo_mass)

    ### PLOTTING 
    fig, ax = plt.subplots(1,2, figsize=(12,6))

    # mass profile plot
    ax[0].plot(r, halo, ':r', label='halo')
    ax[0].plot(r, disk, ':b', label='disk')
    ax[0].plot(r, total, '-k', label='total')
    ax[0].plot(r, hern, '--m', label='Hernquist, a=%d' %a)
    ax[0].set(title='%s mass profile' %galaxy, xlabel='radius [kpc]', 
           ylabel='log(enclosed mass) [Msun]', yscale='log')

    # rotation curve plot
    ax[1].plot(r, v_halo, ':r', label='halo')
    ax[1].plot(r, v_disk, ':b', label='disk')
    ax[1].plot(r, v_total, '-k', label='total')
    ax[1].plot(r, v_hern, '--m', label='Hernquist, a=%d' %a)
    ax[1].set(title='%s rotation curve' %galaxy, xlabel='radius [kpc]', 
           ylabel='velocity [km/s]')

    # special case for M33 (no bulge)
    if galaxy != 'M33':
        bulge = gal.MassEnclosed(3, r)
        v_bulge = gal.CircularVelocity(3, r)
        ax[0].plot(r, bulge, ':g', label='bulge')
        ax[1].plot(r, v_bulge, ':g', label='bulge')

    ax[0].legend(loc='lower right')
    ax[1].legend(loc='lower right')
    fig.savefig('%s_massprofile.png' %galaxy)


if __name__ == '__main__':
    # values for functions
    r = np.arange(0.1, 30.1, 0.1)   # min max step

    # thanks avichal for the a values
    plotmass('MW', r, 60)
    plotmass('M31', r, 60)
    plotmass('M33', r, 25)