
# Homework 4
# Center of Mass Position and Velocity
# Solutions: G.Besla, R. Li, H. Foote

# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from Homeworks.Homework2.ReadFile import Read

class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component (x)
        b : `float or np.ndarray of floats`
            second vector component (y)
        c : `float or np.ndarray of floats`
            third vector component (z)
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        # computing generic COM
        # total mass
        m_tot = np.sum(m)
        # xcomponent Center of mass
        a_com = np.sum(a * m) / m_tot
        # ycomponent Center of mass
        b_com = np.sum(b * m) / m_tot
        # zcomponent Center of mass
        c_com = np.sum(c * m) / m_tot
        
        # return the 3 components separately
        return a_com, b_com, c_com
       
    
    
    def COM_P(self, delta, volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        volDec: 'float'
            defines amount by which RMAX is decreased (RMAX/volDec)
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at reduced radius defined by volDec                                                          
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            # to check this                                                                                               
            # print ("CHANGE = ", change)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of volDec                                                                 
            r_max /= volDec
            # check this.                                                                                              
            # print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position   
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM position vector
        return np.round(p_COM, 2) * u.kpc
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        xV = np.abs(self.x*u.kpc - x_COM)
        yV = np.abs(self.y*u.kpc - y_COM)
        zV = np.abs(self.z*u.kpc - z_COM)
        rV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius
        indexV = np.where(rV < rv_max)
        
        # determine the velocity and mass of those particles within the max radius
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        v_COM = np.array([vx_COM, vy_COM, vz_COM])

        # return the COM vector
        # set the correct units using astropy
        # round all values               
        return np.round(v_COM, 2) * u.km / u.s                                                                         
     
    
# ANSWERING QUESTIONS
#######################
if __name__ == '__main__' : 

    # Create a Center of mass object for the MW, M31 and M33
    MW_COM = CenterOfMass("/home/cworley/Documents/400b/MW_000.txt", 2)
    M31_COM = CenterOfMass("/home/cworley/Documents/400b/M31_000.txt", 2)
    M33_COM = CenterOfMass("/home/cworley/Documents/400b/M33_VLowRes/M33_000.txt", 2)

    # below gives you an example of calling the class's functions
    # MW:   store the position and velocity COM
    MW_COM_p = MW_COM.COM_P(0.1, 2)
    MW_pmag = np.sqrt(np.sum(np.square(MW_COM_p)))
    print('Milky Way COM position components:', MW_COM_p)
    print('MW COM p magnitude:', np.round(MW_pmag,2))
    MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
    print('Milky Way COM velocity components:', MW_COM_v)
    MW_vmag = np.sqrt(np.sum(np.square(MW_COM_v)))
    print('MW COM v magnitude:', np.round(MW_vmag,2), '\n')

    # M31 position and velocity COM
    M31_COM_p = M31_COM.COM_P(0.1, 2)
    M31_pmag = np.sqrt(np.sum(np.square(M31_COM_p)))
    print('M31 COM position components:', M31_COM_p.value)
    print('M31 COM p magnitude:', np.round(M31_pmag,2))
    M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
    print('M31 COM velocity components:', M31_COM_v)
    M31_vmag = np.sqrt(np.sum(np.square(M31_COM_v)))
    print('M31 COM v magnitude:', np.round(M31_vmag,2), '\n')

    # M33 position and velocity COM
    M33_COM_p = M33_COM.COM_P(0.1, 2)
    M33_pmag = np.sqrt(np.sum(np.square(M33_COM_p)))
    print('M33 COM position components:', M33_COM_p)
    print('M33 COM p magnitude:', np.round(M33_pmag,2))
    M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
    print('M33 COM velocity components:', M33_COM_v)
    M33_vmag = np.sqrt(np.sum(np.square(M33_COM_v)))
    print('M33 COM v magnitude:', np.round(M33_vmag,2), '\n')

    # separation between MW and M31
    p_sep = MW_COM_p - M31_COM_p
    v_sep = MW_COM_v - M31_COM_v
    # magnitude of separation
    p_mag = np.sqrt(np.sum(np.square(p_sep)))
    v_mag = np.sqrt(np.sum(np.square(v_sep)))
    print('MW and M31 separation:', np.round(p_mag, 3))
    print('MW and M31 velocity difference:', np.round(v_mag, 3))

    # separation between M33 and M31
    p_sep2 = M33_COM_p - M31_COM_p
    v_sep2 = M33_COM_v - M31_COM_v
    # magnitude of separation
    p_mag2 = np.sqrt(np.sum(np.square(p_sep2)))
    v_mag2 = np.sqrt(np.sum(np.square(v_sep2)))
    print('M33 and M31 separation:', np.round(p_mag2, 3))
    print('M33 and M31 velocity difference:', np.round(v_mag2, 3), '\n')

    # q4
    '''
        The iterative process to determine the COM is important because
        particles are flung out to farther distances / higher velocities
        from gravitational interactions with the other galaxies, so
        iteratively calculating the COM will reduce the effect of these
        outliers on the COM position and velocity.
    '''

