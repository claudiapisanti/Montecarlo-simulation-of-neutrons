"""
Utils for montecarlo program.
"""

# FUNCTIONS
import numpy as np
from random import random
import constants as c
import shutil
import sys

import time

from tqdm import tqdm


def random_rescale(val_max, val_min = 0):
    """
    random_rescale(val_max, val_min = 0)
    
    Return a random number between the interval of interest (val_min, val_max).

    Parameters:
    -----------
    val_max : float or int
        Upper boundary of the interval

    val_min: float or int, optional
        Lower boundary of the interval. Default is 0.

    Returns:
    -------
    val : float
        Random number between interval boundaries

    """
    val = random() * (val_max - val_min) + val_min
    return val

def source_position_est(face):
    """
    Define the position of the source in a rectangular surface. 
    
    Parameters:
    ----------
    face: int
        The face of the rectangle where the source must be placed

    Returns:
    -------
    pos: array
        Position of the source (array [x,y,z] ). 
        Two elements of the array are random between the position of the array, 
        while the third is selected in order to place the source on one of the face of the rectangle of interest.


    ---------
    The faces of the rectangle are numbered in the following way:

    In the reference system x,y,z as in the graph,
        1 --> Top
        2 --> Right
        3 --> Front
        4 --> Left
        5 --> Right
        6 --> Botton

    See figure for more clarification.

    
    In the reference:
        

               ________
              /       /|
             /   1   / | <-- 5 (backside)
          z /_______/  |                        
            |       | 2| y 
        4-> |   3   |  /
            |       | /  <-- 6 (bottom)
            |______ |/
           O         x
    """

    x_source = random_rescale(c.pos_max[0])
    y_source = random_rescale(c.pos_max[1])
    z_source = random_rescale(c.pos_max[2])
    
    # initialize position
    pos = [0,0,0]
    
    if (face == 1): pos = np.array([x_source, y_source, c.pos_max[2]])
    if (face == 2): pos = np.array([c.pos_max[0], y_source, z_source])
    if (face == 3): pos = np.array([x_source, c.pos_min[1], z_source])
    if (face == 4): pos = np.array([c.pos_min[0], y_source, z_source])
    if (face == 5): pos = np.array([x_source, c.pos_max[1], z_source])
    if (face == 6): pos = np.array([x_source, y_source, c.pos_min[2]])

    return pos



def get_cs(E, cs_table):   
    """
    Get cross section given Energy.


    Parameters:
    ----------
    E: float or int

    cs_table: DataFrame
        7 columns DataFrame table with all the cross section:

        E	total_cs	cs_el	cs_inel	cs_h_tot	cs_h_el	cs_h_inel	cs_c_tot	cs_c_el	cs_c_inel

        The function is specific for hydrocarbons, or molecules only constituted of carbon(C) and hydrogen(H).   
        
        Energy must be  in eV
        cross sections must be in barn (b)



    Returns:
    -------
    cs: array
        For the determined energy, cross section form cs_Table, converted in cm^2

        [ total cross section, 
        total elastic cross section, 
        total inelastic cross section, 
        total proton corss section, 
        proton elastic cross section, 
        proton inelastic cross section, 
        carbon total corss section, 
        carbon elastic cross section, 
        carbon inelastic cross section ]

    l: array
        Mean free path. The order of the array is the same of cs. In cm^2.

    p: array
        Array probability of having proton vs carbon scattering, and, 
        for the single element (proton and carbon), of having elastic vs inelastic scattering.

        [ prob of having carbon scattering, 
        prob of having proton scattering, 
        probability of having hydrogen elastic scattering, 
        probability of having hydrogen inelastic scattering,
        probability of having carbon elastic scattering, 
        probability of having carbon inelastic scattering ]

    """

    index = find_nearest(cs_table['E'], E)


    cs_tot = cs_table['cs_tot'][index]
    cs_el = cs_table['cs_el'][index]
    cs_inel = cs_table['cs_inel'][index]

    cs_h_tot = cs_table['cs_h_tot'][index]
    cs_h_el = cs_table['cs_h_el'][index]
    cs_h_inel = cs_table['cs_h_inel'][index]

    cs_c_tot = cs_table['cs_c_tot'][index]
    cs_c_el = cs_table['cs_c_el'][index]
    cs_c_inel = cs_table['cs_c_inel'][index]


    cs = [cs_tot, cs_el, cs_inel, 
          cs_h_tot, cs_h_el, cs_h_inel, 
          cs_c_tot, cs_c_el, cs_c_inel] # cs[0] = total, cs[1] = elastic, cs[2] = inelastic
    cs = np.dot(cs, 10e-24) # convert cs from barn to cm^2
    l = np.dot(cs, c.n) # lambda

    p_carbon = 9 * cs[6] / cs[0]
    p_proton = 10 * cs[3] / cs[0]

    # for hydrogen (H)
    p_h_elastic = cs[4] / cs[3] # probability of having elastic scattering
    p_h_inelastic = cs[5] / cs[3] # probability of having inelastic scattering with respect to total probability of having scattering 

    # for carbon(C)
    p_c_elastic = cs[7] / cs[6] # probability of having elastic scattering
    p_c_inelastic = cs[8] / cs[6] # probability of having inelastic scattering with respect to total probability of having scattering 

    p = [p_carbon, p_proton,
         p_h_elastic, p_h_inelastic,
         p_c_elastic, p_c_inelastic]

    return cs, l, p




def find_nearest(array,value):
    """
    This function find theindex of nearest value of a given one in an array.
    The function finds the upper value. 

    Parameters:
    ----------
    array: array
        list of values from which the nearest value of 'value' must be found

    value: int or float
        the value to be found

    Returns:
    -------
    idx: int
        the index of the value to be found.

    """

    idx = np.searchsorted(array, value, side="right") 

    assert idx < len(array), "Index value exceeds array length."

    return idx 
   

def face_func():
    """
    Randomly select a face of the rectangular extended source, given the cumulativefunction (c.face_prob_cum).

    Returns:
    -------
    out: int
        index of the face. See source_position_est() for face number legend.

    """

    r = random()
    index = find_nearest(c.face_prob_cum, r)
    return int(index +1)



def from_sph_coord_to_xyz(r,phi,theta,x0= 0.,y0=0.,z0=0.):
    """
    Conversion form spherical coordinates to cartesian coordinates

    Parameters:
    ----------
        r: float or int
            from 0 to infinity [0, inf ]
        phi: float or int
            in radiant, boundaries: [0, pi]

        theta: float or int
            in radiant, boundaries: [0, 2*pi]
        x0: float o int
            initial coordinates x. Deaflut = 0.
        y0: float or int
            initial coordinates y. Deaflut = 0.
        z0: float or int
            initial coordinates z. Deaflut = 0.

    Returns:
    -------
    pos: array  
        position in x,y,z coordinates


    """

    # NB: angle is in radiants
    
    x1 = x0 + r * np.cos(phi) * np.sin(theta)
    y1 = y0 + r * np.sin(phi) * np.sin(theta)
    z1 = z0 + r * np.cos(theta)
    pos = np.array([x1,y1,z1])
    return pos

def scattering_angle(E, A):
    """
    scattering for a given energy. 
    For the maximum loss of energy with respect to the different atom (H or C), 
    ((A - 1)/(A+1))**2*E has been used. (For neutrons)
    
    For more information on the formula see:
    Leo, W. R. (2012). Techniques for nuclear and particle physics experiments: a how-to approach. 
    Springer Science & Business Media.

    Parameters:
    ----------
    E: float
        energy of the neutron before the scattering

    A: int
        atomic number of the nucleus in which neutron is inpinging. 
        (eg. A = 12 for carbon, A = 1 for proton)
        

    Returns:
    -------
    theta_scat: float
        angle of scattering of the neutron. Properly the angle between neutron direction
        before and after scattering.

    E_new: float
        Energy of the neutron after the scattering.


    """

    E_new = random_rescale(E, ((A - 1)/(A+1))**2*E) # formula on LEO -- the energy of the scattered neutron is limited in the range between DeltaE and E0                 
    cos_theta_CM = ((E_new) / E * (A+1)**2 - A**2 - 1)/(2*A) # see LEO
    theta_scat = np.arccos((A * cos_theta_CM + 1)/ np.sqrt(A**2 + 1 + 2*A*cos_theta_CM)) # see LEO

    
    return theta_scat, E_new


def merge_tmp_tables(table_name : str, tmp_tables_list : list):
    """
    Function that merges together the temporary files created by the 
    different processes into a single final table
    Inputs:
        table_name : filename of the final table
        tmp_tables_list : list of temporary tables filenames
    """

    tmp_tables_list.sort()

    with open(table_name,'a') as table:
        for f in tmp_tables_list:
            with open(f,'r') as tmp:
                shutil.copyfileobj(tmp, table)



# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# single event
def event_func(i, cs_table): # [cs_table, k]
    """run events"""
    N_i = int(c.N / c.n_processes)

    for w in tqdm(range(N_i)): # run a certain number of events
        e = w + (i * N_i) +1 # number of the event
        
        

        # create temporary file
        step_name = "tmp_step%d.txt" % (i)
        event_name = "tmp_event%d.txt" % (i)


        with open(step_name, 'a+') as step, open(event_name, 'a+') as event:

            
            if(c.type_source == 'EST'): # get initial position
                face = face_func()
                pos_source = source_position_est(face)
            elif(c.type_source == 'SPH'): # get initial position
                face = 0
                phi_source = random_rescale(np.pi)
                theta_source = random_rescale(2*np.pi)
                pos_source = from_sph_coord_to_xyz(c.sph_radius,phi_source,theta_source,c.pos_max[0]/2.,c.pos_max[1]/2.,c.pos_max[2]/2.) # it is centered in the centre of the system
            elif(c.type_source == 'PUNT'): 
                face = 0
                pos_source = c.pos_source
            
            start2 = time.time()
            if(c.En_type == 'UNIF') :
                E = random_rescale(c.E_max, c.E_min)
                cs, l, p = get_cs(E,cs_table)
            elif(c.En_type == 'MONO'):
                E = c.E_mono
                cs, l, p = get_cs(c.E_mono, cs_table)

            end2 = time.time()
            # print('MONO --', end2 - start2)


            start = time.time()
            j = 0
            p_type = 0.
            x0 = pos_source[0]
            y0 = pos_source[1]
            z0 = pos_source[2]
            pos = np.array([0.,0.,0.]) # initialize
            
            step.write(f'{e}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\tsource\t{E}\n') # source position for each event

            DeltaE = 0
            theta_scat = 0
            while( (pos >= c.pos_min).all() & (pos <= c.pos_max).all()):

                # it is easier to work in spherical coordinates

                # initial direction
                if j == 0:
                    phi = random_rescale(2*np.pi) # phi is inside [0 ,2*pi ]
                    theta = random_rescale(np.pi)   # thetha is inside [] 0 e pi]
                else: 
                    rand = random_rescale(1, -1)
                    phi = phi + r * np.sin(theta_scat)* np.cos(r)  
                    theta = theta + r * np.sin(theta_scat) * np.sin(r) 
                

                p_interaction = random() # probability of interaction (cs_total). It is necessary for the measure of the lenght claculated form the particle 
                r = - l[0] * np.log(1-p_interaction) # by using the BEER LAMBERT law, i can calucate how long the particle is travelling
                
                # translate sphericla coordinates into cartesian coordinates
                pos = from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0) # position of interaction (x,y,z)
                
                # check if the particle is inside the scintillatro
                if(  (pos >= c.pos_min ).all()  & (pos <= c.pos_max).all() ):
                            
                    # check if neutron interact with carbon or proton
                    p_atom = random()
                    if (p_atom <= p[0]): # interact with carbon
                        # elastic or inelastic interaction?
                        p_type = random() 
                    
                        if(p_type <= p[4] ): # if elastic scattering with carbon
                            A = 12
                            
                            # calculate energy loss and so theta scattering of the neutron
                            theta_scat, E = scattering_angle(E, A)
                            if E < c.MeV1: break # energy threshold (for the energy resolution of my detector)

                            step.write(f'{e}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\t{E}\n')
                            
                            j = j+1 # next step
                            

                        
                        else: # if inelastic scattering with carbon
                            step.write(f'{e}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\t{E}\n')
                            event.write(f'{e}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
                            j = j+1 # step successivo

                            break

                    else: # interact with proton
                        # check type of interaction (elastic) or inelastic
                        p_type = random() 
                    
                        if(p_type <= p[2] ): # if elastic scattering with proton
                            A = 1
                            
                            # calculate energy loss and so theta scattering of the neutron
                            theta_scat, E = scattering_angle(E, A)
                            if E < c.MeV1: break # energy thershold (for the energy resolution of my detector)

                            step.write(f'{e}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\t{E}\n')
                            j = j+1 # next step

                        
                        else: # if inelastic scattering with proton
                            step.write(f'{e}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\t{E}\n')
                            event.write(f'{e}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
                            
                            break # I'm interested only in multiple scattering
                        
                        

                else: # if particle exit the scintillator
                    if ((pos!=pos_source).all()):
                        event.write(f'{e}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\n')
            
                x0,y0,z0 = pos



