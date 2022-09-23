"""
Utils for montecarlo program.
"""

# FUNCTIONS
import numpy as np
from random import random, seed
import shutil
import csv
from tqdm import tqdm


def random_rescale(rand, val_max, val_min = 0):
    """
    random_rescale(rand, val_max, val_min = 0)
    
    Return a random number between the interval of interest (val_min, val_max).

    Parameters:
    -----------
    rand: float
        Value between 0 and 1

    val_max : float or int
        Upper boundary of the interval

    val_min: float or int, optional
        Lower boundary of the interval. Default is 0.
    

    Returns:
    -------
    val : float
        Random number between interval boundaries

    """
    val = rand * (val_max - val_min) + val_min
    return val

def source_position_est(r1,r2,r3, face, pos_max, pos_min):
    """
    source_position_est(r1,r2,r3, face, pos_max, pos_min)

    Define the position of the source in a rectangular surface. 
    
    Parameters:
    ----------
    r1: float 
        random number between 0 and 1

    r2: float 
        random number between 0 and 1

    r3: float 
        random number between 0 and 1


    face: int
        The face of the rectangle where the source must be placed
    
    pos_max: array
        Dimension of the plastic scintillator

    pos_min: array
        (0,0,0) array, minimum coordinates of ths scintillator

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

    
    x_source = random_rescale(r1, pos_max[0])
    y_source = random_rescale(r2, pos_max[1])
    z_source = random_rescale(r3, pos_max[2])
    
    # initialize position
    pos = [0,0,0]
    
    if (face == 1): pos = np.array([x_source, y_source, pos_max[2]])
    if (face == 2): pos = np.array([pos_max[0], y_source, z_source])
    if (face == 3): pos = np.array([x_source, pos_min[1], z_source])
    if (face == 4): pos = np.array([pos_min[0], y_source, z_source])
    if (face == 5): pos = np.array([x_source, pos_max[1], z_source])
    if (face == 6): pos = np.array([x_source, y_source, pos_min[2]])

    return pos



def get_cs(E, cs_table, n):   
    """
    get_cs(E, cs_table, n)

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

    n: float
        Molecular density of the material 



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
    l = cs[0] * n  # lambda

    p_carbon = 9 * cs[6] / cs[0]
    p_proton = 10 * cs[3] / cs[0]

    # for hydrogen (H)
    p_h_elastic = cs[4] / cs[3] # probability of having elastic scattering
    p_h_inelastic = cs[5] / cs[3] # probability of having inelastic scattering with respect to total probability of having scattering 

    # for carbon(C)
    p_c_elastic = cs[7] / cs[6] # probability of having elastic scattering
    p_c_inelastic = cs[8] / cs[6] # probability of having inelastic scattering with respect to total probability of having scattering 

    p = [p_carbon*p_c_elastic, 
         p_carbon*p_c_inelastic,
         p_proton*p_h_elastic,
         p_proton*p_h_inelastic]

    # get cumulative
    p_cumulative = np.zeros(len(p))

    for i in range(len(p)):
        p_cumulative[i] = sum(p[:(i+1)]) / sum(p)


    return cs, l, p_cumulative




def find_nearest(array,value):
    """
    find_nearest(array,value)

    This function find the index of nearest value of a given one in an array.
    The function finds the upper value. 

    Parameters
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
   

def face_func(r, face_prob_cum):
    """
    face_func(r, face_prob_cum)

    Randomly select a face of the rectangular extended source, given the cumulative function (face_prob_cum).

    Parameters:
    ----------
    r: float
        random number between 0 and 1

    face_prob_cum: array
        a 6 dimension array with the cumulative distribution probability that the source is fount in one of the faces of the scintillator.
        It is calculated in file make_constants.py 

    Returns:
    -------
    out: int
        index of the face. See source_position_est() for face number legend.

    """

    index = find_nearest(face_prob_cum, r)
    return int(index +1)



def from_sph_coord_to_xyz(r,phi,theta,x0= 0.,y0=0.,z0=0.):
    """
    from_sph_coord_to_xyz(r,phi,theta,x0= 0.,y0=0.,z0=0.)

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

def scattering_angle(E, A, rand):
    """
    scattering_angle(E, A, rand)
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

    rand: float
        random number between 0 and 1.


    Returns:
    -------
    theta_scat: float
        angle of scattering of the neutron. Properly the angle between neutron direction
        before and after scattering.

    E_new: float
        Energy of the neutron after the scattering.


    """
    E_new = random_rescale(rand, E, ((A - 1)/(A+1))**2*E) # formula on LEO -- the energy of the scattered neutron is limited in the range between DeltaE and E0                 
    cos_theta_CM = ((E_new) / E * (A+1)**2 - A**2 - 1)/(2*A) # see LEO
    theta_scat = np.arccos((A * cos_theta_CM + 1)/ np.sqrt(A**2 + 1 + 2*A*cos_theta_CM)) # see LEO

    
    return theta_scat, E_new


def get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3,r4):
    """
    get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3,r4)
    
    Get source position given configuration informations (the input parameters of the source)

    Parameters:
    ----------
    type_source: string
        Source geometry (PUNT, EST or SPH).    

    source_params: 
        Give some addictional specific to source geometry.

            - if type_source is pointlike --> source position (3 elements array)
            - if type_Source is extended --> probability cumulative distribution of the elements (6 elements array)
            - if type_source is spherical --> the rsdius of the spheres (float)

    pos_max: array
        Scintillator dimension, line 5 of macro file.
    
    pos_min: array
        Mimimum position of the source ([0,0,0]) (array).
    r1: float
        random number between 0 and 1.

    r2: float
        random number between 0 and 1.

    r3: float
        random number between 0 and 1.

    r4: float
        random number between 0 and 1.


    Returns:
    -------
    pos_source: array
        position (x,y,z) of the source.

    """
    if(type_source == 'EST'): # get initial position
            face = face_func(r4, source_params)
            pos_source = source_position_est(r1,r2,r3, face, pos_max, pos_min)
    elif(type_source == 'SPH'): # get initial position
            phi_source = random_rescale(r1, np.pi)
            theta_source = random_rescale(r2, 2*np.pi)
            pos_source = from_sph_coord_to_xyz(source_params,phi_source,theta_source,pos_max[0]/2.,pos_max[1]/2.,pos_max[2]/2.) # it is centered in the centre of the system
    elif(type_source == 'PUNT'): 
            pos_source = source_params

    return pos_source

def get_initial_energy(En_type, Energy, cs_table, n, E_rand):
    """
    get_initial_energy(En_type, Energy, cs_table, n, E_rand)

    Get initial energy given input information defined in the macro file (eg. UNIF, MONO, Energy)

    Parameters:
    ----------
    En_type: string
        MONO for monoenergetic or UNIF for uniform distribution
    Energy:
        Value of energy or energy ranges of the source energy (between 1 and 98.1 MeV)
            - if the energy type is monoenergetic --> single value (float)
            - if the energy type is uniform --> array with minimum and maximum energy array([float, float])

    cs_table: 
        cross section table.

    n: float
        molecular density calculated for the plastic scintillator.

    E_rand: float
        random number between 0 and 1.
    
    Returns:
    -------
    E: float
        Energy of the neutron 
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
    if(En_type == 'UNIF') :
        E = random_rescale(E_rand, Energy[1], Energy[0])
        cs, l, p = get_cs(E,cs_table, n)
    elif(En_type == 'MONO'):
        E = Energy
        cs, l, p = get_cs(Energy, cs_table, n)

    return E, l, p



# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-


def merge_tmp_tables(table_name : str, tmp_tables_list : list):
    """
    Function that merges together the temporary files created by the 
    different processes into a single final table
    Parameters:
    ----------
        table_name : string
            filename of the final table.

        tmp_tables_list : list of string
            list of temporary tables filenames. (eg. 'tmp_file_i.txt')
    Returns:
    -------

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
# all events
def event_func(i, cs_table, data): # [cs_table, k]
    """
    event_func(i, cs_table, data)

    Main function. This function group into itself N_i single events that ill be compiled in the single process.

    Parameters:
    ----------
    i: int
        Number of the process (for multiprocessing)

    cs_table: dataframe
        cross section table.
    
    data: dictionary
        Dictionary with all information required for the montecarlo.

        {'n_processes': 
            number of processes used for multiprocessing (int) 
        'n':
            molecular density calculated for the plastic scintillator (float)
        'seed': 
            random seed (int), line 1 of macro file
        'N': 
            number of events (int), line 2 of macro file
        'En_type': 
            energy type (UNIF or MONO) (str), line 3 of macro file
        'Energy': 
            Value of energy or energy ranges of the source energy (between 1 and 98.1 MeV)
                - if the energy type is monoenergetic --> single value (float)
                - if the energy type is uniform --> array with minimum and maximum energy array([float, float])
        'pos_min': 
            Mimimum position of the source ([0,0,0]) (array)
        'pos_max': 
            Scintillator dimension, line 5 of macro file
        'type_source': str
            Source geometry (PUNT, EST or SPH) (str)
        'source_params': 
            Give some addictional specific to source geometry.

                - if type_source is pointlike --> source position (3 elements array)
                - if type_Source is extended --> probability cumulative distribution of the elements (6 elements array)
                - if type_source is spherical --> the rsdius of the spheres (float)
        }

        

    Returns:
    -------
    None

    The function produces two files: example.txt and step.txt that get track of the information of respectively each event and step information (position, energy, etc).
    See README.md file for further informations.    
    """

    # my data
    N = data['N'] # number of events
    n_processes = data['n_processes'] # number of processes
    myseed = data['seed']
    seed(myseed + i) # 'i' is used in order to get different random number for each process.

    # number of event for the i_th process
    N_i = int(N / n_processes)

    # create temporary file
    

    step_list = [] # ['event','step','x','y','z','interaction','Energy']
    event_list = [] # (columns = ['event','last_step','x','y','z'])

    for w in tqdm(range(N_i)): # run a certain number of events
        step_list, event_list = event(i, cs_table, data, N_i, w,step_list, event_list)


    step_name = "tmp_step%d.txt" % (i)
    event_name = "tmp_event%d.txt" % (i)


    with open(step_name,'w') as s:
        writer = csv.writer(s)
        writer.writerows(step_list)

    with open(event_name, 'w') as e:
        writer = csv.writer(e)
        writer.writerows(event_list)



def event(i, cs_table, data, N_i, w, step_list, event_list):
    """
    event(i, cs_table, data, N_i, w, step_list, event_list)
    Single event. Simulate the single neutron from its creation to its end. 
    The information on the event (final position, number of interactions) 
    and on the interactions (steps) the particle had on the scintillator 
    will be saved as a list and returned as an output.

    Parameters:
    ----------
    i: int
        Number of the process (for multiprocessing)

    cs_table:
        cross section table.
    
    data: dictionary
        Dictionary with all information required for the montecarlo.

        {'n_processes': 
            number of processes used for multiprocessing (int) 
        'n':
            molecular density calculated for the plastic scintillator (float)
        'seed': 
            random seed (int), line 1 of macro file
        'N': 
            number of events (int), line 2 of macro file
        'En_type': 
            energy type (UNIF or MONO) (str), line 3 of macro file
        'Energy': 
            Value of energy or energy ranges of the source energy (between 1 and 98.1 MeV)
                - if the energy type is monoenergetic --> single value (float)
                - if the energy type is uniform --> array with minimum and maximum energy array([float, float])
        'pos_min': 
            Mimimum position of the source ([0,0,0]) (array)
        'pos_max': 
            Scintillator dimension, line 5 of macro file
        'type_source': str
            Source geometry (PUNT, EST or SPH) (str)
        'source_params': 
            Give some addictional specific to source geometry.

                - if type_source is pointlike --> source position (3 elements array)
                - if type_Source is extended --> probability cumulative distribution of the elements (6 elements array)
                - if type_source is spherical --> the rsdius of the spheres (float)
        }

    N_i: int
        number of events of th esimulation (for event counting)
    w: int
        number of i_th event that has been simulated (for event counting)
        event_number = w + (i * N_i) + 1 

    step_list: list
        history of all the stpes of the previous particles in the simulation.
        In the list are written:
            [event_number, step_number, x,y,z, type_of_interaction, Energy]
        
        where x,y,z is the position of the interaction.



    event_list: list
        history of all the events of the previous particles in the simulation.
        In the list are written:
        
            [event,last_step,x,y,z]

        where x,y,z is the position where the particle stop (exit form the scintillator or reached energy lower than 1 MeV)


    Returns:
    -------
    step_list: list
        History of all the stpes of the previous particles in the simulation plus the steps of the current event.
        In the list are written:
            [event_number, step_number, x,y,z, type_of_interaction, Energy]
        
        where x,y,z is the position of the interaction.


    event_list: list
        history of all the events of the previous particles in the simulation plus the current event.
        In the list are written:
        
            [event,last_step,x,y,z]

        where x,y,z is the position where the particle stop (exit form the scintillator, reached energy lower than 1 MeV or did inelastic scattering)

    """
    e = w + (i * N_i) + 1 # number of the event

    MeV = 1000000 # 1 MeV

    n = data['n'] # molecular density
    En_type = data['En_type'] # energy type (PUNT or EST)
    Energy = data['Energy'] # source energy
    pos_min = data['pos_min'] # position of the scintillator
    pos_max = data['pos_max'] # scintillator dimension
    type_source = data['type_source'] # source geometrical distribution
    source_params = data['source_params'] # other params for the geometrical characterizzation of the source

    # get initial params
    r1 = random()
    r2 = random()
    r3 = random()
    r4 = random()
    r5 = random()
    pos_source = get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3, r4)
    E, l, p = get_initial_energy(En_type, Energy, cs_table, n, r5)

    j = 0
    x0 = pos_source[0]
    y0 = pos_source[1]
    z0 = pos_source[2]
    position = np.array([0.,0.,0.]) # initialize
    
    step_list.append([e,j,x0,y0,z0,'source',E/MeV])
    theta_scat = 0
    r = 0.
    # j_th step
    while( (position >= pos_min).all() & (position <= pos_max).all()):
        
        # 1) probability of interaction

        # it is easier to work in spherical coordinates
        if j == 0: # initial direction
            phi = random_rescale(r1, 2*np.pi) # phi is inside [0 ,2*pi ]
            theta = random_rescale(r2, np.pi)   # thetha is inside [] 0 e pi]
        else: # subsequent direction
            phi = phi + r * np.sin(theta_scat)* np.cos(r)  
            theta = theta + r * np.sin(theta_scat) * np.sin(r) 

        p_interaction = random() # probability of interaction (cs_total). It is necessary for the measure of the lenght calculated form the particle 
        r = - l * np.log(1-p_interaction) # by using the BEER LAMBERT law, i can calucate how long the particle is travelling
        
        # translate spherical coordinates into cartesian coordinates
        pos = from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0) # position of interaction (x,y,z)

        # 2) type of interaction
        r = random()

        if (r < p[0]): #elastic scattering with carbon
            A = 12
            rand = random()
            # calculate energy loss and so theta scattering of the neutron
            theta_scat, E = scattering_angle(E, A, rand)
            if E < MeV: break # energy threshold (for the energy resolution of my detector)
            j = j+1 # next step
            step_list.append([e,j,pos[0],pos[1],pos[2],'elastic',E/MeV])
        
        elif(r < p[1]): # inelastic scattering with carbon
            j = j+1 # next step
            step_list.append([e,j,pos[0],pos[1],pos[2],'inelastic',E/MeV])
            break # I'm interested onlys in multiple scattering

        elif(r < p[2]): # elastic scattering with hydrogen
            A = 1  
            # calculate energy loss and so theta scattering of the neutron
            rand = random()
            theta_scat, E = scattering_angle(E, A, rand)
            if E < MeV: break # energy thershold (for the energy resolution of my detector)
            j = j+1 # next step
            step_list.append([e,j,pos[0],pos[1],pos[2],'elastic',E/MeV])

        else: # inelastic scattering with hydrogen
            step_list.append([e,j,pos[0],pos[1],pos[2],'inelastic',E/MeV])
            break # I'm interested onlys in multiple elastic scattering
            
        position = pos
        x0,y0,z0 = pos

    # when particle is outside the scintillator
    event_list.append([e,j,pos[0],pos[1],pos[2]])

    return step_list, event_list
