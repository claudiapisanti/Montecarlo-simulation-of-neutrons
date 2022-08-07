import numpy as np
import multiprocessing as mp
#########################################################

def get_pos(f):
    """
    This file get form the macro file a triad of coordinates (x,y,z)

    Parameters:
    ----------
    f: file
        File must be already opened for reading (f = open('file.txt','r')) 
    
    Returns:
    ------
    pos: array 
        Position array of data from input file (x,y,z coordinates)

    """

    string = f.readline() 
    x = float(string.split(' ')[0])
    y = float(string.split(' ')[1])
    z = float(string.split(' ')[2])
    pos = np.array([x, y, z])

    return pos

def get_single_val(f):
    """
    Get single val from file sting. the separation is a void space ' '

    Parameters:
    ----------
    f: file
        File must be already opened for reading (f = open('file.txt','r')) 
    
    Returns:
    ------
    datum: string 
        Single datum form a file line
        
    """
    line = f.readline() 
    datum = line.split(' ')[0]

    return datum


def get_molecular_density(rho, Mmol):
    """
    This function returns the density of the material used. (Polivyniltoluene)
    
    Returns:
    ------
    n: float
        Density of molecules in the target
    """
    Na = 6.0221409e23 # Avogadro's number
    n = Na*rho/Mmol # density of molecules in the target

    return n 



def make_dictionary(f): 
    """
    This function create a dictionary form an input macro file (see. example.txt)

    Parameters:
    ----------
    f: file 
        File must be already opened for reading (f = open('file.txt','r')) 
    
    Returns:
    ------
    dictionary: dictionary

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
        'type_source': 
            Source geometry (PUNT, EST or SPH) (str)
        'source_params': 
            Give some addictional specific to source geometry.

                - if type_source is pointlike --> source position (3 elements array)
                - if type_Source is extended --> probability cumulative distribution of the elements (6 elements array)
                - if type_source is spherical --> the rsdius of the spheres (float)
        }

        
    """

    dictionary = {}

    MeV = 1000000

    # count number of processes available
    n_processes = mp.cpu_count()
    if n_processes != 0 : n_processes = n_processes - 1 # in order not to slow the computer too much
    dictionary['n_processes'] = n_processes

    # my material is polyvinyl toluene
    rho = 1.023 #g/cm^3 density                
    Mmol = 118. # g/mol molar mass
    n = get_molecular_density(rho, Mmol)
    dictionary['n'] = n

    # LINE 1 --> SEED (random seed)
    seed = int(get_single_val(f))
    dictionary['seed'] = seed
    
    # LINE 2 --> N (number of particles)
    N = int(get_single_val(f))
    dictionary['N'] = N

    # LINE 3 --> En_type (energy type: MONO UNIF)
    En_type = get_single_val(f) 
    dictionary['En_type'] = En_type

    # LINE 4 --> source energy
    if (En_type == 'MONO'): # monoenergetic
        Energy = float(get_single_val(f)) * MeV
    elif(En_type == 'UNIF'): # uniform distribution
        line = f.readline() 
        Energy = np.array([float(line.split(' ')[0]), float(line.split(' ')[1])]) * MeV# [pos_min, pos_max]
    dictionary['Energy'] = Energy

    # LINE 5 --> Scintillator parameters + minimum position (0,0,0)
    pos_max = get_pos(f)
    pos_min = np.array([0,0,0])
    dictionary['pos_min'] = pos_min
    dictionary['pos_max'] = pos_max

    # LINE 6 --> Source Geometry (PUNT, EST, SPH)
    type_source = get_single_val(f)
    dictionary['type_source'] = type_source

    # LINE 7 --> Source geometry parameters
    if (type_source == 'PUNT'): # pointlike
        source_params = get_pos(f) # in PUNT, the source parameter is the position of the source
    elif(type_source == 'EST'): # extended
        face_prob = np.zeros(6)
        face_prob_cum = np.zeros(6)
        face_prob_line = f.readline() 
        
        for i in range(6): # get faces probability distribution from file
            face_prob[i] = float(face_prob_line.split(' ')[i])


        for i in range(6): # calculate fsce probability cumulative
            face_prob_cum[i] = sum(face_prob[:(i+1)]) / sum(face_prob) # cumulative function  of face distribution probability (it is a stepped function)

        source_params = face_prob_cum

    elif(type_source == 'SPH'): # spherical
        source_params = float(get_single_val(f)) # in EST, the source parameter is the radius of the scintillator
    dictionary['source_params'] = source_params

    return dictionary
