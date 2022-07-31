import numpy as np
import multiprocessing as mp
import sys
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





# #########################################################
# # limit value for energy (maximum and minimum value of cs.txt)
# Maximum_energy = 98100000
# Minimum_energy = 1100000

# MeV1 = 1000000 # = 1 MeV


# # inizializing
# E_mono = 0
# E_max = 0
# E_min = 0
# face_prob_cum = 0
# sph_radius = 0
# pos_source = [0.,0.,0.]



# f = open('example.txt')

# mydict = make_dictionary(f)
# print(mydict)


# ####################################################
# # number of particles generated
# line = f.readline() # line 0 -- number of particles

# N = int(line.split(' ')[0])

# assert type(N) == int, 'number of particles should be an integer number. Correct physical_caracteristics.txt'




# # select energy distribution type (MONO, UNIF) = monoenergetic, uniformly distributed in a range
# line  = f.readline() # line 1 --- type of source
# En_type = line.split(' ')[0]
# if(En_type == 'MONO'): # monoenergetic source
#     line = f.readline() # in eV # line 2 --- energy of source
#     E_mono = float(line.split(' ')[0]) # get energy of particles form physical_characteristics file

#     # check if E_mono is inside the proper range
#     assert E_mono < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Correct physical_caracteristics.txt"
#     assert E_mono > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Correct physical_caracteristics.txt"





# elif(En_type == 'UNIF'):
#     line = f.readline() # in eV  # line 2 -- energy of source
#     E_min = float(line.split(' ')[0])
#     E_max = float(line.split(' ')[1])

#     # check if E_max is inside the proper range
#     assert E_max < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Please, correct physical_caracteristics.txt"
#     assert E_max > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Please, correct physical_caracteristics.txt"
#     # check if E_min is inside the proper range
#     assert E_max < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Please, correct physical_caracteristics.txt"
#     assert E_max > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Please, correct physical_caracteristics.txt"

#     assert E_max > E_min, "Maximum energy should be larger than minimu energy. Please, correct physical_caracteristics.txt"

# else:
#     print('Error: energy type not valid. Please, correct physical_caracteristics.txt')
#     quit()




# # get rectangle dimensions (scintillator)
# pos_max = get_pos(f) # read line 3 --- geometrical boundaries of the scintillator
# pos_min = np.array([0,0,0])


# ###########################################################
# ###########################################################
# # define type of source distribution (PUNT, SPH, EST)= point, spherical, extended
# type_source = f.readline().split(' ')[0] # line 4 --- type of source : EST, PUNT, SPH



# # souce position

# if (type_source == 'PUNT'): # pointlike
#     # get source position
#     pos_source = get_pos(f) # line 5 --- position of point source
#     face = 0

#     while((pos_source >= pos_max ).any()):
#         print ('errore: il la sorgente è fuori dallo scintillatore. Correggi il file.txt')
#         quit()

# elif(type_source == 'EST'): #  extended source (= a rectanglular surface around the scintillator)

#     # choice of the face in which the particle is generated
#     face_prob = np.zeros(6)
#     face_prob_cum = np.zeros(6)
#     face_prob_line = f.readline() # line 5 --- face probabilities

#     for i in range(6):
#        face_prob[i] = float(face_prob_line.split(' ')[i])


#     for i in range(6):
#         face_prob_cum[i] = sum(face_prob[:(i+1)]) / sum(face_prob) # cumulative function  of face distribution probability (it is a stepped function)

# elif(type_source == 'SPH'): # spherical
#     sph_radius = float(f.readline().split(' ')[0]) # (cm) 
#     face = 0 # for step.txt and event.txt filling 



# # -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# # -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# # -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# # # WRITE ON CONSTANTS.PY
# # w.write("import numpy as np\n\n")

# # w.write(f"N = {N}\n")
# # w.write(f"n_processes = {n_processes}\n")
# # w.write(f"type_source = '{type_source}'\n")
# # w.write(f"En_type = '{En_type}'\n")
# # w.write(f'E_mono = {E_mono}\n')
# # w.write(f'E_max = {E_max}\n')
# # w.write(f'E_min = {E_min}\n')
# # w.write(f'n = {n}\n')

# # if type_source == 'EST': w.write(f'face_prob_cum = np.array([{face_prob_cum[0]},{face_prob_cum[1]},{face_prob_cum[2]},{face_prob_cum[3]},{face_prob_cum[4]},{face_prob_cum[5]}])\n')
# # else:  w.write(f'face_prob_cum = {face_prob_cum}\n' )

# # w.write(f'pos_max = np.array([{pos_max[0]},{pos_max[1]},{pos_max[2]}])\n')
# # w.write(f'pos_min = np.array([{pos_min[0]},{pos_min[1]},{pos_min[2]}])\n')
# # w.write(f'sph_radius = {sph_radius}\n')
# # w.write(f'pos_source = np.array([{pos_source[0]},{pos_source[1]},{pos_source[2]}])\n')
# # w.write(f'MeV1 = {MeV1}')


# # w.close()
# f.close()

