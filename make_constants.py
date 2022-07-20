import numpy as np
import multiprocessing as mp
import sys
#########################################################

def get_pos(f):
    """
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

    


#########################################################
# get file
macro_file = sys.argv[1]

# open files
w = open('constants.py', 'w')
f = open(macro_file, "r") 

# limit value for energy (maximum and minimum value of cs.txt)
Maximum_energy = 98100000
Minimum_energy = 1100000

# inizializing
E_mono = 0
E_max = 0
E_min = 0
face_prob_cum = 0
sph_radius = 0
pos_source = [0.,0.,0.]





Na = 6.0221409e23 # Avogadro's number
MeV1 = 1000000 # = 1 MeV

# my material is polyvinyl toluene
rho = 1.023 #g/cm^3 density                
Mmol = 118. # g/mol molar mass
n = Na*rho/Mmol # density of molecules in the target

# number of particles generated
line = f.readline() # line 0 -- number of particles
N = int(line.split(' ')[0])
assert type(N) == int, 'number of particles should be an integer number. Correct physical_caracteristics.txt'

# count number of processes available
n_processes = mp.cpu_count()
if n_processes != 0 : n_processes = n_processes - 1 # in order not to slow the computer too much



# select energy distribution type (MONO, UNIF) = monoenergetic, uniformly distributed in a range
line  = f.readline() # line 1 --- type of source
En_type = line.split(' ')[0]
if(En_type == 'MONO'): # monoenergetic source
    line = f.readline() # in eV # line 2 --- energy of source
    E_mono = float(line.split(' ')[0]) # get energy of particles form physical_characteristics file

    # check if E_mono is inside the proper range
    assert E_mono < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Correct physical_caracteristics.txt"
    assert E_mono > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Correct physical_caracteristics.txt"





elif(En_type == 'UNIF'):
    line = f.readline() # in eV  # line 2 -- energy of source
    E_min = float(line.split(' ')[0])
    E_max = float(line.split(' ')[1])

    # check if E_max is inside the proper range
    assert E_max < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Please, correct physical_caracteristics.txt"
    assert E_max > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Please, correct physical_caracteristics.txt"
    # check if E_min is inside the proper range
    assert E_max < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Please, correct physical_caracteristics.txt"
    assert E_max > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Please, correct physical_caracteristics.txt"

    assert E_max > E_min, "Maximum energy should be larger than minimu energy. Please, correct physical_caracteristics.txt"

else:
    print('Error: energy type not valid. Please, correct physical_caracteristics.txt')
    quit()




# get rectangle dimensions (scintillator)
pos_max = get_pos(f) # read line 3 --- geometrical boundaries of the scintillator
pos_min = np.array([0,0,0])


###########################################################
###########################################################
# define type of source distribution (PUNT, SPH, EST)= point, spherical, extended
type_source = f.readline().split(' ')[0] # line 4 --- type of source : EST, PUNT, SPH



# souce position

if (type_source == 'PUNT'): # pointlike
    # get source position
    pos_source = get_pos(f) # line 5 --- position of point source
    face = 0

    while((pos_source >= pos_max ).any()):
        print ('errore: il la sorgente è fuori dallo scintillatore. Correggi il file.txt')
        quit()

elif(type_source == 'EST'): #  extended source (= a rectanglular surface around the scintillator)

    # choice of the face in which the particle is generated
    face_prob = np.zeros(6)
    face_prob_cum = np.zeros(6)
    face_prob_line = f.readline() # line 5 --- face probabilities

    for i in range(6):
       face_prob[i] = float(face_prob_line.split(' ')[i])


    for i in range(6):
        face_prob_cum[i] = sum(face_prob[:(i+1)]) / sum(face_prob) # cumulative function  of face distribution probability (it is a stepped function)

elif(type_source == 'SPH'): # spherical
    sph_radius = float(f.readline().split(' ')[0]) # (cm) 
    face = 0 # for step.txt and event.txt filling 



# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# WRITE ON CONSTANTS.PY
w.write("import numpy as np\n\n")

w.write(f"N = {N}\n")
w.write(f"n_processes = {n_processes}\n")
w.write(f"type_source = '{type_source}'\n")
w.write(f"En_type = '{En_type}'\n")
w.write(f'E_mono = {E_mono}\n')
w.write(f'E_max = {E_max}\n')
w.write(f'E_min = {E_min}\n')
w.write(f'n = {n}\n')

if type_source == 'EST': w.write(f'face_prob_cum = np.array([{face_prob_cum[0]},{face_prob_cum[1]},{face_prob_cum[2]},{face_prob_cum[3]},{face_prob_cum[4]},{face_prob_cum[5]}])\n')
else:  w.write(f'face_prob_cum = {face_prob_cum}\n' )

w.write(f'pos_max = np.array([{pos_max[0]},{pos_max[1]},{pos_max[2]}])\n')
w.write(f'pos_min = np.array([{pos_min[0]},{pos_min[1]},{pos_min[2]}])\n')
w.write(f'sph_radius = {sph_radius}\n')
w.write(f'pos_source = np.array([{pos_source[0]},{pos_source[1]},{pos_source[2]}])\n')
w.write(f'MeV1 = {MeV1}')


w.close()
f.close()

