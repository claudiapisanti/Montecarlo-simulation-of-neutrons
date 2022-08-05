from errno import ENETRESET
from typing import TYPE_CHECKING
from functions import random_rescale, source_position_est, get_cs, find_nearest, face_func, from_sph_coord_to_xyz, scattering_angle
from functions import get_source_position, get_initial_energy, merge_tmp_tables
from functions import event_func, event
from make_constants import get_pos, get_single_val, get_molecular_density, make_dictionary
import numpy as np
import constants as c
import os
import pandas as pd
from random import seed

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# UTILS

def round_data(test_val, significative):
    """
    Round data to a certain significative number (significative)

    """

    for i in range(len(test_val)):
        test_val[i] = round(test_val[i], significative)
    
    return test_val

def my_dataframe_and_energy():
    cs = {'E': [0.9,1.1], 
          'cs_tot': [2,3],
          'cs_el' : [3,4],
          'cs_inel' : [4,5],
          'cs_h_tot' : [5,6],
          'cs_h_el': [5,6],
          'cs_h_inel': [6,7],
          'cs_c_tot' : [7,8],
          'cs_c_el' : [9,10],
          'cs_c_inel' : [10,11] }
    E = 1
    cs_table = pd.DataFrame(cs)

    return cs_table, E


# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# TEST FUNCTIONS.PY 
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_random_rescale_1():
    """
    Tests:
        - if rand = 1, val is equal to val max
    """

    rand = 1
    val_min = -1
    val_max = 1
    val = random_rescale(rand, val_max, val_min)
    print(val)
    assert val == val_max
    
def test_random_rescale_2():
    """
    Tests:
        - if rand = 0, val is equal to val min
    """

    rand = 0
    val_min = -1
    val_max = 1
    val = random_rescale(rand, val_max, val_min)
    print(val)

    assert val == val_min

def test_random_rescale_3():
    """
    Tests:
        - if rand = 0.5, val is equal to 0.
    """
    rand = 0.5
    val_min = -1
    val_max = 1
    val = random_rescale(rand, val_max, val_min)
    print(val)

    assert val == 0.0

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_source_position_est_1():
    """
    Tests:
        - if, given an input face, the function place the source in the correct face of the square
    """
    r1 = 1.
    r2 = 1.
    r3 = 1.

    pos_max = np.array([1,1,1])
    pos_min = np.array([0,0,0])

    # check for face 1
    face = 1
    pos = source_position_est(r1,r2,r3,face, pos_max, pos_min)
    assert pos[2] == pos_max[2]

    # check for face 2
    face = 2
    pos = source_position_est(r1,r2,r3,face, pos_max, pos_min)
    assert pos[0] == pos_max[0]

    # check for face 3 
    face = 3
    pos = source_position_est(r1,r2,r3,face, pos_max, pos_min)
    assert pos[1] == pos_min[1]

    # check for face 4
    face = 4
    pos = source_position_est(r1,r2,r3,face, pos_max, pos_min)
    assert pos[0] == pos_min[0]

    # check for face 5
    face = 5
    pos = source_position_est(r1,r2,r3,face, pos_max, pos_min)
    assert pos[1] == pos_max[1]

    # check for face 6
    face = 6
    pos = source_position_est(r1,r2,r3,face, pos_max, pos_min)
    assert pos[2] == pos_min[2]

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_get_cs_1():
    """
    Tests:
        - if, given a proper csv, the cross section in output is as expected
    """
    n =  5.220889949745763e+21
    # import my sample dataframe
    cs_table, E = my_dataframe_and_energy()
    test_cs, test_l, test_p = get_cs(E, cs_table, n)
    
    # round data to the first devimal number
    test_cs = np.dot(test_cs,10**23) # re-convert data into integer (from e-23 to unit) ro to easily round my data
    test_cs = round_data(test_cs, 1)
    

    assert (test_cs == cs_table.iloc[1][1:]).all()


def test_get_cs_2():
    """
    Tests:
        - if, the output of the free mean path (l) is read correctly
    """
    n =  5.220889949745763e+21

    cs_table, E = my_dataframe_and_energy()
    test_cs,test_l,test_p = get_cs(E, cs_table, n)


    test_l = round_data(test_l, 2)

    l = [0.16, 0.21,  0.26,  0.31,  0.31,  0.37, 0.42,  0.52, 0.57] 

    assert (test_l == l).all()

def test_get_cs_3():
    """
    Tests:
        - if, the output of the probability (p) is read correctly
    """
    n =  5.220889949745763e+21

    cs_table, E = my_dataframe_and_energy()
    test_cs,test_l,test_p = get_cs(E, cs_table, n)

    test_p = round_data(test_p, 2)

    p = [24.0, 20.0, 1.0, 1.17, 1.25, 1.38]

    assert test_p == p

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_find_nearest_1():
    """
    Test: 
        - if the output index is correct for a sample array
    """

    array = [0, 0.1, 0.2, 0.3, 0.4]
    value = 0.25

    index = find_nearest(array,value)

    assert index == 3

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_face_func_1():
    """
    Tests:
        - if value correspond to all faces number [1,2,3,4,5,6]
    """
    r = 0.51
    array_cumulative = np.array([0.1,0.2,0.3,0.4,0.5,1.]) # cumulative distribution

    value = face_func(r, array_cumulative)
    faces = [1,2,3,4,5,6] # faces of the rectangle
    assert value == 6

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_from_sph_coord_to_xyz_1():
    """
    Tests:
        - if the function properly convert a spherical coordinate into x,y,z, coordinates
    """

    # sample initial spherical coordinates
    r = 1
    phi = np.pi /6. #30°
    theta = np.pi /6. #30°
    x0 = 0.
    y0 = 0.
    z0 = 0.

    xyz = from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0)
    
    # round data to the 2nd significative cifer
    xyz = round_data(xyz,2)
    
    assert (xyz == [0.43, 0.25, 0.87]).all()

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_scattering_angle_1():
    """
    Test:
        - if rand = 0, E = E_min
    """
    r = 0.
    E0 = 1
    A = 12
    Delta_E = ((A-1)/(A+1))**2*E0

    theta_scat, E = scattering_angle(E0, A, r)

    E = round(E, 2)
    Delta_E = round(Delta_E, 2)
    assert E == Delta_E

def test_scattering_angle_2():
    """
    Test:
        - if rand = 1, E = E_max
    """
    r = 1.
    E0 = 1
    A = 12

    theta_scat, E = scattering_angle(E0, A, r)
    theta_scat = round(theta_scat, 2)

    assert E == E0

def test_scattering_angle_3():
    """
    Test: 
       - if rand = 1, theta scattering = 0
    """
    r = 1.
    E0 = 1
    A = 12

    theta_scat, E = scattering_angle(E0, A, r)
    theta_scat = round(theta_scat, 2)

    assert theta_scat == 0.


def test_scattering_angle_4():
    """
    Test:
       - if rand = 0, theta_scattering = pi
    """
    r = 0.
    E0 = 1
    A = 12

    theta_scat, E = scattering_angle(E0, A, r)
    theta_scat = round(theta_scat, 2)

    assert theta_scat == 3.14

def test_scattering_angle_5():
    """
    Test:
       - if rand = 0.5, theta_scattering is equal to the theoretical one
    """
    r = 0.5
    E0 = 1
    A = 12

    theta_scat, E = scattering_angle(E0, A, r)
    theta_scat = round(theta_scat, 2)

    cos_theta_CM = ((E) / E0 * (A+1)**2 - A**2 - 1)/(2*A) # see LEO
    Theta_sc = np.arccos((A * cos_theta_CM + 1)/ np.sqrt(A**2 + 1 + 2*A*cos_theta_CM))
    Theta_sc = round(Theta_sc, 2)
    
    assert Theta_sc==theta_scat

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_get_source_position_1():
    """
    Test:
       - if the source is 'PUNT', given the proper information the output is as expected.
    """
    type_source = 'PUNT'
    source_params = np.array([1.,1.,1.])
    pos_max = np.array([2.,2.,2.])
    pos_min = np.array([0.,0.,0.])
    r1 = 1.
    r2 = 1.
    r3 = 1.
    r4 = 1.

    pos_source = get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3,r4)

    assert (pos_source == np.array([1.,1.,1.])).all()


def test_get_source_position_2():
    """
    Test:
       - if the source is 'EST', given the proper information the output is as expected.
    """
    type_source = 'EST'
    source_params = np.array([0.1,0.2,0.3,0.4,0.5,1.])
    pos_max = np.array([2.,2.,2.])
    pos_min = np.array([0.,0.,0.])
    r1 = 1.
    r2 = 1.
    r3 = 1.
    r4 = 0.99

    pos_source = get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3,r4)

    assert (pos_source == np.array([2.,2.,0.])).all()

def test_get_source_position_3():
    """
    Test:
       - if the source is 'SPH', given the proper information the output is as expected.
    """
    type_source = 'SPH'
    source_params = 1.
    pos_max = np.array([2.,2.,2.])
    pos_min = np.array([0.,0.,0.])
    r1 = 1.
    r2 = 1.
    r3 = 1.
    r4 = 1.

    pos_source = get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3,r4)
    pos_source = round_data(pos_source, 1)

    assert ( pos_source == np.array([1.,1.,2.])).all()

def test_get_source_position_4():
    """
    Test: 
       - if source is est, the source position is in a sphere of radius as in source params.
    """
    type_source = 'SPH'
    source_params = 1.
    pos_max = np.array([2.,2.,2.])
    pos_min = np.array([0.,0.,0.])
    r1 = 1.
    r2 = 1.
    r3 = 1.
    r4 = 1.

    pos_source = get_source_position(type_source, source_params, pos_max, pos_min, r1,r2,r3,r4)
    pos_source = round_data(pos_source, 1)

    centre = pos_max / 2.
    calculated_radius =(pos_source[0]- centre[0])**2 + (pos_source[1]- centre[1])**2 + (pos_source[2]- centre[2])**2 
    assert calculated_radius == source_params

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_get_initial_energy_1():
    """
    Test:
        - if, given proper informationo, the initial energy is as expecetd
    """

    En_type = 'MONO'
    cs_table, Energy = my_dataframe_and_energy()
    n = 1.
    E_rand = 1.
    E, l, p = get_initial_energy(En_type, Energy, cs_table, n, E_rand)

    l = l*10**23
    l = round_data(l, 1)
    p = round_data(p, 1)

    assert E == Energy, "Check errors in energy calculation"
    assert (l == np.array([3.0, 4.0, 5.0, 6.0, 6.0, 7.0, 8.0, 10.0, 11.0])).all(), "check errors in free mean path (l)"
    assert (p == np.array([24.0, 20.0, 1.0, 1.2, 1.2, 1.4])).all(), "Check errors in probability od having carbon / hydrogeb scattering (p)"

def test_get_initial_energy_2():
    """
    Test:
        - if, given proper information, the initial energy is as expecetd
    """

    En_type = 'UNIF'
    cs_table, Energy = my_dataframe_and_energy()
    Energy = np.array([0.8, 1.])
    n = 1.
    E_rand = 1.
    E, l, p = get_initial_energy(En_type, Energy, cs_table, n, E_rand)

    l = l*10**23
    l = round_data(l, 1)
    p = round_data(p, 1)
    assert E == 1., "Check errors in energy calculation"
    assert (l == np.array([3.0, 4.0, 5.0, 6.0, 6.0, 7.0, 8.0, 10.0, 11.0])).all(), "check errors in free mean path (l)"
    assert (p == np.array([24.0, 20.0, 1.0, 1.2, 1.2, 1.4])).all(), "Check errors in probability od having carbon / hydrogeb scattering (p)"

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_merge_tmp_tables():
    """
    Test:
        - if file correctly merge two tables
    """
    # generate temporaty files
    table_name = 'test.txt'
    tmp_tables_list =['tmp_test_1.txt', 'tmp_test_2.txt']

    for file_name in tmp_tables_list:
        print(file_name)
        file = open(file_name, "w")
        file.write('bla ')
        file.close()

    # my tested function
    merge_tmp_tables(table_name, tmp_tables_list)

    # check if the content of my merged file correspond to what I wrote
    file = open(table_name, "r")
    words = file.readlines()

    # remove temporary files
    for file_name in tmp_tables_list:
        os.remove(file_name)
    
    os.remove(table_name)

    assert words == ['bla bla ']

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
def test_event_func():
    """
    Test: 
        - if the function create a file as expected
    """
    i = 0
    cs_table, Energy = my_dataframe_and_energy()
    data = {
        'n_processes' : 1, # number of processes
        'n':1, # n
        'seed': 42, # seed
        'N':10, # N
        'En_type':'MONO', # En_type
        'Energy': Energy, # Energy
        'pos_min':np.array([0.,0.,0.]), # pos_min
        'pos_max':np.array([2.,2.,2.]), # pos_max
        'type_source':'PUNT', # type_source
        'source_params':np.array([1.,1.,1.]), # pos_min
    }

    # MY FUNCTION TO BE TESTED
    event_func(i, cs_table, data)

    step_name = "tmp_step%d.txt" % (i)
    event_name = "tmp_event%d.txt" % (i)

    # step_file = open(step_name, 'r')
    # event_file = open(event_name, 'r')

    with open(step_name,'r') as s:
        string1 = s.readline()
        print(string1)
        assert string1 == '1,0,1.0,1.0,1.0,source,1\n', "Check problems in step writing"

    with open(event_name, 'r') as e:
        string2 = e.readlines()
        assert string2 == ['1,0,1.0,1.0,1.0\n'], "Check probelms with event writing"


    os.remove(step_name)
    os.remove(event_name)

def test_event():
    """
    Test:
        - if the single event give the expected results.
    """
    
    i = 0
    cs_table, Energy = my_dataframe_and_energy()
    data = {
        'n_processes' : 1, # number of processes
        'n':1, # n
        'seed': 42, # seed
        'N':10, # N
        'En_type':'MONO', # En_type
        'Energy': Energy, # Energy
        'pos_min':np.array([0.,0.,0.]), # pos_min
        'pos_max':np.array([2.,2.,2.]), # pos_max
        'type_source':'PUNT', # type_source
        'source_params':np.array([1.,1.,1.]), # pos_min
    }

    N_i = data['N']
    w = 1



    step_list = []
    event_list = []
    step_list, event_list = event(i, cs_table, data, N_i, w, step_list, event_list)

    assert step_list == [[2, 0, 1.0, 1.0, 1.0, 'source', 1]]
    assert event_list == [[2, 0, 1.0, 1.0, 1.0]]
    

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# MAKE_CONSTANTS FILE
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_get_pos_1():
    """
    Tests:
        - if the function read properly a line on a text.txt
    """

    file = open("test_file.txt", "w")
    
    file.write('1 2 3 #position')
    file.close()

    file = open("test_file.txt", "r")
    pos = get_pos(file)
    assert (pos == np.array([1.,2.,3.])).all()

    file.close()
    os.remove("test_file.txt")

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

def test_get_single_val():
    """
    Test:
        - if the function read properly a value form an input file.
    """

    with open('tmp_file.txt', 'w') as f:
        f.write('bla bla bla')

    with open('tmp_file.txt', 'r') as f:
        val = get_single_val(f)
    os.remove('tmp_file.txt')

    assert val == 'bla'

def test_get_molecular_density():
    """
    Test:
        - if molecular density respect the formula 
    """
    rho = 1.
    Mmol = 1.

    n = get_molecular_density(rho, Mmol)
    n = n * 10 ** (-23)
    n = round(n, 3)

    assert n == 6.022

def test_make_dictionary():
    """
    Test:
        - if, given an input macro file (dictionary_test.txt), the function properly get all the information contained in it.
    """
    
    with open('dictionary_test.txt', 'r') as test:
        dd = make_dictionary(test)

    assert dd['seed'] == 42
    assert dd['N'] == 100
    assert dd['En_type'] == 'MONO'
    assert round(dd['Energy'], 0) == 1000000
    assert (dd['pos_min' ] == np.array([0.,0.,0.])).all()
    assert (dd['pos_max'] == np.array([2.,2.,2.])).all()
    assert dd['type_source'] == 'SPH'
    assert dd['source_params'] == 1.


