from functions import random_rescale, source_position_est, get_cs, find_nearest, face_func, from_sph_coord_to_xyz, scattering_angle
# from make_constants import get_pos
import numpy as np
import constants as c
import os
import pandas as pd
from random import seed





def test_random_rescale_1():
    """
    Tests:
        - if output val is inside the asked range
    """

    val_min = -1
    val_max = 1
    val = random_rescale(val_max, val_min)

    assert val <= val_max
    assert val >= val_min

def test_random_rescale_2():
    """
    Tests:
        - if the value is uniformly distributed betewwn val_min and val_max
    """
    val_min = -1
    val_max = 1
    size = 10 ** 6
    bins = 20
    confidence_interval = 37.6 # for bins = 20
    val = random_rescale(val_max, val_min)
    vals = []
    for i in range(size):
        val = random_rescale(val_max, val_min)
        vals.append(val)

    observed = np.histogram(vals, bins = bins)[0]
    expected = [size/bins] * bins

    chisq = 0

    for i in range(bins):
        chisq += (expected[i] - observed[i])**2  / expected[i]

    assert chisq < confidence_interval # confidence interval of 0.99 (from tables)

def test_random_rescale_3():
    """
    Tests:
	- if, given a fixed seed, the output of random_rescale() is fixed.
    """
    seed (42)
    val_min = -1
    val_max = 1

    val = random_rescale(val_max, val_min)
    
    val = round(val, 4)

    assert val == 0.2789








def test_source_position_est_1():
    """
    Tests:
        - if, given an input face, the function place the source in the correct face of the square
    """
    pos_max = np.array([1,1,1])
    pos_min = np.array([0,0,0])

    # check for face 1
    face = 1
    pos = source_position_est(face, pos_max, pos_min)
    assert pos[2] == pos_max[2]

    # check for face 2
    face = 2
    pos = source_position_est(face, pos_max, pos_min)
    assert pos[0] == pos_max[0]

    # check for face 3 
    face = 3
    pos = source_position_est(face, pos_max, pos_min)
    assert pos[1] == pos_min[1]

    # check for face 4
    face = 4
    pos = source_position_est(face, pos_max, pos_min)
    assert pos[0] == pos_min[0]

    # check for face 5
    face = 5
    pos = source_position_est(face, pos_max, pos_min)
    assert pos[1] == pos_max[1]

    # check for face 6
    face = 6
    pos = source_position_est(face, pos_max, pos_min)
    assert pos[2] == pos_min[2]




# def test_get_pos_1():
#     """
#     Tests:
#         - if the function read properly a line on a text.txt
#     """

#     file = open("test_file.txt", "w")
    
#     file.write('1 2 3 #position')
#     file.close()

#     file = open("test_file.txt", "r")
#     pos = get_pos(file)
#     assert (pos == np.array([1.,2.,3.])).all()

#     file.close()
#     os.remove("test_file.txt")

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

def round_data(test_val, significative):
    """
    Round data to a certain significative number (significative)

    """

    for i in range(len(test_val)):
        test_val[i] = round(test_val[i], significative)
    
    return test_val




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



def test_find_nearest_1():
    """
    Test: 
        - if the output value is an int
    """

    array = [0, 0.1, 0.2, 0.3, 0.4]
    value = 0.25

    index = find_nearest(array,value)

    assert type(index) == np.int64


def test_find_nearest_2():
    """
    Test: 
        - if the output index is correct for a sample array
    """

    array = [0, 0.1, 0.2, 0.3, 0.4]
    value = 0.25

    index = find_nearest(array,value)

    assert index == 3


def test_face_func_1():
    """
    Tests:
        - if output is int
    """
    if c.type_source == 'EST':
        value = face_func()
        print (value)

        assert type(value) == int


def test_face_func_2():
    """
    Tests:
        - if value correspond to all faces number [1,2,3,4,5,6]
    """
    
    if c.type_source == 'EST':
        value = face_func()
        faces = [1,2,3,4,5,6] # faces of the rectangle
    
        assert value in faces


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



def test_scattering_angle_1():
    """
    Test:
        - if energy is in the proper range
    """

    E0 = 1
    A = 12
    Delta_E = ((A-1)/(A+1))**2*E0

    E_min = E0 - Delta_E

    theta_scat, E = scattering_angle(E0, A)

    print(theta_scat, E, E_min )

    assert E <= E0, "Energy should be <= of initial energy"
    assert E >= Delta_E, "Energy should be >= of the minimum energy"


def test_scattering_angle():
    """
    Test:
        - if theta scattering is between 0 and pi radiants
    """
    
    E0 = 1
    A = 12
    Delta_E = ((A-1)/(A+1))**2*E0

    E_min = E0 - Delta_E

    theta_scat, E = scattering_angle(E0, A)

    assert theta_scat < np.pi, "scattering angle should be smallet than pi"
    assert theta_scat > 0, "scattering angle should be larger than 0"

    
