# FUNCTIONS
import numpy as np
from random import random

#%% FUNZIONI da mettere in un altro file

def random_rescale(rescale):
    val = random() * rescale
    return val

def source(face, x_rect,y_rect,z_rect):
    x_source = random() * x_rect
    y_source = random() * y_rect
    z_source = random() * z_rect
    pos = [0,0,0]

    if (face == 1): pos = np.array([x_source, y_source, z_rect])
    if (face == 2): pos = np.array([x_rect, y_source, z_source])
    if (face == 3): pos = np.array([x_source, 0.0, z_source])
    if (face == 4): pos = np.array([0.0, y_source, z_source])
    if (face == 5): pos = np.array([x_source, y_rect, z_source])
    if (face == 6): pos = np.array([x_source, y_source, 0.0])

    return pos


def get_pos(f):
    "get x,y,z coords from file"
    string = f.readline()
    print(string)
    x = float(string.split(' ')[0])
    y = float(string.split(' ')[1])
    z = float(string.split(' ')[2])
    pos = np.array([x, y, z])

    return pos



def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="right")
    return idx
   

def face_func(face_prob_cum):
    r = random()
    index = find_nearest(face_prob_cum, r)
    return int(index +1)

def from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0):
    x1 = x0 + r * np.sin(phi) * np.cos(theta)
    y1 = y0 + r * np.sin(phi) * np.sin(theta)
    z1 = z0 + r * np.cos(phi)
    pos = np.array([x1,y1,z1])
    return pos

