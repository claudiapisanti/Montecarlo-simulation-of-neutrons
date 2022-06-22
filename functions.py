# FUNCTIONS
import numpy as np
from random import random

#%% FUNZIONI da mettere in un altro file

def random_rescale(val_max, val_min = 0):
    'riscalo il mio valore in un determinato range'
    val = random() * (val_max - val_min) + val_min
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
    #print(string)
    x = float(string.split(' ')[0])
    y = float(string.split(' ')[1])
    z = float(string.split(' ')[2])
    pos = np.array([x, y, z])

    return pos

def get_cs(E, n, cs_table):   
    """get csv info given Energy"""

    index = find_nearest(cs_table['energia'], E)

    cs_tot = cs_table['cs_tot'][index]
    cs_el = cs_table['cs_el'][index]
    cs_inel = cs_table['cs_inel'][index]
    cs = [cs_tot, cs_el, cs_inel] # cs[0] = totale, cs[1] = elastico, cs[2] = inelastico
    cs = np.dot(cs, 10e-24) # converto la cs da barn a cm^2
    l = np.dot(cs, n) # calcolo lambda

    p_elastic = cs[1] / cs[0] # probabilità di avere scattering elastico
    p_inelastic = cs[2] / cs[0] # probabilità di avere scattering inelastico rispetto alla probabilità di avere scattering

    return cs, l, p_elastic, p_inelastic




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

