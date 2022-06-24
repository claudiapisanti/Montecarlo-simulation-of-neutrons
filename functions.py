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

    cs_h_tot = cs_table['cs_h_tot'][index]
    cs_h_el = cs_table['cs_h_el'][index]
    cs_h_inel = cs_table['cs_h_inel'][index]

    cs_c_tot = cs_table['cs_c_tot'][index]
    cs_c_el = cs_table['cs_c_el'][index]
    cs_c_inel = cs_table['cs_c_inel'][index]


    cs = [cs_tot, cs_el, cs_inel, 
          cs_h_tot, cs_h_el, cs_h_inel, 
          cs_c_tot, cs_c_el, cs_c_inel] # cs[0] = totale, cs[1] = elastico, cs[2] = inelastico
    cs = np.dot(cs, 10e-24) # converto la cs da barn a cm^2
    l = np.dot(cs, n) # calcolo lambda

    p_carbon = 9 * cs[6] / cs[0]
    p_proton = 10 * cs[3] / cs[0]

    p_h_elastic = cs[4] / cs[3] # probabilità di avere scattering elastico
    p_h_inelastic = cs[5] / cs[3] # probabilità di avere scattering inelastico rispetto alla probabilità di avere scattering

    p_c_elastic = cs[7] / cs[6] # probabilità di avere scattering elastico
    p_c_inelastic = cs[8] / cs[6] # probabilità di avere scattering inelastico rispetto alla probabilità di avere scattering

    p = [p_carbon, p_proton,
         p_h_elastic, p_h_inelastic,
         p_c_elastic, p_c_inelastic]
    return cs, l, p




def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="right")
    return idx
   

def face_func(face_prob_cum):
    r = random()
    index = find_nearest(face_prob_cum, r)
    return int(index +1)

def from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0): 
    # ANGOLO IN RADIANTI!!!!
    
    x1 = x0 + r * np.cos(phi) * np.sin(theta)
    y1 = y0 + r * np.sin(phi) * np.sin(theta)
    z1 = z0 + r * np.cos(theta)
    pos = np.array([x1,y1,z1])
    return pos

def scattering_angle(E, A):
    """calcolo l'angolo di scattering del neutrone rispetto alla direzione di origine.
    Formule by LEO<3"""

    DeltaE = random_rescale(E, ((A - 1)/(A+1))**2*E) # formula che sta sul LEO                     
    cos_theta_CM = ((E - DeltaE) / E * (A+1)**2 - A**2 - 1)/(2*A)
    theta_scat = (A * cos_theta_CM + 1)/ np.sqrt(A**2 + 1 + 2*A*cos_theta_CM)
    return theta_scat, E