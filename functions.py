# FUNCTIONS
from logging import FileHandler
import numpy as np
from random import random
import constants as c

#%% FUNZIONI da mettere in un altro file

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
          cs_c_tot, cs_c_el, cs_c_inel] # cs[0] = totale, cs[1] = elastico, cs[2] = inelastico
    cs = np.dot(cs, 10e-24) # converto la cs da barn a cm^2
    l = np.dot(cs, c.n) # calcolo lambda

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
    return idx
   

def face_func():
    """
    Randomly select a face of the rectangular extended source

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







# def evento(i, step, event, type_source, En_type, 
# E_mono, E_max, E_min, n,
# face_prob_cum, pos_max, pos_min, sph_radius, 
# cs_table, pos_source, face, sub_rect, k,
# MeV1, 
# ):
#     """run sul singolo evento"""
#     print('evento = ', i)
#     if(type_source == 'EST'): # ottengo una posizione iniziale
#         face = face_func(face_prob_cum)
#         pos_source = source(face,pos_max[0],pos_max[1],pos_max[2])
#     elif(type_source == 'SPH'): # ottengo una posizione finale
#         phi_source = random_rescale(np.pi)
#         theta_source = random_rescale(2*np.pi)
#         pos_source = from_sph_coord_to_xyz(sph_radius,phi_source,theta_source,pos_max[0]/2.,pos_max[1]/2.,pos_max[2]/2.) # centrato al centro del sistema
    
#     if(En_type == 'UNIF') :
#         E = random_rescale(E_max, E_min)
#         cs, l, p = get_cs(E,n,cs_table)
#     elif(En_type == 'MONO'):
#         E = E_mono
#         cs, l, p = get_cs(E_mono,n,cs_table)



#     j = 0
#     p_type = 0.
#     x0 = pos_source[0]
#     y0 = pos_source[1]
#     z0 = pos_source[2]
#     pos = [0.,0.,0.] # ?? inizializzo
    
#     step.write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\tsource\n') # sorgente

#     DeltaE = 0
#     theta_scat = 0
#     while( (pos >= pos_min ).all()  & (pos <= pos_max).all()):
       
#         # mi conviene lavorare in coordinate polari

#         # direzione iniziale
#         if j == 0:
#             phi = random_rescale(2*np.pi) # phi è compreso tra 0 e 2*pi 
#             theta = random_rescale(np.pi)   # thetha è compreso tra 0 e pi
#         else: 
#             rand = random_rescale(1, -1)
#             phi = phi + r * np.sin(theta_scat)* np.cos(r)  
#             theta = theta + r * np.sin(theta_scat) * np.sin(r) 
#         # percorso fatto dal neutrone
#         p_interaction = random() # probabilità di interazione con cui calcolare lo spazio
#         r = - l[0] * np.log(1-p_interaction)# sto usando la BEER LAMBERT LAW ma non so se posso usarla per i fotoni # uso lambda della cross section totale
        
#         # traduco le coordinate sferuche in coordinate cartesiane (x,y,z)
#         pos = from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0) # posizione dell'interazione
        
#         # controllo di stare dentro il rettangolo
#         if(  (pos >= pos_min ).all()  & (pos <= pos_max).all() ):
 
#             # segno quante interazioni accadono in un subrettangolo
#             if( (pos >= pos_min ).all() & (pos <= sub_rect).all() ): k = k+1
                
#             # vedo se il protone interagisce con un carbonio o con un protone
#             p_atom = random()
#             if (p_atom <= p[0]): # allora interagisce con il carbonio
#                 # vedo se fa urto elastico o inelastico
#                 p_type = random() 
            
#                 if(p_type <= p[4] ): # ovvero se sono all'interno di uno scattering ELASTICO con il carbonio
#                     A = 12
#                     step.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\n')
                    
#                     # calcolo la perdita di energia e di conseguenza il theta_sc (da capire)

#                     theta_scat, E = scattering_angle(E, A)
#                     if E < MeV1: break # QUESTA SOGLIA VA ABBASSATA (?)
#                     j = j+1 # step successivo

                
#                 else: # ovvero se il neutrone fa scattering inelastico
#                     step.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\n')
#                     event.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
#                     j = j+1 # step successivo

#                     break

#             else: # allora interagisce con il protone
#                 # vedo se fa urto elastico o inelastico
#                 p_type = random() 
            
#                 if(p_type <= p[2] ): # ovvero se sono all'interno di uno scattering elastico con il protone
#                     A = 1
#                     step.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\n')
                    
                    
#                     # calcolo la perdita di eienrgia e di conseguenza il theta_sc (da capire)
#                     theta_scat, E = scattering_angle(E, A)
#                     if E < MeV1: break # QUESTA SOGLIA VA ABBASSATA (?)

#                     j = j+1 # step successivo

                
#                 else: # ovvero se il neutrone fa scattering inelastico
#                     step.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\n')
#                     event.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
                     
#                     j = j+1 # step successivo

#                     break
                
                
#         elif((pos!=pos_source).all()):
#             # print('particle out of detector')
#             event.write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\n')
     
#         x0,y0,z0 = pos
        
#     return k



#### PROVO A VEDERE SE, LEGGENDO DI VOLTA IN VOLTA LE COSE DAL FILE.C MI ESCE MEGLIO IL FILE.

def evento_letto_da_const(i, args): # [step, event, cs_table, k]


    """run sul singolo evento"""
    print('evento = ', i)
    if(c.type_source == 'EST'): # ottengo una posizione iniziale
        face = face_func()
        pos_source = source_position_est(face)
    elif(c.type_source == 'SPH'): # ottengo una posizione finale
        face = 0
        phi_source = random_rescale(np.pi)
        theta_source = random_rescale(2*np.pi)
        pos_source = from_sph_coord_to_xyz(c.sph_radius,phi_source,theta_source,c.pos_max[0]/2.,c.pos_max[1]/2.,c.pos_max[2]/2.) # centrato al centro del sistema
    elif(c.type_source == 'PUNT'): 
        face = 0
        pos_source = c.pos_source
    
    if(c.En_type == 'UNIF') :
        E = random_rescale(c.E_max, c.E_min)
        cs, l, p = get_cs(E,args[2])
    elif(c.En_type == 'MONO'):
        E = c.E_mono
        cs, l, p = get_cs(c.E_mono, args[2])


    
    j = 0
    p_type = 0.
    x0 = pos_source[0]
    y0 = pos_source[1]
    z0 = pos_source[2]
    pos = np.array([0.,0.,0.]) # ?? inizializzo
    
    args[0].write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\tsource\n') # sorgente

    DeltaE = 0
    theta_scat = 0
    while( (pos >= c.pos_min).all() & (pos <= c.pos_max).all()):
       
        # mi conviene lavorare in coordinate polari

        # direzione iniziale
        if j == 0:
            phi = random_rescale(2*np.pi) # phi è compreso tra 0 e 2*pi 
            theta = random_rescale(np.pi)   # thetha è compreso tra 0 e pi
        else: 
            rand = random_rescale(1, -1)
            phi = phi + r * np.sin(theta_scat)* np.cos(r)  
            theta = theta + r * np.sin(theta_scat) * np.sin(r) 
        # percorso fatto dal neutrone
        p_interaction = random() # probabilità di interazione con cui calcolare lo spazio
        r = - l[0] * np.log(1-p_interaction)# sto usando la BEER LAMBERT LAW ma non so se posso usarla per i fotoni # uso lambda della cross section totale
        
        # traduco le coordinate sferuche in coordinate cartesiane (x,y,z)
        pos = from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0) # posizione dell'interazione
        
        # controllo di stare dentro il rettangolo
        if(  (pos >= c.pos_min ).all()  & (pos <= c.pos_max).all() ):
 
            # segno quante interazioni accadono in un subrettangolo
            if( (pos >= c.pos_min ).all() & (pos <= c.sub_rect).all() ): args[3] = args[3] + 1
                
            # vedo se il protone interagisce con un carbonio o con un protone
            p_atom = random()
            if (p_atom <= p[0]): # allora interagisce con il carbonio
                # vedo se fa urto elastico o inelastico
                p_type = random() 
            
                if(p_type <= p[4] ): # ovvero se sono all'interno di uno scattering ELASTICO con il carbonio
                    A = 12
                    args[0].write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\n')
                    
                    # calcolo la perdita di energia e di conseguenza il theta_sc (da capire)

                    theta_scat, E = scattering_angle(E, A)
                    if E < c.MeV1: break # QUESTA SOGLIA VA ABBASSATA (?)
                    j = j+1 # step successivo

                
                else: # ovvero se il neutrone fa scattering inelastico
                    args[0].write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\n')
                    args[1].write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
                    j = j+1 # step successivo

                    break

            else: # allora interagisce con il protone
                # vedo se fa urto elastico o inelastico
                p_type = random() 
            
                if(p_type <= p[2] ): # ovvero se sono all'interno di uno scattering elastico con il protone
                    A = 1
                    args[0].write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\n')
                    
                    
                    # calcolo la perdita di eienrgia e di conseguenza il theta_sc (da capire)
                    theta_scat, E = scattering_angle(E, A)
                    if E < c.MeV1: break # QUESTA SOGLIA VA ABBASSATA (?)

                    j = j+1 # step successivo

                
                else: # ovvero se il neutrone fa scattering inelastico
                    args[0].write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\n')
                    args[1].write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
                     
                    j = j+1 # step successivo

                    break
                
                
        elif((pos!=pos_source).all()):
            # print('particle out of detector')
            args[1].write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\n')
     
        x0,y0,z0 = pos

        
    return args[3]