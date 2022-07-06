"""
AGGIUNGI: 
sorgente sferica isotropa
muoni cosmici a 0slm
lettura da file
"""

# IMPORT LIBRARIES
from random import random, seed
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functions as func
#from multiprocessing import Pool
import time
#from functools import partial
import constants as c
#import threading
seed(1)

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# E_mono = 0
# E_max = 0
# E_min = 0
# face_prob_cum = 0
# sph_radius = 0



# CONSTANTS -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# read info from file
cs_table = pd.read_csv('./cross_sections/cs.txt', sep = '\t')


# f = open("physical_characteristics.txt", "r")



# Na = 6.0221409e23 # Avogadro's number
# MeV1 = 1000000

N = int(input('insert number of particles:'))
# # il mio materiale e il polyvinyl toluene
# # E = 100 # eV energia (è monoenergetico) # !!!!!!!!!!!!
# rho = 1.023 #g/cm^3 densità                 
# Mmol = 118. # g/mol massa molare 
# n = Na*rho/Mmol # densità di particelle target 


# # vedo che energia prendere
# line  = f.readline()
# En_type = line.split(' ')[0]
# if(En_type == 'MONO'): # per una sorgente monoenergetica mi ottengo le cross section dal file.txt
#     line = f.readline() # in eV 
#     E_mono = float(line.split(' ')[0])

# elif(En_type == 'UNIF'):
#     line = f.readline() # in eV 
#     E_min = float(line.split(' ')[0])
#     E_max = float(line.split(' ')[1])

# else:
#     print('errore: enrgy type not valid. Correggi il file di input')
#     quit()




# # GENERO UN RETTANGOLO (lo scintillatore) E UN SUBRETTANGOLO (detectabile) -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# pos_max = func.get_pos(f) # this is the first line of the code
# pos_min = np.array([0,0,0])

# print("\nLo scintillatore ha dimensioni ", 
#     pos_max[0], 'x', pos_max[1],'x', pos_max[2], 'cm')

# # seleziono una zona di interesse dove contare il numero di eventi
# sub_rect = func.get_pos(f) # questa è la seconda riga di codice
  

# if( (sub_rect >= pos_max ).any()):
#         print ('errore: il la sottoregione eccede la regione. Aggiusta il tuo file')
#         quit() # chiudi il programma se i dati di input sono sbagliati
      
# print("\nIl la sottoregione ha coordinate (0.0,", sub_rect[0], '; 0.0,', sub_rect[1],';0.0,', sub_rect[2], ')')


# METTI UN WHILE PER VEDERE CONTINUARE A CHIEDERE INPUT SE IL PROGRAMMA FINO A CHE NON VANNO BENE LE CONDIZIONI SE IL PROGRAMMA 
# SIA PER SUBRECT SIA PER LA SORGENTE

###########################################################
###########################################################
# DEFINISCO IL TIPO DI SORGENTE (puntiforme/estesa/sferica)
#type_source = f.readline().split(' ')[0]


# POSIZIONE SORGENTE

# if (type_source == 'PUNT'):
#     print('puntiforme')
#     # posizione della sorgente
#     pos_source = func.get_pos(f) # questa è la terza riga di codice
#     print(pos_source)
#     face = 0

#     while((pos_source >= pos_max ).any()):
#         print ('errore: il la sorgente è fuori dallo scintillatore. Correggi il file.txt')
#         quit()


#     print("La sorgente è situata in (", pos_source[0],',', pos_source[1], ',', pos_source[2], ')')
# elif(type_source == 'EST'):
#     print('esteso')

#     # SCELTA DELLA FACCIA
#     face_prob = np.zeros(6)
#     face_prob_cum = np.zeros(6)
#     face_prob_line = f.readline()

#     for i in range(6):
#        face_prob[i] = float(face_prob_line.split(' ')[i])


#     for i in range(6):
#         face_prob_cum[i] = sum(face_prob[:(i+1)]) / sum(face_prob) # cumulativa

#     print("------------> ", face_prob, sum(face_prob), face_prob_cum)
# elif(type_source == 'SPH'):
#     print('spherical')
#     sph_radius = float(f.readline().split(' ')[0]) # (cm) # non mi metto esattamente sul bordo perché potrei avere problemi al bordo /2 perché metto il centro della sfera al centro dello scintillatore
#     print('-------------------> '+ str(sph_radius))
#     face = 0 # otw mi da errore



##############################################################
##############################################################
##############################################################
# MONTECARLO

xs = []
ys = []
zs = []


# step = open("step.txt", "w")
# event = open("event.txt", "w")
step, event = func.file_init()
step.write('evento\tstep\tx\ty\tz\tface\tinteraction\n')
event.write('evento\tlast_step\tx\ty\tz\tface\n')
# step.close()
# event.close()

k = 0




start = time.time()
for i in range (0,N): # eventi
    # k = func.evento(i, step, event, 
    #     type_source, En_type, 
    #     E_mono, E_max, E_min, n, 
    #     face_prob_cum, pos_max, pos_min, 
    #     sph_radius, 
    #     pos_source, face, sub_rect, k,
    #     MeV1,
    #     cs_table,
    #     )
    
    k = func.evento_letto_da_const(i,[step,event,cs_table,k])

end = time.time()

print("TIME = ", end - start)

# a_pool = Pool()
# start = time.time()
# evento_i = partial(func.evento_letto_da_const, [step,event,cs_table,k] )
# #evento_i(10)
# a_pool.apply_async(evento_i, range(N) )
# #a_pool.map(evento_i, range(N))
# end = time.time()
# print("TIME = ", end - start)


# pool = Pool()
# start = time.time()
# with Pool(processes=N) as pool:
#     for i in range(N):
#         pool.apply_async(func.evento_letto_da_const, args=([step,event,cs_table,k]))
#     pool.close()
#     pool.join()
    
# end = time.time()
# print('TIME = ', end-start)

start = time.time()

end = time.time()
print("TIME = ", end - start)




#print(evento_i)


step.close()  
event.close()  

print("nel rettangolo(0.0,", c.sub_rect[0], '; 0.0,', c.sub_rect[1],';0.0,', c.sub_rect[2], '), sono avvenute k = ', k, 'interazioni' )

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# TO DO: NEW FILE FOR PLOTS!!!

# IMPORTO LE TABELLE

step = pd.read_csv('step.txt', sep = '\t', index_col = False)
event =  pd.read_csv('step.txt', sep = '\t', index_col = False)

# 


fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')

for i in range(6):
    x = event['x'][event['face'] == i]
    y = event['y'][event['face'] == i]
    z = event['z'][event['face'] == i]
    #c = event['face'][event['face'] == i]
    plt.plot(x,y,z, '.')
        #c = c, cmap = 'viridis' )
        
plt.title('Distribuzione eventi')
ax.view_init(25, 45)

ax.axes.set_xlim3d(left=-30.0, right=30.0) 
ax.axes.set_ylim3d(bottom=-30.0, top=30.0) 
ax.axes.set_zlim3d(bottom=-30.0, top=30.0)
ax.set_xlabel('x', fontsize = 12)
ax.set_ylabel('y', fontsize = 12)
ax.set_zlabel('z', fontsize = 12)
plt.savefig('tot_eventi.png')
plt.close()
    
#%%
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
plt.title("Primi 30 eventi ")
ax.view_init(25, 45)

for i in range (0,30):
    
    x = step['x'][step['evento'] == i]
    y = step['y'][step['evento'] == i]
    z = step['z'][step['evento'] == i]
    plt.plot(x,y,z, alpha = 0.7)
    
ax.axes.set_xlim3d(left=0.0, right=6.0) 
ax.axes.set_ylim3d(bottom=0.0, top=6.0) 
ax.axes.set_zlim3d(bottom=0.0, top=6.0)
ax.set_xlabel('x', fontsize = 12)
ax.set_ylabel('y', fontsize = 12)
ax.set_zlabel('z', fontsize = 12)
plt.savefig('30_eventi.png')
plt.show()
#plt.close()


          
          
