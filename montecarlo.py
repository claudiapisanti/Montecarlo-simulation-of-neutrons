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
seed(1)



# CONSTANTS -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# read info from file
f = open("physical_characteristics.txt", "r")


Na = 6.0221409e23 # Avogadro's number

N = int(input('insert number of particles:'))
# provo un materiale di acqua
E = 100 # eV energia (è monoenergetico)
rho = 1 #g/cm^3 densità
Mmol = 18.01528 # g/mol massa molare


# CROSS SECTIONS, densità di molecole, etc -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# costruisco la mia tabella della cs e mi calcolo 
cs = {'hydrogen': [20.4252, 20.4305, 0.0053], 
      'oxygen':[3.79371, 3.793714, 0.0 ],
      }
cs['water'] = [2*cs['hydrogen'][0] + cs['oxygen'][0],
               2*cs['hydrogen'][1] + cs['oxygen'][1],
               2*cs['hydrogen'][2] + cs['oxygen'][2] ]

cs_df = pd.DataFrame(cs, index = ['elastic', 'total', 'inelastic'])

cs_df = np.dot(cs_df, 10e-24) # converto la cs da barn a cm^2

# cs_wat = 2*c

n = Na*rho/Mmol # densità di particelle target
l = cs_df * n

p_elastic = cs_df[0,2] / cs_df[1,2]
p_inelastic = (cs_df[1,2] - cs_df[0,2]) / cs_df[1,2]

# GENERO UN RETTANGOLO (lo scintillatore) E UN SUBRETTANGOLO (detectabile) -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
print("""
      MONTECARLO METHOD.
      Inserisci coordinate x,y,z in centimetri del volume da considerare: \n
       """)

pos_max = func.get_pos(f) # this is the first line of the code
pos_min = np.array([0,0,0])


print("\nIl rettangolo ha dimensioni ", 
    pos_max[0], 'x', pos_max[1],'x', pos_max[2], 'cm')


print("""
      Ora inserisci le dimensioni della sottoregione dove misurare quante
      interazioni accadono.: \n
      """)
  

sub_rect = func.get_pos(f) # questa è la seconda riga di codice
  

if( (sub_rect >= pos_max ).any()):
        print ('errore: il la sottoregione eccede la regione. Aggiusta il tuo file')
        quit() # chiudi il programma se i dati di input sono sbagliati
    
       
print("\nIl la sottoregione ha coordinate (0.0,", sub_rect[0], '; 0.0,', sub_rect[1],';0.0,', sub_rect[2], ')')


# METTI UN WHILE PER VEDERE CONTINUARE A CHIEDERE INPUT SE IL PROGRAMMA FINO A CHE NON VANNO BENE LE CONDIZIONI SE IL PROGRAMMA 
# SIA PER SUBRECT SIA PER LA SORGENTE

###########################################################
###########################################################
# DEFINISCO IL TIPO DI SORGENTE (puntiforme/estesa/sferica)
type_source = f.readline().split(' ')[0]
print(type_source)





# POSIZIONE SORGENTE


if (type_source == 'PUNT'):
    print('puntiforme')
    print("""
         Inserisci la posizione della sorgente: 
            """)

    pos_source = func.get_pos(f) # questa è la terza riga di codice
    print(pos_source)
    face = 0

    while((pos_source >= pos_max ).any()):
        print ('errore: il la sorgente è fuori dallo scintillatore. Correggi il file.txt')
        quit()


    print("La sorgente è situata in (", pos_source[0],',', pos_source[1], ',', pos_source[2], ')')
elif(type_source == 'EST'):
    print('esteso')
    # SCELTA DELLA FACCIA

    face_prob = np.zeros(6)
    face_prob_cum = np.zeros(6)
    face_prob_line = f.readline()

    for i in range(6):
       face_prob[i] = float(face_prob_line.split(' ')[i])


    for i in range(6):
        face_prob_cum[i] = sum(face_prob[:(i+1)]) / sum(face_prob) # cumulativa

    print("------------> ", face_prob, sum(face_prob), face_prob_cum)
elif(type_source == 'SPH'):
    print('spherical')
    sph_radius = float(f.readline().split(' ')[0]) # (cm) # non mi metto esattamente sul bordo perché potrei avere problemi al bordo /2 perché metto il centro della sfera al centro dello scintillatore
    print('-------------------> '+ str(sph_radius))
    face = 0 # otw mi da errore




##############################################################
##############################################################
##############################################################
# MONTECARLO

xs = []
ys = []
zs = []

step = open("step.txt", "w")
event = open("event.txt", "w")
step.write('evento\tstep\tx\ty\tz\tface\tinteraction\n')
event.write('evento\tx\ty\tz\tface\n')

k = 0
p_type = 0



for i in range (0,N): # eventi
    print('evento = ', i)
    if(type_source == 'EST'): # ottengo una posizione iniziale
        face = func.face_func(face_prob_cum)
        pos_source = func.source(face,pos_max[0],pos_max[1],pos_max[2])
    elif(type_source == 'SPH'): # ottengo una posizione finale
        phi_source = func.random_rescale(np.pi)
        theta_source = func.random_rescale(2*np.pi)
        pos_source = func.from_sph_coord_to_xyz(sph_radius,phi_source,theta_source,pos_max[0]/2.,pos_max[1]/2.,pos_max[2]/2.) # centrato al centro del sistema
        
        
    
    j = 0
    p_type = 0.
    x0 = pos_source[0]
    y0 = pos_source[1]
    z0 = pos_source[2]
    pos = [0.,0.,0.] # ??
    
    step.write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\tsource\n') # sorgente

    while( (pos >= pos_min ).all()  & (pos <= pos_max).all()):
       
        # mi conviene lavorare in coordinate polari
        phi = random() * np.pi # phi è compreso tra 0 e pi 
        theta = random() * 2 * np.pi   # thetha è compreso tra 0 e 2pi
        
        p_interaction = random() # probabilità di interazione con cui calcolare lo spazio
        r = - l[1,2] * np.log(1-p_interaction)# sto usando la BEER LAMBERT LAW ma non so se posso usarla per i fotoni
        
        # traduco le coordinate sferuche in coordinate cartesiane (x,y,z)
        pos = func.from_sph_coord_to_xyz(r,phi,theta,x0,y0,z0)
        
        # controllo di stare dentro il rettangolo
        if(  (pos >= pos_min ).all()  & (pos <= pos_max).all() ):
            
            # print( (pos >= pos_min ).all()  & (pos <= pos_max).all())
            
            
            # segno quante interazioni accadono in un subrettangolo
            if( (pos >= pos_min ).all() & (pos <= sub_rect).all() ): k = k+1
                
            p_type = random() # tipo di interazione
            
            if(p_type <= p_elastic ): 
                # print('evento=',i,'step=', j)
                # print(pos)
                
                step.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\telastic\n')
                
                j = j+1 # step successivo
                
            else:
                # print('particle got absorbed') 
                step.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\tinelastic\n')
                event.write(f'{i}\t{j}\t{pos[0]}\t{pos[1]}\t{pos[2]}\t{face}\n')
                
        elif((pos!=pos_source).all()):
            # print('particle out of detector')
            event.write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\t{face}\n')
     
        x0,y0,z0 = pos

        
    
    
    
step.close()  
event.close()  

print("nel rettangolo(0.0,", sub_rect[0], '; 0.0,', sub_rect[1],';0.0,', sub_rect[2], '), sono avvenute k = ', k, 'interazioni' )

##########################################
##########################################
##########################################
# IMPORTO LE TABELLE

step = pd.read_csv('step.txt', sep = '\t', index_col = False)
event =  pd.read_csv('step.txt', sep = '\t', index_col = False)

# 

# from mpl_toolkits import mplot3d

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

plt.xlabel('x')
plt.ylabel('y')
plt.savefig('tot_eventi.png')
plt.close()
    
#%%

plt.title("Primi 30 eventi ")
for i in range (0,30):
    
    x = step['x'][step['evento'] == i]
    y = step['y'][step['evento'] == i]
    plt.plot(x,y, alpha = 0.7)
    
plt.savefig('20_eventi.png')
plt.xlabel('x')
plt.ylabel('y')
plt.close()


# +380633659226 Anna Grizan a.grizan@kodland.team 
          
          
