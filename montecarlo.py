#%%
"""
rettangolo omogeneo con x,y,z o ingresso a utente
sorgente in 000 o ingresso a utente
sorgente monoenergetica

1)estraggo modulo(lo so già) direzione e verso
    alla fine sta in un altra terna
    
2)
"""
#%% IMPORTO LE LIBRERIE
from random import random, seed
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
seed(1)


Na = 6.0221409e23

N = int(input('insert number of particles:'))
#%% CROSS SECTIONS, densità di molecole, etc
# provo un materiale di acqua
E = 100 # eV energia (è monoenergetico)
rho = 1 #g/cm^3 densità
Mmol = 18.01528 # g/mol massa molare

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

#%% GENERO UN RETTANGOLO (lo scintillatore) E UN SUBRETTANGOLO (detectabile)
print("""
      MONTECARLO METHOD.
      Inserisci coordinate x,y,z in centimetri del volume da considerare: \n
       """)
x_rect = float(input("x: "))
y_rect = float(input("y: "))
z_rect = float(input("z: "))

print("\nIl rettangolo ha dimensioni ", x_rect, 'x', y_rect,'x', z_rect, 'cm')

pos_min = np.array([0,0,0])
pos_max = np.array([x_rect, y_rect, z_rect])

print("""
      Ora inserisci le dimensioni della sottoregione dove misurare quante
      interazioni accadono.: \n
      """)
  
x_subrect = float(input("x: "))
y_subrect = float(input("y: "))
z_subrect = float(input("z: "))

sub_rect = np.array([x_subrect, y_subrect, z_subrect])
  

while( (sub_rect >= pos_max ).any()):
        print ('errore: il la sottoregione eccede la regione. Immetti nuovi valori:')
        
        x_subrect = float(input("x: "))
        y_subrect= float(input("y: "))
        z_subrect = float(input("z: "))
        
        sub_rect = np.array([x_subrect, y_subrect, z_subrect])

    
       
print("\nIl la sottoregione ha coordinate (0.0,", x_subrect, '; 0.0,', y_subrect,';0.0,', z_subrect, ')')


# METTI UN WHILE PER VEDERE CONTINUARE A CHIEDERE INPUT SE IL PROGRAMMA FINO A CHE NON VANNO BENE LE CONDIZIONI SE IL PROGRAMMA 
# SIA PER SUBRECT SIA PER LA SORGENTE

#%% POSIZIONE SORGENTE

print("""
      Inserisci la posizione della sorgente: 
          """)

x_source = float(input("x: "))
y_source = float(input("y: "))
z_source = float(input("z: "))


pos_source = np.array([x_source, y_source, z_source])

while((pos_source >= pos_max ).any()):
        print ('errore: il la sorgente è fuori dallo scintillatore. Immetti nuovi valori:')
        
        x_source = float(input("x: "))
        y_source = float(input("y: "))
        z_source = float(input("z: "))
        
        pos_source = np.array([x_source, y_source, z_source])


print("La sorgente è situata in (", x_source,',', y_source, ',', z_source, ')')

#%% MONTECARLO

xs = []
ys = []
zs = []

step = open("step.txt", "w")
event = open("event.txt", "w")
step.write('evento\tstep\tx\ty\tz\n')
event.write('ciccio')

k = 0
p_type = 0

for i in range (0,N): # eventi
    print('evento = ', i)
    
    j = 0
    p_type = 0.
    x0 = x_source
    y0 = y_source
    z0 = z_source
    pos = [0.,0.,0.]
    step.write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\telastic\n')

    while( (pos >= pos_min ).all()  & (pos <= pos_max).all()):
       
        # mi conviene lavorare in coordinate polari
        phi = random() * np.pi # phi è compreso tra 0 e pi 
        theta = random() * 2 * np.pi   # thetha è compreso tra 0 e 2pi
        
        p_interaction = random() # probabilità di interazione con cui calcolare lo spazio
        r = - l[1,2] * np.log(1-p_interaction)# sto usando la BEER LAMBERT LAW ma non so se posso usarla per i fotoni
        
        # traduco le coordinate sferuche in coordinate cartesiane (x,y,z)
        x1 = x0 + r * np.sin(phi) * np.cos(theta)
        y1 = y0 + r * np.sin(phi) * np.sin(theta)
        z1 = z0 + r * np.cos(phi)
        pos = np.array([x1,y1,z1])
        
        
        # controllo di stare dentro il rettangolo
        if(  (pos >= pos_min ).all()  & (pos <= pos_max).all() ):
            
            # print( (pos >= pos_min ).all()  & (pos <= pos_max).all())
            
            
            # segno quante interazioni accadono in un subrettangolo
            if( (pos >= pos_min ).all() & (pos <= sub_rect).all() ): k = k+1
                
            p_type = random() # tipo di interazione
            
            if(p_type <= p_elastic ): 
                # print('evento=',i,'step=', j)
                # print(pos)
                
                step.write(f'{i}\t{j}\t{x1}\t{y1}\t{z1}\telastic\n')
                
                j = j+1 # step successivo
                
            else:
                # print('particle got absorbed') 
                step.write(f'{i}\t{j}\t{x1}\t{y1}\t{z1}\tinelastic\n')
                event.write(f'{i}\t{j}\t{x1}\t{y1}\t{z1}\n')
                
        elif((pos!=pos_source).all()):
            # print('particle out of detector')
            event.write(f'{i}\t{j}\t{x0}\t{y0}\t{z0}\n')
     
        x0,y0,z0 = x1,y1,z1

        
    
    
    
step.close()  
event.close()  

print("nel rettangolo(0.0,", x_subrect, '; 0.0,', y_subrect,';0.0,', z_subrect, '), sono avvenute k = ', k, 'interazioni' )

#%% IMPORTO LE TABELLE

step = pd.read_csv('step.txt', sep = '\t', index_col = False)
event =  pd.read_csv('step.txt', sep = '\t', index_col = False)

#%% 

# from mpl_toolkits import mplot3d

x = event['x'][:10000]
y = event['y'][:10000]
z = event['z'][:10000]

fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
plt.plot(x,y,z, '.')
plt.title('Distribuzione eventi')
ax.view_init(25, 45)
plt.savefig('tot_eventi.png')
plt.xlabel('x')
plt.ylabel('y')
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
          
          
