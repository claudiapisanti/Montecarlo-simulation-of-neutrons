#import functions as func
from re import sub
import numpy as np

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
    #print(string)
    x = float(string.split(' ')[0])
    y = float(string.split(' ')[1])
    z = float(string.split(' ')[2])
    pos = np.array([x, y, z])

    return pos

    


#########################################################

# open files
w = open('constants.py', 'w')
f = open("physical_characteristics.txt", "r")

# limit value for energy (maximum and minimum value of cs.txt)
Maximum_energy = 98100000
Minimum_energy = 1100000

# inizializzo
E_mono = 0
E_max = 0
E_min = 0
face_prob_cum = 0
sph_radius = 0
pos_source = [0.,0.,0.]





Na = 6.0221409e23 # Avogadro's number
MeV1 = 1000000 # = 1 MeV

# il mio materiale e il polyvinyl toluene
# E = 100 # eV energia (è monoenergetico) # !!!!!!!!!!!!
rho = 1.023 #g/cm^3 densità                 
Mmol = 118. # g/mol massa molare 
n = Na*rho/Mmol # densità di particelle target 


# vedo che energia prendere
line  = f.readline() # line 1 --- type of source
En_type = line.split(' ')[0]
if(En_type == 'MONO'): # per una sorgente monoenergetica mi ottengo le cross section dal file.txt
    line = f.readline() # in eV # line 2 --- energy of source
    E_mono = float(line.split(' ')[0])

    # check if E_mono is inside the proper range
    assert E_mono < Maximum_energy, "Energy is larger than maximum energy given in cs.txt. Correci physical_caracteristics.txt"
    assert E_mono > Minimum_energy, "Energy is smaller than minimun enrgy given in cs.txt. Correci physical_caracteristics.txt"

    



elif(En_type == 'UNIF'):
    line = f.readline() # in eV  # line 2 ---energy of source
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




# GENERO UN RETTANGOLO (lo scintillatore) E UN SUBRETTANGOLO (detectabile) -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
pos_max = get_pos(f) # read line 3 --- geometrical boundaries of the scintillatro
pos_min = np.array([0,0,0])

print("\nLo scintillatore ha dimensioni ", 
    pos_max[0], 'x', pos_max[1],'x', pos_max[2], 'cm')

# seleziono una zona di interesse dove contare il numero di eventi
sub_rect = get_pos(f) # read line 4 --- geometrical boundaries for sub-rectangle of interest
  

if( (sub_rect >= pos_max ).any()):
        print ('errore: il la sottoregione eccede la regione. Aggiusta il tuo file')
        quit() # chiudi il programma se i dati di input sono sbagliati
      
print("\nIl la sottoregione ha coordinate (0.0,", sub_rect[0], '; 0.0,', sub_rect[1],';0.0,', sub_rect[2], ')')


# METTI UN WHILE PER VEDERE CONTINUARE A CHIEDERE INPUT SE IL PROGRAMMA FINO A CHE NON VANNO BENE LE CONDIZIONI SE IL PROGRAMMA 
# SIA PER SUBRECT SIA PER LA SORGENTE

###########################################################
###########################################################
# DEFINISCO IL TIPO DI SORGENTE (puntiforme/estesa/sferica)
type_source = f.readline().split(' ')[0] # line 5 --- type of source : EST, PUNT, SPH



# POSIZIONE SORGENTE

if (type_source == 'PUNT'):
    print('puntiforme')
    # posizione della sorgente
    pos_source = get_pos(f) # line 6 --- position of point source
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
    face_prob_line = f.readline() # line 6 --- face probbilities

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

k = 0
w.write("import numpy as np\n\n")

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
#w.write(f'face = {face}\n')
w.write(f'sub_rect = np.array([{sub_rect[0]},{sub_rect[1]},{sub_rect[2]}])\n')
#w.write(f'k = {k}\n')
w.write(f'MeV1 = {MeV1}')


w.close()
f.close()

