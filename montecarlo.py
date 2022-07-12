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
import multiprocessing as mp
import time
import constants as c
import os
import glob

# fix seed for random
seed(1)

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-

# CONSTANTS -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# read info from file
cs_table = pd.read_csv('./cross_sections/cs.txt', sep = '\t')

N = int(input('insert number of particles:'))

##############################################################
##############################################################
##############################################################
# MONTECARLO
n_processes = mp.cpu_count()
print('Using ',n_processes, ' processes.')


# start = time.time()
# for i in range (0,N): # eventi
#     #k = func.evento_letto_da_const(i,[step,event,cs_table,k])
#     func.evento_letto_da_const(i, cs_table)

# end = time.time()

# print("TIME = ", end - start)



start = time.time()
if __name__ == "__main__":
    # Start processes in asyncronous way
    with mp.Pool(processes=n_processes) as pool:
        for i in range(N):
            pool.apply_async(func.evento_letto_da_const, args = (i, cs_table))


end = time.time()
print("TIME = ", end - start)


# merge all temporary files in a singel file
step_name = "step.txt"
event_name = "event.txt"

# write indentation
step = open(step_name, 'w')
event = open(event_name, 'w')
# write column names (only the first time)
step.write('evento\tstep\tx\ty\tz\tface\tinteraction\n')
event.write('evento\tlast_step\tx\ty\tz\tface\n')
step.close()
event.close()

# select all tempoirary files
tmp_step_list = glob.glob("tmp_step*.txt")
tmp_event_list = glob.glob("tmp_event*.txt")


func.merge_tmp_tables(step_name, tmp_step_list)
func.merge_tmp_tables(event_name, tmp_event_list)

# delete temporary file
for tmp_step in tmp_step_list:
    os.remove(tmp_step)


for tmp_event in tmp_event_list:
    os.remove(tmp_event)


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


          
          
