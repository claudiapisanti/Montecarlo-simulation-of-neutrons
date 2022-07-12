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
# read info from file and get number of particles form bash
cs_table = pd.read_csv('./cross_sections/cs.txt', sep = '\t')

N = int(input('insert number of particles:'))

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
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
print(len(tmp_event_list))


func.merge_tmp_tables(step_name, tmp_step_list)
func.merge_tmp_tables(event_name, tmp_event_list)

# delete temporary file
for tmp_step in tmp_step_list:
    os.remove(tmp_step)


for tmp_event in tmp_event_list:
    os.remove(tmp_event)


