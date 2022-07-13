# IMPORT LIBRARIES
from random import random, seed
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
import os
import glob

# generate from physical characteristics a file of constants. This will reduce 
os.system('python3 make_constants.py')

# my functions
import constants as c # I want to first generate the file, and then use it
import functions as func


# fix seed for random
seed(1)

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# read info from file and get number of particles form bash
cs_table = pd.read_csv('./cross_sections/cs.txt', sep = '\t')

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# MONTECARLO


# # WITHOUT MULTIPROCESSING
# start = time.time()
# for i in range (0,c.n_processes): # eventi
#     func.event_func(i, cs_table)

# end = time.time()

# print("TIME no MP = ", end - start)
# print('\n')



# WITH MULTIPROCESSING
if __name__ == "__main__":

    start = time.time()

    # Start processes in asyncronous way
    with mp.Pool(processes=c.n_processes) as pool:

        for i in range(c.n_processes):
            pool.apply_async(func.event_func, args = (i, cs_table))
        
        pool.close()
        pool.join()
        

    end = time.time()
    print("TIME MP= ", end - start)


# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# FILE MERGING
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

# merge al temporary files in a single file
func.merge_tmp_tables(step_name, tmp_step_list)
func.merge_tmp_tables(event_name, tmp_event_list)

# delete temporary file
for tmp_step in tmp_step_list:
    os.remove(tmp_step)

for tmp_event in tmp_event_list:
    os.remove(tmp_event)


