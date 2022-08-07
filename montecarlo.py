# IMPORT LIBRARIES
from random import seed
import pandas as pd
import multiprocessing as mp
import time
import os
import glob
import sys

# my function
import functions as func
import make_constants as mkc


# check argparse
if len(sys.argv) < 2:
    print ("Using macro file 'example.txt'")
    macro_path = 'example.txt'
elif len(sys.argv)==2:
    macro_path = sys.argv[1]
    print(f"Using file '{sys.argv[1]}'")


macro = open(macro_path, 'r')
data = mkc.make_dictionary(macro)

n_processes = data['n_processes'] # get number of processes from dictionary

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# read info from file and get number of particles form bash
cs_table = pd.read_csv('./cross_sections/cs.txt', sep = '\t')

# -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
# MONTECARLO


# # # WITHOUT MULTIPROCESSING
# start = time.time()
# for i in range (0,n_processes): # eventi
#     func.event_func(i, cs_table, data)

# end = time.time()

# print("TIME no MP = ", end - start)
# print('\n')



# WITH MULTIPROCESSING
if __name__ == "__main__":

    start = time.time()

    # Start processes in asyncronous way
    with mp.Pool(processes=n_processes) as pool:

        for i in range(n_processes):
            pool.apply_async(func.event_func, args = (i, cs_table, data))
        
        pool.close()
        pool.join()
        

    end = time.time()
    print("Simulation finished in = ", round(end - start, 3), f"seconds with {n_processes} processes")


    # -:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
    # FILE MERGING
    # merge all temporary files in a singel file
    step_name = "step.txt"
    event_name = "event.txt"

    # write indentation
    step = open(step_name, 'w')
    event = open(event_name, 'w')
    # write column names (only the first time)
    step.write('event,step,x,y,z,interaction,Energy\n')
    event.write('event,last_step,x,y,z\n')
    step.close()
    event.close()

    # select all temporary files
    tmp_step_list = glob.glob("tmp_step*.txt")
    tmp_event_list = glob.glob("tmp_event*.txt")

    # merge al temporary files in a single file
    func.merge_tmp_tables(step_name, tmp_step_list)
    func.merge_tmp_tables(event_name, tmp_event_list)

    # delete temporary file
    for tmp_step in tmp_step_list:
        os.remove(tmp_step)

    for tmp_event in tmp_event_list:
        os.remove(tmp_event)


