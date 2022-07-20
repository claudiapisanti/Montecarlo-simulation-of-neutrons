import pandas as pd
import numpy as np

# apro i file
C_elastic = pd.read_csv("./C_elastic.csv")
H_elastic = pd.read_csv("./H_elastic.csv")
H_total = pd.read_csv("./H_total.csv")
C_total = pd.read_csv("./C_total.csv")

# generate a file of the cross section of the scintillator 
cs_file = open("./cs.txt", "w")
cs_file.write("E\tcs_tot\tcs_el\tcs_inel\tcs_h_tot\tcs_h_el\tcs_h_inel\tcs_c_tot\tcs_c_el\tcs_c_inel\n")

energies = range(1000000, 100000000 ,100000)
for i in range(1,len(energies)):
    # H totale
    if H_total[ (H_total['Incident energy ']>=energies[i -1]) &  (H_total['Incident energy ']<energies[i])].empty != True:
        subtab = H_total[ (H_total['Incident energy ']>=energies[i-1]) &  (H_total['Incident energy ']<energies[i])]
        cs_h_total = np.mean(subtab[' σ(E)'])
    else: continue # go on with the loop (skip this step)

    # C totale
    if C_total[ (C_total['Incident energy ']>=energies[i -1]) &  (C_total['Incident energy ']<energies[i])].empty != True:
        subtab = C_total[ (C_total['Incident energy ']>=energies[i-1]) &  (C_total['Incident energy ']<energies[i])]
        # print(subtab)
        cs_c_total = np.mean(subtab[' σ(E)'])
    else: continue # go on with the loop (skip this step)

    # H elastico
    if H_elastic[ (H_elastic['Incident energy ']>=energies[i -1]) &  (H_elastic['Incident energy ']<energies[i])].empty != True:
        subtab = H_elastic[ (H_elastic['Incident energy ']>=energies[i-1]) &  (H_elastic['Incident energy ']<energies[i])]
        cs_h_elastic = np.mean(subtab[' σ(E)'])
    else: continue # go on with the loop (skip this step)

    # C totale
    if C_elastic[ (C_elastic['Incident energy ']>=energies[i -1]) &  (C_elastic['Incident energy ']<energies[i])].empty != True:
        subtab = C_elastic[ (C_elastic['Incident energy ']>=energies[i-1]) &  (C_elastic['Incident energy ']<energies[i])]
        #print(subtab)
        cs_c_elastic = np.mean(subtab[' σ(E)'])
    else: continue # go on with the loop (skip this step)

    cs_tot = 10 * cs_h_total + 9 * cs_c_total
    cs_el = 10 * cs_h_elastic + 9 * cs_c_elastic
    

    cs_file.write(f'{energies[i]}\t{cs_tot}\t{cs_el}\t{cs_tot - cs_el}\t{cs_h_total}\t{cs_h_elastic}\t{cs_h_total - cs_h_elastic}\t{cs_c_total}\t{cs_c_elastic}\t{cs_c_total - cs_c_elastic}\n')

