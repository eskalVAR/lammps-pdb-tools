#!/usr/bin/env python
# coding: utf-8
# Contact map plotter using PDB files created by lammps2pdb.pl
# Copyright (c) 2022, Eryk Skalinski
# eskalinski@protonmail.com

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math 
import altair as alt
import csv
import gc


# In[2]:


bad_words = ['REMARK', 'END', 'CRYST1']

count=1
pdb_dict_list = []
with open('./2beg_assoc_ngo_trj.pdb') as pdb_file:
    total = len(open('./2beg_assoc_ngo_trj.pdb').readlines())
    skip = total*0.005
    for i,line in enumerate(pdb_file):
        if not any(bad_word in line for bad_word in bad_words):
            pdb_dict={
                "header":line[0:6].strip(),
                "id":line[6:12].strip(),
                "Atom_id":line[12:16].strip(),
                "res_id":line[17:20].strip(),
                "res_number":line[22:26].strip(),
                "x":line[30:38].strip(),
                "y":line[38:46].strip(),
                "z":line[46:54].strip(),
                "q_num1":line[54:60].strip(),
                "q_num2":line[60:66].strip(),
                "protein":line[66:].strip()  
            }
            count = count + 1
            if (count > skip) :
                count = 0
                print(i*100/total, end ='\r')
            
            pdb_dict_list.append(pdb_dict)
print("")
print("Data loaded.", end='\n')


# In[3]:


pdb = pd.DataFrame.from_dict(pdb_dict_list)
gc.collect()
dtype0= {
    "header":"str",
    "id":"uint16",
    "Atom_id":"str",
    "res_id":"str",
    "res_number":"uint16",
    "x":"float32",
    "y":"float32",
    "z":"float32",
    "q_num1":"float16",
    "q_num2":"float16",
    "protein":"str" 
}
print("Assigning types.")
pdb= pdb.astype(dtype0)
print("Pre-proc finished.")


# In[ ]:


#print(pdb)


# 
# 
# pdb = pd.read_fwf(
#     "./2beg_association_model_no_go/results/20A_500_spun/2beg_refined.pdb",
#     header=None,
#     infer_nrows=1500,
# )
# col_names = [
#     "header",
#     "id",
#     "Atom_id",
#     "res_id",
#     "res_number",
#     "x",
#     "y",
#     "z",
#     "q_num1",
#     "q_num2",
#     "protein",
# ]
# pdb.columns = col_names
# 

# In[ ]:


res_num = int(pdb.iloc[::-1].iloc[0][4])
print('Residues : {:d}'.format(res_num))
atom_num = pdb.iloc[::-1].iloc[0][1]
print('Number of atoms : {:d}'.format(atom_num))
snapshot_num = pdb[pdb.id == atom_num].shape[0]
print('Number of snapshots : {:d}'.format(snapshot_num))
overall = pd.DataFrame()

print("Starting main loop.")


# In[ ]:


count=0
indy=0
step_counter = 1
skip=snapshot_num*0.005
out_dict_list=[]
plt.figure(figsize=(12,8))
for snapshot in range(snapshot_num): # Snap_shot_number
    out_dict_snapshot_list=[]
    raw_pdb = pdb[(atom_num*(step_counter-1)):(atom_num*step_counter)]
    CA_list = raw_pdb.query("Atom_id == 'CA'")
    CA = list(zip(CA_list["res_number"].values, CA_list["x"].values,CA_list["y"].values,CA_list["z"].values))
    for res_num_i, x_i, y_i, z_i  in CA:
        for res_num_j, x_j, y_j, z_j in CA:
            calc_dist = math.dist([x_i, y_i, z_i], [x_j, y_j, z_j])
            if calc_dist < 6:            
                out_dict={
                    "index": indy,
                    "model" : step_counter,
                    "i" : res_num_i,
                    "j" : res_num_j,
                    "distance" : calc_dist
                }
                out_dict_snapshot_list.append(out_dict)
                out_dict_list.append(out_dict)            
                indy=indy+1
                #Data = Data.append({'res_num_i': res_num_i, 'res_num_j': res_num_j ,
                #'i': [x_i, y_i, z_i], 'j': [x_j, y_j, z_j], 'dist': calc_dist},
                #ignore_index=True)   
             #   overall = overall.append(
              #      {'model' : step_counter, 'i': res_num_i, 'j': res_num_j,'distance': calc_dist},
               #      ignore_index=True)
               # overall.to_csv('plotting_material.csv')
    step_counter += 1
    count += 1
    if (count > skip) :
        count = 0
        x = [element["i"] for element in out_dict_snapshot_list]
        y = [element["j"] for element in out_dict_snapshot_list]
        plt.scatter(x,y, marker="^", cmap="viridis", color="k")
        plt.xlabel("Residue index Number", size = 15)
        plt.ylabel("Residue index Number",  size = 15)
        plt.gca().invert_yaxis()
        id_name = str(step_counter)
        id_name=id_name.zfill(5)
        plt.title("C"+"\u03B1"+"-"+"C"+"\u03B1"+" Protein Contact Map" + str(step_counter), size = 25)
        name =  "p_contact_map" + id_name
        plt.savefig(name)
        plt.cla()
        print(step_counter*100/snapshot_num, end ='\r') 
    
print("")
print("Finished processing.")
print("Writing output...")

with open('out.csv', 'w', encoding='utf8', newline='') as output_file:
    fc = csv.DictWriter(output_file, 
        fieldnames=out_dict_list[0].keys(),

                       )
    fc.writeheader()
    fc.writerows(out_dict_list)
print("Done.")
                    


# In[ ]:





# In[ ]:




