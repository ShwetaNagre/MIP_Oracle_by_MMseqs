# -*- coding: utf-8 -*-
"""
@author: Shweta Nagre
"""
import pandas as pd
from Bio.Seq import Seq

df = pd.read_excel (r"Passable MIPs.xlsx")
ligarm = df["Ligation Arm"].tolist()
extarm = df["Extension Arm"].tolist()
target = df["Target region"].tolist()
accid = df["Def Line"].tolist()
org = df["Organism"].tolist()
fname="whole_region"
tag2=[]
cnt2=0
#number the sequences and add the organisms to the definition line
# Process each definition line and organism
for k in accid:
    k = k.rstrip("\n")
    
    # Remove the space before the organism name within square brackets and replace spaces inside the square brackets with underscores
    if ' [' in k:
        k = k.replace(' [', '[')
        
    # Replace spaces within the square brackets with underscores
    start = k.find('[')
    end = k.find(']')
    if start != -1 and end != -1:
        k = k[:start+1] + k[start+1:end].replace(" ", "_") + k[end:]
    
    # Replace spaces with underscores in organism names
    cleaned_org = org[cnt2].replace(" ", "_")
    
    # Create the new tag with the correct format
    tag2.append(k + "_whole_region_" + str(cnt2) + "|" + cleaned_org)
    cnt2 += 1
ofile = open(fname+".txt", "w")
cnt=0
#create a FASTA file with new definition lines and the arm1+target+arm2 sequence
for i in range(len(accid)):
    arm1put=Seq(target[i])
    larm=Seq(ligarm[i]).complement()
    earm=Seq(extarm[i]).reverse_complement()
    cnt+=1
    ofile.write(str(tag2[i]) + "\n"+str(larm)+str(arm1put)+str(earm)+"\n")

ofile.close()
