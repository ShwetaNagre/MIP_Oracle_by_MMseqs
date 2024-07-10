# -*- coding: utf-8 -*-
"""
@author: Shweta Nagre
"""
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import re
import pandas as pd
import os
import sys
import configparser
fname = sys.argv[1]

#read in the length of the target region and the length of the extension and ligation arms
config = configparser.ConfigParser()
config.read_file(open(r'config.txt'))
rep = config.get('My Section', 'Number of repeats',fallback='No Value Entered')
armlengthel = config.get('My Section', 'length of MIP ext and lig arm',fallback='No Value Entered')
targetlen = config.get('My Section', 'length of target region',fallback='No Value Entered')

alllen=int(targetlen)+int(armlengthel)+int(armlengthel)


#check definiton line for organism from excel file Organism Dictionary
#following this if organism not found ask for user input to identify the Organism
#full scientific name of the organism should be entered. Ex. Streptomyces lavendulae

df=pd.read_excel('Organism Dictionary.xlsx', index_col=None) 
org=df['Organism'].tolist()
accid=[]
seqs=[]

#Read in the sequences and definition lines from the FASTA file
with open(fname+".fasta") as fh:
    for line in fh:
        if line.startswith(">"):
            accid.append(line)
fasta_sequences = SeqIO.parse(open(fname+".fasta"),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    seqs.append(sequence)
bad_chars = [';', ':', '!', "*", '_']
cnt=0
cntno=0
org_notfound=[]
organism=[]
#Search for the organism against the dictionary file
for i in accid:
    for a in bad_chars:
        i=i.replace(a," ")
    #check for any special characters and remove them
    if any(o in i for o in org):
        matched = [o for o in org if o in i]
        organism.append(matched[0])
        cnt+=1
    else:
        cntno+=1
        org_notfound.append(i)
if cntno != 0:
    for i in org_notfound:
        print("An organism could not be detected in the following definition line:")
        print(i)
        print("Please enter the full scientific name of the organism:")
        org1=input()
        organism.append(org1)


arm1=[]
arm2=[]
v4seq=[]
temp=[]
gccontent=[]
temp1=[]
gccontent1=[]
stretch=[]
stretch1=[]

# v4 region is  576-682 numbering based on the E. coli system of nomenclature (Brosius et al., 1978) 

target=[]
mseq=[]
revarm1=[]
revarm2=[]
revtarget=[]
accid2=[]
organism2=[]
extstart=[]
extend=[]
ligationstart=[]
ligationend=[]

#generate arms for both the forward and reverse strand
for s in seqs:
    indi=seqs.index(s)
    seqmain=s
    acid=accid[indi]
    org2=organism[indi]
    out =[]
    for i in range(len(s)):
        for j in range(i + 1, len(s) + 1):
           if (len(s[i:j]) == alllen): 
                out.append(s[i:j])
                extstart.append(i)
                extend.append(i+int(armlengthel))
                ligationstart.append(j-int(armlengthel))
                ligationend.append(j)
    lenarm=int(armlengthel)
    for g in out:
        #EDIT HERE
        organism2.append(org2)
        accid2.append(acid)
        mseq.append(seqmain)
        a1=Seq(g[:lenarm])
        a2=Seq(g[-lenarm:])
        #EDIT HERE
        arm1.append(a1.complement())
        arm2.append(a2.reverse_complement())
        target.append(g[lenarm:-lenarm])
        srev2=a2.complement()
        revarm1.append(a1)
        revarm2.append(srev2.reverse_complement())
        
temprev=[]
gccontentrev=[]
stretch1rev=[]
temprev1=[]
gccontentrev1=[]
stretch1rev1=[]

#get temperature, gc content, and check for polyNs
for a in arm1:
    myseq = a
    inn=arm1.index(a)
    temp.append('%0.2f' % mt.Tm_NN(myseq, nn_table=mt.DNA_NN2))
    gccontent.append(GC(myseq))
    a=str(a)
    if re.search(r"G{rep,}", a) or re.search(r"A{rep,}", a) or re.search(r"C{rep,}", a) or re.search(r"T{rep,}", a):
        stretch.append('yes')
    else:
        stretch.append('no')
for b in arm2:    
    myseq = b
    temp1.append('%0.2f' % mt.Tm_NN(myseq, nn_table=mt.DNA_NN2))
    gccontent1.append(GC(myseq))
    b=str(b)
    if re.search(r"G{rep,}", b) or re.search(r"A{rep,}", b) or re.search(r"C{rep,}", b) or re.search(r"T{rep,}", b):
        stretch1.append('yes')
    else:
        stretch1.append('no')
for c in revarm1:    
    myseq = c
    temprev1.append('%0.2f' % mt.Tm_NN(myseq, nn_table=mt.DNA_NN2))
    gccontentrev1.append(GC(myseq))
    c=str(c)
    if re.search(r"G{rep,}", c) or re.search(r"A{rep,}", c) or re.search(r"C{rep,}", c) or re.search(r"T{rep,}", c):
        stretch1rev1.append('yes')
    else:
        stretch1rev1.append('no')

for d in revarm2:
    myseq = d
    
    temprev.append('%0.2f' % mt.Tm_NN(myseq, nn_table=mt.DNA_NN2))
    gccontentrev.append(GC(myseq))
    d=str(d)
    if re.search(r"G{rep,}", d) or re.search(r"A{rep,}", d) or re.search(r"C{rep,}", d) or re.search(r"T{rep,}", d):
        stretch1rev.append('yes')
    else:
        stretch1rev.append('no')
for i in target:
    revtarget.append(Seq(i).complement())

#output all MIPs possible for the FASTA file input by the user
df=pd.DataFrame({"Def Line":pd.Series(accid2), "Main Sequence":pd.Series(mseq),"Target region":pd.Series(target),"Ligation Arm":pd.Series(arm1),"TM 1":pd.Series(temp),"GC content 1":pd.Series(gccontent),"Continuous stretch":pd.Series(stretch),"Extension Arm":pd.Series(arm2),"TM 2":pd.Series(temp1),"GC content 2":pd.Series(gccontent1),"Continuous stretch 2":pd.Series(stretch1),"Reverse Strand Target region":pd.Series(revtarget),"Ligation Arm(Rev)":pd.Series(revarm1),"TM Rev":pd.Series(temprev1),"GC content Rev":pd.Series(gccontentrev1),"Continuous stretch Rev":pd.Series(stretch1rev1),"Extension Arm Rev":pd.Series(revarm2),"TM 2 Rev":pd.Series(temprev),"GC content 2 Rev":pd.Series(gccontentrev),"Continuous stretch 2 Rev":pd.Series(stretch1rev),"Extension Start":pd.Series(extstart),"Extension End":pd.Series(extend),"Ligation Start":pd.Series(ligationstart), "Ligation End":pd.Series(ligationend),"Organism":pd.Series(organism2)})

df.to_excel(fname+"_MIPs.xlsx", index=False)
