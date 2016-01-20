# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 12:44:22 2015

@author: Emanuel
"""

#from Bio import Entrez
#from Bio import SeqIO
#Entrez.email = "emanuel_queiroga1@hotmail.com" 

#Gives the names of NCBI's databases
#handle = Entrez.einfo()
#result = handle.read()
#print(result)

#handle = Entrez.efetch(db="nucleotide", term="NC_000919", rettype="gb", retmode="text")
#record = SeqIO.read(handle, "genbank")
#handle.close()
#print(record)

import os
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "riotforce@sapo.pt"     # Always tell NCBI who you are
filename = "g6full2.gbk"
#if not os.path.isfile(filename):
#    # Downloading...
#    net_handle = Entrez.efetch(db="nucleotide",id="15638995",rettype="gb", retmode="file")
#    out_handle = open(filename, "w")
#    print(net_handle.read())
#    out_handle.write(net_handle.read())
#    out_handle.close()
#    net_handle.close()
#    print("Saved")
#
#print("Parsing...")
#record = SeqIO.read(filename, "genbank")
#print(record)


net_handle = Entrez.efetch(db="nucleotide",id="15638995",rettype="gb", property='complete genome', retmode="text")
line=""
line = net_handle.readline()
while (line ):
    print(line)
    line = net_handle.readline()
net_handle.close()