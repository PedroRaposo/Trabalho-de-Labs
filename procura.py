# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 19:40:28 2015

@author: Emanuel
"""

import os
from Bio import SeqIO
from Bio import Entrez
from Bio import SwissProt
import re
Entrez.email = "emanuel_queiroga1@hotmail.com"     # Always tell NCBI who you are
genoma = "genoma.gb"
ls = []    
file = open(genoma)
for line in file: # abrir o ficheiro GenBank
    ls.append(line)
#print(ls)