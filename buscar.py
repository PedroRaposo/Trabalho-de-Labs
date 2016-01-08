import os
from Bio import SeqIO
from Bio import Entrez
from Bio import SwissProt
import inspect

import re
Entrez.email = "mafaldasap@gmail.com"     # Always tell NCBI who you are
handle4 = Entrez.efetch(db="nucleotide", id="15638995", rettype="gb", strand=1, seq_start=639851, 
                       seq_stop=765100)
record4 = SeqIO.read(handle4, "genbank")
filename2 = "g6.gb"
save_file2 = open(filename2, "w")
save_file2.write(handle4.read())
save_file2.close()
handle4.close()
#print (record4.seq)  # recupera a nossa parte do genoma a estudar

features = record4.features
for feature in features:
#feature = features[0]
    print("Location: %s" % feature.location)
    print("Type: %s" % feature.type)
#sig= inspect.signature(feature)
#print(sig);

    print("Location Operator: %s" % feature.location_operator)
    print("Strand: %s" % feature.strand)
    print("ID: %s" % feature.id)
    print("Ref: %s" % feature.ref)
    print("Ref_db: %s" % feature.ref_db)
    print("Qualifiers: %s" % feature.qualifiers)
    print("----------------------------------------------------------------------")

#print(inspect.getdoc(feature))

#print(for method in dir(feature) if callable(getattr(feature, method)))