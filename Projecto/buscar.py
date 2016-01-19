import os
from Bio import SeqIO
from Bio import Entrez
from Bio import SwissProt
import re
import inspect
import Gene

import re

#Get data from NCBI

Entrez.email = "mafaldasap@gmail.com"     # Always tell NCBI who you are
handle4 = Entrez.efetch(db="nucleotide", id="15638995", rettype="gb", strand=1, seq_start=639851, 
                       seq_stop=765100)

handle5 = Entrez.efetch(db="nucleotide", id="15638995", rettype="fasta", strand=1, seq_start=639851, 
                       seq_stop=765100)

#---------------------------------------------------------------------------

#print(str(recordF.seq)[122609:124403])
filename2 = "g6.gb"
fastaFile = "g6.fasta"
save_file2 = open(filename2, "w")
save_file2.write(handle4.read())
save_file2.close()
save_fileF = open(fastaFile, "w")
save_fileF.write(handle5.read())
save_fileF.close()
handle4.close()
handle5.close()

#------------------------------------------------------------
#record4 = SeqIO.read(handle4, "genbank")
#recordF = SeqIO.read(handle5, "fasta")
#-------------------------------------------------------------

record4 = SeqIO.read(filename2, "genbank")
recordF = SeqIO.read(fastaFile, "fasta")
#Save information in variables
fastaSeq= str(recordF.seq)

#Qualifiers expressions
locus_tag="locus_tag"
old_locus_tag="old_locus_tag"
db_xref="db_xref"
product="product"
#only for trna
anticodon="anticodon"
#only for CDS
accession="protein_id"
translation="translation"
ec="EC_number"
tc="TC_number"

#-----------------------------------------------------

#print (record4.seq)  # recupera a nossa parte do genoma a estudar
#print(inspect.getdoc(record4))

cdsList=[]
trnaList=[]

features = record4.features
for feature in features:
    
    f= str(feature.location) 
    if(str(feature.type)!="gene" and str(feature.type)!="source"):        
        #l=re.search("\[[0-9]+:[0-9]+\]",f)
        #print("Location: %s" % l.group())
        locations=re.findall("[0-9]+",f)
        l1=int(locations[0])+639851
        l2=int(locations[1])+639850
        loc="["+str(l1)+":"+str(l2)+"]"
#        print("Type: %s" % feature.type)
#        print("Strand: %s" % feature.strand)
        qualifiers= feature.qualifiers
        #print(feature.qualifiers)
#        print("Locus_tag :%s" % qualifiers.get(locus_tag))
#        print("db_xref :%s" % qualifiers.get(db_xref))
#        print("product :%s" % qualifiers.get(product))
        if(str(feature.type)=="tRNA"):
#            print("anticodon :%s" % qualifiers.get(anticodon))
            trna= Gene.MyTRNA(str(feature.type),str(feature.strand),loc,str(qualifiers.get(locus_tag)),str(qualifiers.get(old_locus_tag)), str(qualifiers.get(db_xref)),str(qualifiers.get(product)),str(qualifiers.get(anticodon)),fastaSeq[int(locations[0]):int(locations[1])] )
            trnaList.append(trna)            
            #print(trna)
            #print("----------------------------------------------------------------------")
        else:
            cDS= Gene.MyCDS(str(feature.type),str(feature.strand),loc,str(qualifiers.get(locus_tag)),str(qualifiers.get(old_locus_tag)), str(qualifiers.get(db_xref)),str(qualifiers.get(product)),str(qualifiers.get(accession)),str(qualifiers.get(translation)),fastaSeq[int(locations[0]):int(locations[1])],str(qualifiers.get(ec)),str(qualifiers.get(tc)) )
            cdsList.append(cDS)
            #cDS.blast()
            
            #cDS.uniprotSearch()  
#            print("-------------------------------------------")
            #print(cDS)
#            print("---------------------------------------------")
            #print("----------------------------------------------------------------------")
#            print("accession :%s" % qualifiers.get(accession))
#            print("translation :%s" % qualifiers.get(translation))
#        print(fastaSeq[int(locations[0]):int(locations[1])])
        
    
    
    
    #print("Qualifiers: %s" % feature.qualifiers)
    
#cdsList[0].uniprotSearch()
cdsList[0].blast()

#print(inspect.getdoc(feature))

#print(for method in dir(feature) if callable(getattr(feature, method)))
#sig= inspect.signature(feature)
#print(sig);
#feature = features[0]