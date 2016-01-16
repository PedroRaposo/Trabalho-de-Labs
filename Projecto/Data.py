# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 16:10:14 2016

@author: Emanuel
"""
import os
from Bio import SeqIO
from Bio import Entrez
from Bio import SwissProt
import re
import inspect
import Gene
import re

#Script responsible to get the regist from NCBI trough it's Gi
#its also responssible to parse the file and get all the relevant information
#saving the information in the proper objects defined in the Gene script which 
#are responsible for processing the information

class Data:
    def __init__(self, Name="g6"):
        self.fileName=Name #name of genbank and fasta files
        
        
    #Method responsible with NCBI communication: getting the information and saving it in files     
    def writeFiles(self,database,gi,strandIn,email,seqStart=None,seqStop=None):
        Entrez.email = email     # Its necessary to tell NCBI who we are
        #Using the Entrez import we fetch the genbank and fasta files of the method parameters
        handleGb = Entrez.efetch(db=database, id=gi, rettype="gb", strand=strandIn, seq_start=seqStart, 
                               seq_stop=seqStop)
        
        handleFasta = Entrez.efetch(db=database, id=gi, rettype="fasta", strand=strandIn, seq_start=seqStart, 
                               seq_stop=seqStop)
                               
        self.seqStart= seqStart  # needed to normalize the location index 
        
#       //-------------------------------------------------------------------------------------------------------\\
        #saving in variables the file names with the appropriate extenssion
        gb_filename =self.fileName + ".gb"
        fasta_filename = self.fileName+".fasta"
        #end
        #writing the information from the handles into files
        save_fileGb = open(gb_filename, "w")
        save_fileGb.write(handleGb.read())
        save_fileGb.close()
        save_fileF = open(fasta_filename, "w")
        save_fileF.write(handleFasta.read())
        save_fileF.close()
        #end
        #closing handles
        handleGb.close()
        handleFasta.close()
        #end
        #method end
        
    def getData(self):
        #initialize array variables to store results
        self.cdsList = []
        self.trnaList = []
        #saving in variables the file names with the appropriate extenssion
        gb_filename =self.fileName + ".gb"
        fasta_filename = self.fileName+".fasta"
        #end
        #Initialize record variables from files
        recordGb = SeqIO.read(gb_filename, "genbank")
        recordF = SeqIO.read(fasta_filename, "fasta")
        #end
        
        #Save fasta information in string
        fastaSeq= str(recordF.seq)
        
        #Qualifiers expressions needed for the get value from qualifiers's dictionary
#        //------------------------------------------\\
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
#       //-------------------------------------------\\
        #end
       
        features = recordGb.features # get features from genbank record
        for feature in features:     # go trough all features to retrive information      
            f= str(feature.location) # store the feature location
            if(str(feature.type)!="gene" and str(feature.type)!="source"):  # we dont need to process gene and source      
                locations=re.findall("[0-9]+",f)#getting the numbers from the location with format: [NUMBER:NUMBER]
                #store the locations and normalize them with the total index (index of full sequence)                 
                l1=int(locations[0])+self.seqStart
                l2=int(locations[1])+(self.seqStart-1)
                #end
                loc="["+str(l1)+":"+str(l2)+"]" #store location with desired format
                qualifiers= feature.qualifiers #store qualifiers
                if(str(feature.type)=="tRNA"):
                    #crete MyTRNA object and add it to the trnaList
                    trna= Gene.MyTRNA(str(feature.type),str(feature.strand),loc,str(qualifiers.get(locus_tag)),str(qualifiers.get(old_locus_tag)), str(qualifiers.get(db_xref)),str(qualifiers.get(product)),str(qualifiers.get(anticodon)),fastaSeq[int(locations[0]):int(locations[1])] )
                    self.trnaList.append(trna)
                    #end
                else:
                    # create MyCDS object and add it to cdslist
                    cDS= Gene.MyCDS(str(feature.type),str(feature.strand),loc,str(qualifiers.get(locus_tag)),str(qualifiers.get(old_locus_tag)), str(qualifiers.get(db_xref)),str(qualifiers.get(product)),str(qualifiers.get(accession)),str(qualifiers.get(translation)),fastaSeq[int(locations[0]):int(locations[1])],str(qualifiers.get(ec)),str(qualifiers.get(tc)) )
                    self.cdsList.append(cDS)
                    #end
        #end method
    
    #method responsible for pinting informtion                 
    def printData(self):
        #go though the list and print it's information
        print("CDS and repetition zones:")
        print("\\-------------------------------------------------------------------------------------\\")
        print(" ")
        indice=0
        for cds in self.cdsList:
            print("GENE number ",indice,": ")
            print(cds)
            print()
            indice+=1
        for tr in self.trnaList:
            print("GENE number ",indice,": ")
            print(tr)
            print()
            indice+=1
        #end method
        
