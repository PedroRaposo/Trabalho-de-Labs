# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:45:06 2016

@author: Emanuel
"""
import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio import SwissProt
from Bio.ExPASy import ScanProsite
import Bio
import inspect
import json
import requests
import re
#import bioservices
#import uniprot
import pprint
#import urllib
import time
import sys
#import httplib
#from bioservices import UniProt



class MyRep:
    def __init__(self,typeG,strand,location,seq):
        self.type=typeG
        self.strand = strand
        self.location=location
        self.seq=seq


class MyTRNA:
     def __init__(self,typeG,strand,location,locus_tag,old_locus_tag,geneID,product,anticodon,seq):
         self.type=typeG
         self.strand = strand
         self.location=location
         self.locus_tag = locus_tag
         self.old_locus_tag=old_locus_tag
         self.geneID=geneID
         self.product = product
         self.notes=""
         self.seq=seq
         self.anticodon=anticodon
         
         
     def __str__(self):
         res= " Type: "+self.type+"\n Strand: "+self.strand+"\n Location: "+ self.location+ "\n Locus_Tag: "+ self.locus_tag+ "\n Old_Locus_Tag: "+ self.old_locus_tag+ "\n geneID: "+ self.geneID + "\n Name: "+ self.product+ "\n Notes: "+ self.notes+ "\n Seq: " + self.seq+ "\n Anticodon: "+ self.anticodon 
         return res
         
class MyCDS:
    def __init__(self,typeG,strand,location,locus_tag,old_locus_tag,db_xref,product,accession,translation,seq, ec=None, tc=None,up=1 ):
        self.type=typeG
        self.strand = strand
        self.location=location
        self.locus_tag = locus_tag
        self.old_locus_tag= old_locus_tag
        self.db_xref=db_xref
        self.product = product
        self.accession=accession
        self.translation=translation
        self.notes=""
        self.seq=seq
        self.ec = ec
        self.tc = tc
        tag= re.findall("[0-9]+",db_xref)
        if(len(tag)>1):        
            self.gi = tag[0]
            self.geneID = tag[1]
        else:
            self.geneID = tag[0]
        tag= re.search("[^\['].+[^\]']",self.old_locus_tag)
        self.old_locus=tag.group()
        self.similarity=""
        self.catalytic_Activity=""
        self.function=""
        self.subunit=""
        self.cofactor=""
        self.seqCaution=""
        self.subcelularLoc=""
        self.review=""
        self.goList=[]
        self.uniprotSearch(up)
        
    def __str__(self):
         uniprotInfo=""
         if(self.catalytic_Activity !=""):
             uniprotInfo+="\n Catalytic Activity: " + self.catalytic_Activity
         if(self.cofactor !=""):
             uniprotInfo+="\n Cofactor: " + self.cofactor
         if(self.function !=""):
             uniprotInfo+="\n Function: " + self.function
         for go in self.goList:
             uniprotInfo+="\n " + str(go)
         if(self.molecular_weight !=""):
             uniprotInfo+="\n Molecular weight: " + str(self.molecular_weight)
         for pA in self.prot_accessions :
             uniprotInfo+="\n Uniprot Accession: " + pA
         if(self.review !=""):
             uniprotInfo+="\n Grade of Revision: " + self.review
         if(self.similarity !=""):
             uniprotInfo+="\n Similarities: " + self.similarity
         if(self.subcelularLoc !=""):
             uniprotInfo+="\n  " + self.subcelularLoc
         if(self.seqCaution !=""):
             uniprotInfo+="\n Area of caution: " + self.seqCaution
         if(self.subunit !=""):
             uniprotInfo+="\n " + self.subunit
        
             
         res= " Type: "+self.type+"\n Strand: "+self.strand+"\n Location: "+ self.location+ "\n Locus_Tag: "+ self.locus_tag+"\n Old_Locus_Tag: "+ self.old_locus_tag+ "\n Db_xref: "+ self.db_xref + "\n Name: "+ self.product+ "\n Notes: "+ self.notes+ "\n Seq: " + self.seq+ "\n Accession: "+ self.accession+ "\n Translation: " + self.translation + "\n EC_number: " + self.ec+ "\n TC_number: "+ self.tc+uniprotInfo
         return res
 
   
    #Method uniprotSearch, add relevant information to the MyCDS object from
    #the uniprot website. The method functions by downloading a txt obtained
    #from a query to uniprot's REST service, afterwards the file is stored and
    #parsed. A flag update is used to control the download file by need.     
    def uniprotSearch(self,update=0):
        tag= re.search("[^\['].+[^\]']",self.old_locus_tag)
          
        #If it´s neccessary to download the txt file from uniprot
        if(update==1):
            base_url= "http://www.uniprot.org/uniprot"
            payload={'query':'gene:'+tag.group(), 'format':'txt'}
    
            result= requests.get(base_url,params=payload)
            if result.ok:
                cont = str(result.content)
                advance=0
                flag=1
                linha=""
                #writes the file with the format needed for SwissProt.parse()
                with open(tag.group()+".dat", "w") as text_file:
                    for word in cont:
                        if(advance<2):
                            advance+=1
                        else:
                            if(word=='n' and flag==2):
                                print(linha, file=text_file)
                                flag=1
                                linha=""
                            else:
                                if(flag==2 and word!='n'):
                                    flag=1
                                    linha+="\\"+word
                        
                                if(word=='\\' and flag==1):
                                    flag+=1
                        
                                if(word!='\\' and flag==1):
                                    linha+= word
       
            else:
                #error message
                print("Communication to unipot failed on ",self.geneID)
            
        for record in SwissProt.parse(open(tag.group()+".dat")):
                #store List of the accession numbers, e.g. [‘P00321’]                 
                self.prot_accessions = record.accessions
                #store versified information
                for comment in record.comments:
                    if(comment.find("FUNCTION:")>=0):
                        self.function+=comment
                    if(comment.find("SUBCELLULAR LOCATION:")>=0):
                        self.subcelularLoc = comment
                    if(comment.find("SIMILARITY:")>=0):
                        self.similarity += comment
                    if(comment.find("CATALYTIC ACTIVITY:")>=0):
                        self.catalytic_Activity+=comment
                    if(comment.find("SUBUNIT:")>=0):
                        self.subunit+=comment
                    if(comment.find("SEQUENCE CAUTION")>=0):
                        self.seqCaution= comment
                    if(comment.find("COFACTOR:")>=0):
                        self.cofactor = comment
                #print(record.comments)
                #print(record.molecule_type) // maybe not relevant
                #print(record.organelle) // maybe not relevant
                #print(record.features) //USEFULNESS NOT ANALYSED YET
                #print(record.description)// maybe not relevant
                self.review = record.data_class 
                self.molecular_weight = record.seqinfo[1]
                #fetch GO 
                for reference in record.cross_references:
                    if(reference[0]=="GO"):
                        self.goList.append(reference)
                
                
                
        
    def blast(self,update=0,database="nt"):
        header = ">"+self.old_locus_tag
        name= self.old_locus_tag+".fasta"
        tag = re.search("[^\['].+[^\]']",self.db_xref)
        prot= tag.group().replace(":","-")
        nameProt= prot+"_prot.fasta"
        headerProt= ">"+prot
        #save seqs as fasta--------------------------------
        #seq nucl
        if(update==1):
            with open(name, "w") as text_file:
                print(header, file=text_file)
                print(self.seq,file=text_file)
        #seq prot
            if(self.translation!=None):
                translation= re.search("[^\['].+[^\]']" , self.translation).group()
                with open(nameProt, "w") as text_file:
                    print(headerProt, file=text_file)
                    print(translation,file=text_file)
        
        #blast-----------------------------------------
        
        #record=SeqIO.read(open(name), format="fasta")

        #result_handle=NCBIWWW.qblast("blastn", database, record.format("fasta"))
        
        #save_file = open(self.old_locus_tag+"_"+database+".xml", "w")
        #save_file.write(result_handle.read())
        #save_file.close()
        #result_handle.close()
        #open XML-------------------------------------------------
        
        result_handle = open(self.old_locus_tag+"_"+database+".xml")
        blast_record = NCBIXML.read(result_handle)
        print("Global results: ")
        print("Database:")
        print(blast_record.database)
        print("Substitution matrix: ")
        print(blast_record.matrix)
        print("Gap penalties: ")
        print(blast_record.gap_penalties)
        print("----------------------------------------------- ")
        print(" ")
        #print("------------------------------------------------- ")
        for i in range(len(blast_record.alignments)) :
            alignment = blast_record.alignments[i]
            print("Alignment number :", i+1 )
            print("ascenssion: ")
            print(alignment.accession)
            print("")
            for i2 in range(len(alignment.hsps)):
                alignment_hsp = alignment.hsps[i2]
                print("Hsp ", i2+1, ":")
                print("E-value: ")
                print(alignment_hsp.expect)
                print("Size: ")
                print(alignment_hsp.align_length)
                print("")
            print("-------------next-----------------------------")
#        ac_numb = []
#        organism = []
#        evalues = []
#        for blast_record in blast_records:
#        	for alignment in blast_record.alignments:
#        			organism.append(alignment.hit_id)
#        			ac_numb.append(alignment.accession)
#        			for hsp in alignment.hsps:
#        				evalues.append(hsp.expect)
#            
#        print ("Acession numbers:\n", str(ac_numb), "\n\n", "Organismos:\n", str(organism), "\n\n", "E-values:\n" , str(evalues))
#        
#        result2 = open("my_blast_sc.xml")
#        blast_sc = NCBIXML.read(result2)
#        blast1=blast_sc.alignments[0]  
#        hsp=blast1.hsps[0]
#        
        
        
        #==============================================================================
#        print ("Alinhamento com menor E-value:",hsp.expect, "\n")
#        print ("Accession:", blast1.accession, "\n")
#        print ("Hit ID:", blast1.hit_id, "\n")
#        print ("Definição:", blast1.hit_def, "\n")
#        print ("Comprimento do alinhamento:", blast1.length, "\n")
#        print ("Numero de HSPs:", len(blast1.hsps), "\n")
        #==============================================================================
        
                
                
                
                
                 
                 
                 
         