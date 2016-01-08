# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:45:06 2016

@author: Emanuel
"""
import os
from Bio import SeqIO
from Bio import Entrez
from Bio import SwissProt

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
    def __init__(self,typeG,strand,location,locus_tag,old_locus_tag,db_xref,product,accession,translation,seq):
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
        
    def __str__(self):
         res= " Type: "+self.type+"\n Strand: "+self.strand+"\n Location: "+ self.location+ "\n Locus_Tag: "+ self.locus_tag+"\n Old_Locus_Tag: "+ self.old_locus_tag+ "\n Db_xref: "+ self.db_xref + "\n Name: "+ self.product+ "\n Notes: "+ self.notes+ "\n Seq: " + self.seq+ "\n Accession: "+ self.accession+ "\n Translation: " + self.translation 
         return res
         
    def uniprotSearch(self):
        Entrez.email = "mafaldasap@gmail.com"     # Always tell NCBI who you are
        handle5 = Entrez.efetch(db="nucleotide", id="15638995", rettype="fasta", strand=1, seq_start=639851, 
                       seq_stop=765100)
        recordF = SeqIO.read(handle5, "fasta")
        
         
         
         
         