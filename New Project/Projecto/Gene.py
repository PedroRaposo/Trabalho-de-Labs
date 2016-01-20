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
import Connector
#import bioservices
#import uniprot
import pprint
#import urllib
import time
import sys
#import httplib
#from bioservices import UniProt

class MyHsp:
    def __init__(self,gi,accession,bitScore,cover,e_value,ident,isPdb=False):
        self.gi=gi
        self.accession=accession
        self.bitScore=bitScore
        self.cover = cover        
        self.e_value= e_value
        self.ident=ident
        self.pdbLink=""
        self.function=""
        self.isPdb=isPdb
        if(isPdb==True):
            self.ncbiLink = "http://www.ncbi.nlm.nih.gov/protein/"+gi
            Entrez.email = "emanuel_queiroga1@hotmail.com"    # Its necessary to tell NCBI who we are
            #Using the Entrez import we fetch the genbank file
            if(os.path.isfile("homologous_"+gi+".gb")!=True):
                handleGb = Entrez.efetch(db="protein", id=gi, rettype="gb")
                save_fileGb = open("homologous_"+gi+".gb", "w")
                save_fileGb.write(handleGb.read())
                save_fileGb.close()
                handleGb.close()
            recordGb = SeqIO.read("homologous_"+gi+".gb", "genbank")
#            self.function=recordGb.annotations            
            procura = re.search("class:",recordGb.annotations["db_source"])
            num= procura.span()
            res = ""
            sent = True
            seq2 = recordGb.annotations["db_source"][int(num[1]+1):len(recordGb.annotations["db_source"])]
            for i in seq2:
                if i == ";":
                    sent = False
                if sent == True:
                    res += i
            self.function=res
    
    def __str__(self):
        res= " Gi: " + self.gi + " \n Accession: "+ self.accession + " \n Score: "+ str(self.bitScore) + "\n E-value: " + str(self.e_value) + "\n Identities: " + str(self.ident)
        if(self.isPdb == True):
            res+= "\n Function: "+ self.function + "\n Link to PDB: "+ self.ncbiLink
        return res
    
class MyRep:
    def __init__(self,typeG,strand,location,seq):
        self.type=typeG
        self.strand = strand
        self.location=location
        self.seq=seq
#        c = Connector.Conn()
#        s="'"+self.type+"','"+self.strand+"','"+self.location+"'"
#        c.insertRep(s)
#        
        
    def __str__(self):
        res= " Type: "+self.type+"\n Strand: "+self.strand+"\n Location: "+ self.location + "\n Seq: " + self.seq
        return res        

class MyTRNA:
     def __init__(self,typeG,strand,location,locus_tag,old_locus_tag,geneID,product,anticodon,seq):
         self.type=typeG
         self.strand = strand
         self.location=location
         self.locus_tag = locus_tag
         self.old_locus_tag=old_locus_tag
         tag= re.findall("[0-9]+",geneID)
         if(len(tag)>1):        
            self.gi = tag[0]
            self.geneID = tag[1]
         else:
            self.geneID = tag[0]
            self.gi="NULL"
         #self.geneID=geneID
         self.product = product
         self.notes=""
         self.seq=seq
         self.anticodon=anticodon
         s="'"+self.geneID+"', '"+self.type +"', '"+self.strand+"' ,'"+self.location+"' , '"+self.locus_tag.replace("'", " ")+"' , '" + self.old_locus_tag.replace("'", " ")+"' , '"+self.product+"' , NULL , '"+self.anticodon+"'"
         c = Connector.Conn()
         c.insertRNA(s)
         
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
        if(self.ec=="None"):
            self.ec="NULL"
        if(self.tc=="None"):
            self.tc="NULL"
        
        tag= re.findall("[0-9]+",db_xref)
        if(len(tag)>1):        
            self.gi = tag[0]
            self.geneID = tag[1]
        else:
            self.geneID = tag[0]
            self.gi="NULL"
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
        self.prot_accessions=[]
        self.uniprotSearch(up)
        self.goString=""
        self.uniA=""
        if(len(self.goList)==0):
            self.goString="NULL"
        for go in self.goList:
            self.goString+= str(go).replace("'", " ") +";"
        for a in self.prot_accessions:
            self.uniA+= a+";"
        self.blastPInfo=None
        self.hits=[]
        #c = Connector.Conn()
        #c.insertCDS(self.sqlStringCDS(),self.sqlStringUni())
        #s=""
        
        #self.blast()
        #self.parseALLBlast()
#        self.hitsToFile()
        
        
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
         if(self.similarity !="NULL"):
             uniprotInfo+="\n Similarities: " + self.similarity
         if(self.subcelularLoc !=""):
             uniprotInfo+="\n  " + self.subcelularLoc
         if(self.seqCaution !=""):
             uniprotInfo+="\n Area of caution: " + self.seqCaution
         if(self.subunit !=""):
             uniprotInfo+="\n " + self.subunit
        
             
         res= " Type: "+self.type+"\n Strand: "+self.strand+"\n Location: "+ self.location+ "\n Locus_Tag: "+ self.locus_tag+"\n Old_Locus_Tag: "+ self.old_locus_tag+ "\n Db_xref: "+ self.db_xref + "\n Name: "+ self.product+ "\n Notes: "+ self.notes+ "\n Seq: " + self.seq+ "\n Accession: "+ self.accession+ "\n Translation: " + self.translation + "\n EC_number: " + self.ec+ "\n TC_number: "+ self.tc+uniprotInfo
         return res
         
    def sqlStringCDS(self):
        qcds="'"+self.geneID+"' , "
        qcds+="'"+self.type.replace("'", " ")+"' ,'"+self.strand+"' , '"+ self.location.replace("'", " ")+ "' , '"+ self.locus_tag+"' , '"+ self.old_locus_tag[2:-2]+ "' , "
        if(self.gi!="NULL"):
            gi="'"+self.gi.replace("'", " ")+"'" 
        else:
            gi = self.gi
            
        if(self.tc!="NULL"):
            tc="'"+self.tc.replace("'", " ")+"'" 
        else:
            tc = self.tc
        if(self.ec!="NULL"):
            ec="'"+self.ec.replace("'", " ")+"'" 
        else:
            ec = self.ec
            
        qcds+=gi+" ,  '"+ self.product.replace("'", " ")+ "' ,  '"+ self.accession.replace("'", " ")+ "' ,"
            
        qcds+=ec+" , "+tc+ ", NULL, NULL, NULL" 
        return qcds
        
    def sqlStringUni(self):    
        uniprotInfo="'"+self.geneID+"' ,"
        
        if(self.catalytic_Activity !=""):
             uniprotInfo+="'" + self.catalytic_Activity.replace("'", " ")+ "' , "
        else:
             uniprotInfo+=" NULL ,"
             
        if(self.cofactor !=""):
             uniprotInfo+="'" + self.cofactor.replace("'", " ")+ "' , "
        else:
             uniprotInfo+=" NULL ,"
             
        if(self.function !=""):
             uniprotInfo+="'" + self.function.replace("'", " ")+ "' , "
        else:
             uniprotInfo+=" NULL ,"
             
        if(self.goString !="NULL"):
             uniprotInfo+="'" + self.goString+ "' , "
        else:
            uniprotInfo+=" NULL ,"

        uniprotInfo+="'" + str(self.molecular_weight)+ "' , "       
        uniprotInfo+="'" + self.uniA.replace("'", " ")+ "' , "
        uniprotInfo+="'" + self.review+ "' , "
        if(self.similarity !=""):
             uniprotInfo+="'" + self.similarity.replace("'", " ")+ "' , "
        else:
             uniprotInfo+=" NULL ,"
        if(self.subcelularLoc !=""):
             uniprotInfo+="'" + self.subcelularLoc.replace("'", " ")+ "' , "
        else:
             uniprotInfo+=" NULL ,"
             
        if(self.seqCaution !=""):
             uniprotInfo+="'" + self.seqCaution.replace("'", " ")+ "' , "
        else:
             uniprotInfo+=" NULL ,"
             
        if(self.subunit !=""):
             uniprotInfo+="'" + self.subunit.replace("'", " ")+ "'"
        else:
             uniprotInfo+=" NULL"
        
        return uniprotInfo
        
        
        
        
 
   
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
                
                
                
        
    def blast(self,update=0,databaseP="pdb",database="nt"):
        header = ">"+self.old_locus_tag
        name= self.old_locus_tag+".fasta"
        prot= self.geneID
        nameProt= prot+"_prot.fasta"
        headerProt= ">"+prot
        #save seqs as fasta--------------------------------
        #seq nucl

#        with open(name, "w") as text_file:
#            print(header, file=text_file)
#            print(self.seq,file=text_file)
        #seq prot
#        if(self.translation!="None"):
#            #translation= re.search("[^\['].+[^\]']" , self.translation).group()
#            with open(nameProt, "w") as text_file:
#                print(headerProt, file=text_file)
#                print(self.translation,file=text_file)
        
        #blast-----------------------------------------
        if(os.path.isfile(prot+"_"+databaseP+".xml")!=True):
            record=SeqIO.read(open(name), format="fasta")
    
            result_handle=NCBIWWW.qblast("blastn", database, record.format("fasta"))
            
            save_file = open(self.old_locus_tag+"_"+database+".xml", "w")
            save_file.write(result_handle.read())
            save_file.close()
            result_handle.close()
        
#        #blast protif(os.path.isfile(prot+"_"+databaseP+".xml")!=True):
#        if(os.path.isfile(prot+"_"+databaseP+".xml")!=True):   
#            if(self.translation!="None"):
#                record=SeqIO.read(open(nameProt), format="fasta")
#        
#                result_handle=NCBIWWW.qblast("blastp", databaseP, record.format("fasta"))
#                time.sleep(3)
#                if(result_handle.readable()):
#                    save_file = open(prot+"_"+databaseP+".xml", "w")
#                    save_file.write(result_handle.read())
#                    save_file.close()
#                    result_handle.close()
#                else:
#                    print("Error in ", nameProt)
        
        #open XML-------------------------------------------------
    def parseBlast(self,database="pdb", molecule="prot", e_param=0.05, i_param=30):
        if(molecule=="prot"):
            name = self.geneID+"_"+database+".xml"
            tam = len(self.translation)
        if(molecule=="nucl"):
            name= self.old_locus_tag+"_"+database+".xml"
            tam = len(self.seq)
        if(os.path.isfile(name)==True):    
            result_handle = open(name)
            blast_record = NCBIXML.read(result_handle)
    #        print("Global results: ")
    #        print("Database:")
    #        print(blast_record.database)
    #        print("Substitution matrix: ")
    #        print(blast_record.matrix)
    #        print("Gap penalties: ")
    #        print(blast_record.gap_penalties)
    #        print("----------------------------------------------- ")
    #        print(" ")
            #print("------------------------------------------------- ")
            for i in range(len(blast_record.alignments)) :
                alignment = blast_record.alignments[i]
                hitAccession = alignment.accession;
                tag = re.search("[0-9]+",alignment.hit_id)
                hitGi = tag.group() 
                self.hits.append([])
                for i2 in range(len(alignment.hsps)):               
                    alignment_hsp = alignment.hsps[i2]
                    rang = (alignment_hsp.align_length/tam)*100
                    if(molecule=="prot"):
                       if(alignment_hsp.expect<=e_param and alignment_hsp.identities>=i_param):
                           hsp= MyHsp(hitGi,hitAccession,alignment_hsp.bits,rang,alignment_hsp.expect,alignment_hsp.identities,True)
                           self.hits[i].append(hsp)
                    else:
                        if(alignment_hsp.expect<=e_param and alignment_hsp.identities>=i_param):
                            hsp= MyHsp(hitGi,hitAccession,alignment_hsp.bits,rang,alignment_hsp.expect,alignment_hsp.identities)
                            self.hits[i].append(hsp)
            
            self.hitsToFileNT()
            

            
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
        
    def hitsToFile(self):
        hitIndice = 1
        hspIndice = 1
        if(len(self.hits)>0):
            if(len(self.hits[0])==0 ):
                with open("failedBlastPwith_0.05_30.txt","a") as text_file:
                    print(self.geneID,file=text_file )
                text_file.close()
            else:            
                name = "Blast_Results_"+ self.geneID+ ".txt"
                with open(name, "w") as text_file:
                    for hit in self.hits:
                        hspIndice = 1
                        if(len(hit)>0):
                            ni= "Hit num: "+str(hitIndice)
                            print(ni, file=text_file)
                            print("",file=text_file)
                            for hsp in hit:
                                n="Hsp num :"+str(hspIndice)
                                print(n, file=text_file)
                                print(hsp, file=text_file)
                                print("------------------------------------------------------------------",file=text_file)
                                print("",file=text_file)
                                hspIndice+= 1    
                        hitIndice+=1
                        print("",file=text_file)
                text_file.close()
        else:
            with open("failedBlastP2.txt","a") as text_file:
                print(self.geneID,file=text_file )
            text_file.close()
            
    def hitsToFileNT(self):
        hitIndice = 1
        hspIndice = 1
        if(len(self.hits)>0):
            if(len(self.hits[0])==0 ):
                with open("failedBlastNTwith_0.05_30.txt","a") as text_file:
                    print(self.geneID,file=text_file )
                text_file.close()
            else:            
                name = "BlastNT_Results_"+ self.geneID+ ".txt"
                with open(name, "w") as text_file:
                    for hit in self.hits:
                        hspIndice = 1
                        if(len(hit)>0):
                            ni= "Hit num: "+str(hitIndice)
                            print(ni, file=text_file)
                            print("",file=text_file)
                            for hsp in hit:
                                n="Hsp num :"+str(hspIndice)
                                print(n, file=text_file)
                                print(hsp, file=text_file)
                                print("------------------------------------------------------------------",file=text_file)
                                print("",file=text_file)
                                hspIndice+= 1    
                        hitIndice+=1
                        print("",file=text_file)
                text_file.close()
        else:
            with open("failedBlastNT2.txt","a") as text_file:
                print(self.geneID,file=text_file )
            text_file.close()
            
    def parseALLBlast(self, e_param=0.05, i_param=30):
        databaseP="pdb"
        moleculeP="prot" 
        databaseN="nt"
        moleculeN="nucl"
        name1 = self.geneID+"_"+databaseP+".xml"
        tam1 = len(self.translation)
        
        name2= self.old_locus_tag+"_"+databaseN+".xml"
        tam2 = len(self.seq)
        
        if(os.path.isfile(name1)==True):
            name=name1
            tam=tam1
            molecule=moleculeP
            
        
        if(os.path.isfile(name2)==True):
            name=name2
            tam=tam2
            molecule=moleculeN
        
        
        if(os.path.isfile(name)==True):    
            result_handle = open(name)
            blast_record = NCBIXML.read(result_handle)
    #        print("Global results: ")
    #        print("Database:")
    #        print(blast_record.database)
    #        print("Substitution matrix: ")
    #        print(blast_record.matrix)
    #        print("Gap penalties: ")
    #        print(blast_record.gap_penalties)
    #        print("----------------------------------------------- ")
    #        print(" ")
            #print("------------------------------------------------- ")
            for i in range(len(blast_record.alignments)) :
                alignment = blast_record.alignments[i]
                hitAccession = alignment.accession;
                tag = re.search("[0-9]+",alignment.hit_id)
                hitGi = tag.group() 
                self.hits.append([])
                for i2 in range(len(alignment.hsps)):               
                    alignment_hsp = alignment.hsps[i2]
                    rang = (alignment_hsp.align_length/tam)*100
                    if(molecule=="prot"):
                       if(alignment_hsp.expect<=e_param and alignment_hsp.identities>=i_param):
                           hsp= MyHsp(hitGi,hitAccession,alignment_hsp.bits,rang,alignment_hsp.expect,alignment_hsp.identities,True)
                           self.hits[i].append(hsp)
                    else:
                        if(alignment_hsp.expect<=e_param and alignment_hsp.identities>=i_param):
                            hsp= MyHsp(hitGi,hitAccession,alignment_hsp.bits,rang,alignment_hsp.expect,alignment_hsp.identities)
                            self.hits[i].append(hsp)
            if(molecule=="nucl"):
                self.hitsToFileNT()
            else:
                self.hitsToFile()
        
                
                
                 
                 
                 
         