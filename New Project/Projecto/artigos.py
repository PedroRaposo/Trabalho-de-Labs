# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 19:36:47 2016

@author: Emanuel
"""

from Bio import Entrez
from Bio import SeqIO


Entrez.email = "mafaldasap@gmail.com"     # Always tell NCBI who you are
handle = Entrez.egquery(term="Treponema pallidum subsp. pallidum")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="pubmed":
        print(row["Count"]) # conta o número de artigos para Treponema pallidum subsp. pallidum
        

handle2 = Entrez.esearch(db="pubmed", term="Treponema pallidum subsp. pallidum", retmax=159)
record2 = Entrez.read(handle2)
IdList = record2["IdList"] # cria lista com IDs dos artigos
print(IdList)
IdList20 = record2["IdList"][:20]
print(IdList20) # cria lista com IDs dos 20 primeiros artigos



from Bio import Medline
handle3 = Entrez.efetch(db="pubmed", id=IdList20, rettype="medline",
                           retmode="text")
records = Medline.parse(handle3)
records = list(records)



#for record in records:
#    print("ID:", record.get("PMID", "?"))
#    print("Title:", record.get("TI", "?"))
#    print("Authors:", record.get("AU", "?"))
#    print("Journal Title:", record.get("JT", "?"))
#    print("Date of Publication:", record.get("DP", "?"))
#    print("Issue:", record.get("IP", "?"))
#    print("Volume:", record.get("VI", "?"))
#    print("Pagination:", record.get("PG", "?"))
#    print("Abstract:", record.get("AB", "?"))
#    print("")
i=1
with open("pubFiles.txt", "w") as text_file:
    for record in records:
        idP= record.get("PMID", "?")
        title =  record.get("TI", "?")
        authors = record.get("AU", "?")
        jornalT = record.get("JT", "?")
        date =  record.get("DP", "?")
        issue =  record.get("IP", "?")
        volume =  record.get("VI", "?")
        pagination =  record.get("PG", "?")
        abstract =  record.get("AB", "?")
        header = '<p class="text-justify">['+str(i)+']   PDMI: '+idP+' ,  <strong>'+str(authors).replace("[","").replace("]","").replace("'","")+'</strong> , <strong class="text-primary">'+title+'</strong>'
        print(header, file=text_file)
        sub = ' '+jornalT + ', issue: '+issue+ ' volume - '+ volume + ' pp. '+pagination+ ' ('+date+')</p>' 
        print(sub, file=text_file)
        ab= '<p class="text-justify"> <strong class="text-primary">Abstract: </strong>'+abstract+'</p>'        
        print(ab, file=text_file)        
        print("",file=text_file)
        i+=1
    text_file.close()
    
    
filename = "pubmed_lista.txt"
save_file = open(filename, "w")
for line in records:
    save_file.write(str(line))
save_file.close()