# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 16:08:28 2016

@author: Emanuel
Coments in this file will be of this type #--
Coments with only # are methods that can be commented/uncommented according to necessity
Later we may group methods in tasks like: retrieve all files, update all etc...
"""
import Data
import Gene


if __name__ == '__main__':
    ##Initial information
    filename="g6_Treponema"    
    db="nucleotide"
    gi="15638995"
    strand=1 
    email="emanuel_queiroga1@hotmail.com"
    seq_start=639851 
    seq_stop=765100
    
    ##--end 
    print("Work status:")
    #--First step of the work get the files from the Treponema	 pallidum	
    #-- subsp.	 pallidum	 str.	 Nichols  
    #--If the u already got the files and dont need to update them comment the following lines
    data= Data.Data(filename)
    print("getting the files of Treponema pallidum's genoma")
    data.writeFiles(db,gi,strand,email,seq_start,seq_stop)
    print("files retrived successfully!")
    #--end--
    #--Next step retrieve the information necessary to process from those files
    print("retrieving information from Treponema files")
    data.getData(0)
    print("information retrieved successfully!")
    #--if neededfor debug print the information in data comment line if it´s not needed
    #data.printData()
    #print only rna retaded genes
    #data.printRNA()
    #print only CDS
    #data.printCDS()
    #print only repeat resgions
    #data.printRepRegion()
    #print CDS tofile
    #data.CDStoFile()
    #make blast
    #data.makeBlastAll()
    
    

   
    
    
    
