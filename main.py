# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 16:08:28 2016

@author: Emanuel
The main class that controls a big chunk of python scripts. (Files Data , Gene and Connector, an image labelled "Architecture.png" 
has a simple graphic description of this classes and their relation with other parts of this project )
Comments of procedure and information in this file will be of this type #--
Comments with only # are methods that can be commented/uncommented according to necessity
There are some files that we only need to print once, for debug or make intermediary analyses or for small scripts not controlled from this main file.
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
    print("files retrieved successfully!")
    #--end--
    #--Next step retrieve the information necessary to process from those files
    print("retrieving information from Treponema files")
	#--Data is where objects from classes MyCDS , MyTRNA and MyRep are created and stored in lists.  
    data.getData(0)
    print("information retrieved successfully!")
    #--if needed for debug print the information in data comment line if it´s not needed
    #--data.printData()
    #--print only rna reladed genes
    #data.printRNA()
    #--print only CDS
    #data.printCDS()
    #--print only repeat regions
    #data.printRepRegion()
    #--print CDS tofile
    #data.CDStoFile()
	#
    #--make blast, blast can be made from multiple methods present in class gene, to suit needs.
	#--Blast can take a while to finish and it isn't good to overload the server in which the blast is requested, so we decided to make a
    #-- a semi-automatic process where we would change which class calls the blast and what method to call, the method below is an example we used.
    data.makeBlastAll()
    
    

   
    
    
    