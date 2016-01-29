# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 17:18:22 2016
This is a small file with a class called Conn, responsible for inserting data in a MySql database.
We used pymysql to make this communication. The database and the tables where this data is inserted, was created with the help from Visual Paradigm.
Later on more tables were created using host's domain tools, such as phpMyAdmin. One good feature it has that helped us is being able to import csv files.(Excell files can be saved with this format)
@author: Emanuel
"""
import pymysql

class Conn:
   #-- The function categorization part was not successfully implemented, nonetheless the tables were created and the dictionary was there to help with the inserts
   #-- The parameters for the init are requested by pymysql to establish a connection
   def __init__(self,host="localhost",port=3306, user="root", db="lb_db"):
       self.host = host
       self.port=port
       self.user=user
       self.db = db
       self.conn = pymysql.connect(host=self.host, port=self.port, user=self.user, db=self.db)
       self.functions={"Metabolism":1, "Carriers": 2 , "Regulation":3 , "Signaling": 4, "Movement":5, "Processing":6 , "Protein_S_P":7 ,"RNA":8, "Other":9, "Unknown":10  }
       
       #-- the methods bellow are INSERT query's with the language syntax, receiving a the "missing" part of the phrase (data to be inserted).
   def insertCDS(self, qCDS, qUNI):
       cur = self.conn.cursor()
       cur.execute("INSERT INTO cds VALUES ("+qCDS+")")
       self.conn.commit()
       cur.close()
       self.conn.close()
       self.conn = pymysql.connect(host=self.host, port=self.port, user=self.user, db=self.db)
       cur = self.conn.cursor()
       cur.execute("INSERT INTO uniprot VALUES ("+qUNI+")")
       self.conn.commit()
       cur.close()
       self.conn.close()
       
       
   def insertRNA(self, qRNA):
       cur = self.conn.cursor()
#       print(qRNA)
#       print()
      
       cur.execute("INSERT INTO rna VALUES ("+qRNA+")")
       self.conn.commit()
       cur.close()
       self.conn.close()
       
   def insertRep(self, qREP):
       cur = self.conn.cursor()
       cur.execute("INSERT INTO rep_region VALUES ("+qREP+")")
       self.conn.commit()
       cur.close()
       self.conn.close()
       
   def insertBlast(self, qBlast):
       cur = self.conn.cursor()
#       print(qBlast)
#       print()
       cur.execute("INSERT INTO blastResult (type , gi , accession , score ,cover,e_value, identities , query , matchq , subj ,subjS, subjE, function , link , geneID ) VALUES ("+qBlast+")")
       self.conn.commit()
       cur.close()
       self.conn.close()
       
       
       
       
       
       