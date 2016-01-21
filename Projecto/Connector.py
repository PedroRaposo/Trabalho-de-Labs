# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 17:18:22 2016

@author: Emanuel
"""
import pymysql

class Conn:
   
   def __init__(self,host="localhost",port=3306, user="root", db="lb_db"):
       self.host = host
       self.port=port
       self.user=user
       self.db = db
       self.conn = pymysql.connect(host=self.host, port=self.port, user=self.user, db=self.db)
       
       
   def insertCDS(self, qCDS, qUNI):
       cur = self.conn.cursor()
       print(qCDS)
       print(qUNI)
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
       
       