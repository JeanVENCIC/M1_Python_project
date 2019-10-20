# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:52:41 2019

@author: Pierre
"""

"""
fonction pour gérer les fichiers données sur dokeos
"""


import os

slash = '/'
if(os.name=='nt'):
    slash='\\'

#si utilisation des fichiers de la profs:correction
def correction(file):
    
    f=open(file,"r")
    out = open("correction.txt", "w")
    
    lines = f.readlines()
    
    for line in lines:
        if(line[0]!='#'):
           line = line.strip("\n")
           tmp=line.split("\t")
           tmp2 = tmp[0].split("_")
           out.write(tmp2[0] + "\t" + tmp[1] + "\t" + tmp[2] + "\t" +tmp[3] +"\n")
        else:
            out.write(line)
            
    f.close()
    out.close()
    
    
def rename(directoryIn):
    
    print("lancement du rename\n")
    
    pathInitial = os.getcwd()
    listTmp= []
    listFile=os.listdir(directoryIn)
    os.chdir(directoryIn)
    
    for file in listFile:
        try:
            rename = file.split("_")[0] + ".pdb"
            
            os.rename(file,rename)
        except:
            os.remove(file)
    os.chdir(pathInitial)    
        
    print("\nfin de la rename\n")