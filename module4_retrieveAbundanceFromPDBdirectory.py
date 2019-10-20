#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:25:38 2019

@author: jean
"""
import sys, os, warnings

###########################
#
#      fonctions
#
###########################

def AbundanceFromPDBdirectory(pdb_files, correlation_file):
    
    f = open(correlation_file, "r")
    lines = f.readlines()
    f.close()
    
    print('#gene_OLN\tPDB_ID\tabundance')
    
    for line in lines:
        line = line.strip()
        
        if(line[0] == '#'):
           if(line != '#gene_OLN\tPDB_ID\tabundance'):
              warnings.warn("\ncorrelation_file given in module4_retrieveAbundanceFromPDBdirectory.py has not usual header :\n\t\t#gene_OLN\tPDB_ID\tabundance\n", Warning)
    
        else:
            if(line.split('\t')[1] in pdb_files):
                print(line)
###########################
#
#      MAIN
#
###########################

if ((len(sys.argv) != 3) or ("-h" in sys.argv)) :
    print("\nUsage : ", sys.argv[0], "\n\t<directory path where pdb files are stored>\n\t<correlation file : gene / PDB / abundance>\n")
else :
    directory_path = sys.argv[1]
    try:
        os.path.isdir(directory_path)
    except:
        print('Unable to find given directory :', directory_path)
        sys.exit()        
    
    correlation_file = sys.argv[2]
    try:
        f = open(correlation_file, "r")
        f.close()
    except:
        print('Unable to read file :', correlation_file)
        sys.exit()
    
    pdb_files = os.listdir(directory_path)
    pdb_list = []
    
    for tmp_pdb in pdb_files:
        pdb_list.append(tmp_pdb[:tmp_pdb.find(".")])

    AbundanceFromPDBdirectory(pdb_list, correlation_file)