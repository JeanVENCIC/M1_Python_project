#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:37:32 2019

@author: jean
"""

import sys, os, urllib.request, time

###########################
#
#      fonctions
#
###########################

def Download_PDBlist(pdb_list, directory_path):
    
    files_already_downloaded = os.listdir(directory_path)
    nb_treatedPDB = 0
    nb_pdbToDownload = len(pdb_list)
    
    for tmp_pdb in pdb_list:
        nb_treatedPDB += 1
        print('Treating ', tmp_pdb, '.pdb ...\t[', nb_treatedPDB, '/', nb_pdbToDownload, ']', sep ='')
        
        if not (tmp_pdb+'.pdb' in files_already_downloaded):
            try :
                urllib.request.urlretrieve(('https://files.rcsb.org/download/' + tmp_pdb + '.pdb'), (directory_path + '/' + tmp_pdb + '.pdb'))
                time.sleep(1)
            except:
                print('Unable to retrieve ' + tmp_pdb + '.pdb file from https://www.rcsb.org/structure/' + tmp_pdb + '#')
                continue
            
        else:
            print(tmp_pdb + '.pdb already in directory')

###########################
#
#      MAIN
#
###########################

if ((len(sys.argv) != 3) or ("-h" in sys.argv)) :
    print("Usage :", sys.argv[0], "\n\t<infile with pdb files to download>\n\t<directory path where pdb files will be downloaded>\n")
else :

    infile_pdb_list = sys.argv[1]
    try:
        f = open(infile_pdb_list, "r") 
        f.close()
    except:
        print('Unable to read file :', infile_pdb_list)
        sys.exit()


    directory_path = sys.argv[2]
    try:
        os.path.isdir(directory_path)
    except:
        print('Unable to find given directory :', directory_path)
        sys.exit()
        

    lines = f.readlines()
    f.close()
    
    pdb_list = []
        
    for line in lines :
        if not (line[0] == "#"):
            line = line.strip()
            
            pdb_list.append(line)
    
    Download_PDBlist(pdb_list, directory_path)
        