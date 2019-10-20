#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:07:11 2019

@author: jean
"""
import sys

###########################
#
#      fonctions
#
###########################

def PDBabundanceDispo(infile_genes, infile_abundance, flag_construct_correlation):
    
    pdb_genes = {}
    pdb_outlist = []
    
    #ouverture du fichier
    f = open(infile_genes, "r")
    lines = f.readlines()
    f.close()
    
    
    for line in lines:
        if not line[0] == "#":
            line = line.strip()
    
        tmp_gene_name = line.split("\t")[0]
        tmp_gene_PDB = line.split("\t")[1]
    
        pdb_genes[tmp_gene_name] = tmp_gene_PDB
    
    f = open(infile_abundance, "r")
    lines = f.readlines()
    f.close()
    
    if(flag_construct_correlation == 1):
        correlation_file = open("correlation_file.txt", "w")
        correlation_file.write('#gene_OLN\tPDB_ID\tabundance\n')

    for line in lines:
        line = line.strip()
        
        if not (line[0] == "#"):
            tmp_gene = line.split("\t")[1].split(".")[1]
            
            for tmp_gene_name in pdb_genes:
                
                if (tmp_gene == tmp_gene_name) :
                    
                    pdb_outlist.append(pdb_genes[tmp_gene_name])
                    
                    if(flag_construct_correlation == 1):
                        tmp_abundance = line.split("\t")[2]
                        correlation_file.write(tmp_gene_name + '\t' + pdb_genes[tmp_gene_name] + '\t' + tmp_abundance + '\n')
    
    return pdb_outlist


###########################
#
#      MAIN
#
###########################

if ((len(sys.argv) != 4) or ("-h" in sys.argv)) :
    print("\nUsage :", sys.argv[0], "\n\t<infile with gene/pdb correlation>\n\t<infile with gene/abundance correlation>\n\t<Y/N> (wether to construct or not a correlation file with gene / PDB_ID / Abundance)\n")
else :
    
    infile_genes = sys.argv[1]
    try:
        f = open(infile_genes, "r")
        f.close()
    except:
        print('Unable to read file :', infile_genes)
        sys.exit()        
    
    infile_abundance = sys.argv[2]
    try:
        f = open(infile_abundance, "r")
        f.close()
    except:
        print('Unable to read file :', infile_abundance)
        sys.exit()        
    
    
    if(sys.argv[3] in ["YES", "yes", 'y', 'Y']) :
        flag_construct_correlation = 1
            
    elif(flag_construct_correlation in ["NO", "no", "N", "n"]) :
        flag_construct_correlation = 0
    
    else :
        print("invalid argument syntax: flag_construct_correlation must be 'Y' or 'N'") 
    
    print('#PDB_ID')
    for tmp_pdb in PDBabundanceDispo(infile_genes, infile_abundance, flag_construct_correlation):
        print(tmp_pdb)