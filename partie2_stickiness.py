#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 15:26:39 2019

@author: jean vencic
"""

import sys, PDBTools, numpy, stickiness, CV

###########################
#
#      fonctions
#
###########################

def printUsage():
    print("\nUsage :\n", sys.argv[0], "<input PDB file> <output PDB file> <stickiness / AA correlation file> <rc threshold used to compute CV value> <CV threshold to determine surface residues>\n")

###########################
#
#      MAIN
#
###########################
if("-h" in sys.argv):
    printUsage()
    
#Call management
if(len(sys.argv) == 6):
    stickiness_file = sys.argv[3]
    rc_threshold = sys.argv[4]
    CV_threshold = sys.argv[5]


else:
    printUsage()
    sys.exit()
    

#Argument Management
print('Gestion des arguments ...\n')
#Test du fichier d'input
infile = sys.argv[1]
try:
    f = open(infile, "r")
    f.close()
except:
    print('Unable to read file :', infile)
    sys.exit()

#Test du fichier d'output
outfile = sys.argv[2]
try:
    f = open(outfile, "w")
    f.close()
except:
    print('Unable to write in file :', outfile)
    sys.exit()

#Test du fichier stickiness.txt
try:
    dicoSticky= stickiness.setStickinessAA(stickiness_file)
except:
    print('Unable to read file :', stickiness_file)
    sys.exit()

#Test de la variable rc_threshold
try:
    rc_threshold = float(rc_threshold)
except:
    print("rc_treshold given as argument cannot be casted to float")
    sys.exit()

#Test de la variable CV_threshold
try:
    CV_threshold = float(CV_threshold)
except:
    print("CV_treshold given as argument cannot be casted to float")
    sys.exit()


#Conversion du fichier pdb en dicoPDB et initialisation
dicoPDB = PDBTools.parsePDBMultiChains(infile)
    
    
#Utilisation de la méthode CV pour extraire les résidus de surface dans un array
print('Détermination des résidus de surface ...\n')
surface_res_info = numpy.array([[],[],[],[],[],[],[]])#chain, res_name, coor_x, coor_y, coor_z, stikiness, cluster_num       

listXYZ_allatoms = CV.allatoms_coor_list(dicoPDB)

for chain in dicoPDB["chains"]:
    
    for residue in dicoPDB[chain]["reslist"]:
        
        #calcul des coordonnées XYZ du barycentre du résidu
        coor_bary = PDBTools.bary_res(dicoPDB[chain][residue])
        
        #format redondant spécifique à la fct CV_res()
        coor_bary_CV = [[coor_bary[0]], [coor_bary[1]], [coor_bary[2]]]

        tmp_CV = CV.CV_res(coor_bary_CV, listXYZ_allatoms, rc_threshold)
        
        #Si le résidu remplie la condition de CV_threshold, il est considérer comme de surface et ses informations sont stockés dans le numpy.array surface_res_info
        if(tmp_CV <= CV_threshold):
            surface_res_info = numpy.c_[surface_res_info, [chain, residue, coor_bary[0], coor_bary[1], coor_bary[2], 0, 0]]
        
        #Si le résidu ne rempliela condition de CV_threshold, il est considérer comme de coeur et on ne le prend pas en compte, on affecte 0 à sa valeur de bfactor
        else:
            for atom in dicoPDB[chain][residue]["atomlist"]:
                dicoPDB[chain][residue][atom]["bfactor"] = 0

#Calcul Stickiness pour les résidus de surface
print('Calcul de la Stickiness pour chaque résidu de surface ...\n')

min_stickiness = min(dicoSticky.values())
max_stickiness = max(dicoSticky.values())

max_color = 10
min_color = 100

for tmp_res_index in range(len(surface_res_info[0])):
    
    tmp_resname = dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]]["resname"]
    
    try:
        tmp_stickiness = dicoSticky[tmp_resname[0:3]]
        tmp_stickiness_normalised = (max_color - min_color)*(tmp_stickiness - min_stickiness) / (max_stickiness - min_stickiness) + min_color
        
        for atom in dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]]['atomlist']:
                dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]][atom]["bfactor"] = tmp_stickiness_normalised
        
    except:
        #Si on ne trouve pas de valeur de stickiness dans le dicoSticki, on ne conserve pas le résidu pour la suite du programme et on affecte 0 à sa valeur de bfactor
        for atom in dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]]['atomlist']:
            dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]][atom]["bfactor"] = 0

#Ecriture du fichier pdb output
print('Ecriture du fichier pdb d\'output ...\n')
PDBTools.writePDB(dicoPDB, outfile, True)