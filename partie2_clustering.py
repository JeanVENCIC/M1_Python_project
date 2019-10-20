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
    print("\nUsage :\n\nCompact call :\n", sys.argv[0], "<input PDB file> <output PDB file> <file containing script arguments>")
    print("\nComplete call :\n", sys.argv[0], "\n\t<input PDB file>\n\t<output PDB file>\n\t<stickiness / AA correlation file>\n\t<rc threshold used to compute CV value>\n\t<CV threshold to determine surface residues>\n\n\tClustering method parameters :\n\t<significative distance threshold for clustering> (in Angstrom)\n\t<significative stickiness threshold for clustering>\n\t<maximum number of cluster>\n\t<maximum number of residues by cluster>\n")


###########################
#
#      MAIN
#
###########################
if("-h" in sys.argv):
    printUsage()
    
#Compact call management
elif(len(sys.argv) == 4):

    argfile = sys.argv[3]

    try:
        f = open(argfile, "r")
    except:
        print('Unable to read argfile :', argfile)
        sys.exit()
    
    lines = f.readlines()
    f.close()

    #recherche du nombre de lignes de commentaire avant les arguments
    cpt_comlines = 0
    for line in lines:
        if (line[0] == "#"):
            cpt_comlines += 1
        else:
            break
    #extraction des arguments
    stickiness_file = lines[0 + cpt_comlines].strip().split("=")[1].strip()
    rc_threshold = lines[1 + cpt_comlines].strip().split("=")[1].strip()
    CV_threshold = lines[2 + cpt_comlines].strip().split("=")[1].strip()
    dist_threshold = lines[3 + cpt_comlines].strip().split("=")[1].strip()
    sticky_threshold = lines[4 + cpt_comlines].strip().split("=")[1].strip()
    nb_max_cluster = lines[5 + cpt_comlines].strip().split("=")[1].strip()
    nb_max_res = lines[6 + cpt_comlines].strip().split("=")[1].strip()   

#Complete call management
elif(len(sys.argv) == 10):
    rc_threshold = sys.argv[3]
    stickiness_file = sys.argv[4]
    CV_threshold = sys.argv[5]
    dist_threshold = sys.argv[6]
    sticky_threshold = sys.argv[7]
    nb_max_cluster = sys.argv[8]
    nb_max_res = sys.argv[9]

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

#Test du fivhier stickiness.txt
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

#Test de la variable dist_threshold
try:
    dist_threshold = float(dist_threshold)
except:
    print("dist_treshold given as argument cannot be casted to float")
    sys.exit()

#Test de la variable sticky_threshold
try:
    sticky_threshold = float(sticky_threshold)
except:
    print("sticky_treshold given as argument cannot be casted to float")
    sys.exit()

#Test de la variable nb_max_cluster
try:
    nb_max_cluster = int(nb_max_cluster)
except:
    print("nb_max_cluster given as argument cannot be casted to int")
    sys.exit()
    
if(nb_max_cluster < 0):
    print("nb_max_cluster must be > 0")
    sys.exit()

#Test de la variable nb_max_res
try:
    nb_max_res = int(nb_max_res)
except:
    print("nb_max_res given as argument cannot be casted to int")
    sys.exit()

if(nb_max_res <= 0):
    print("cluster_nb_res_min must be >= 1")
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
index_to_keep = []

for tmp_res_index in range(len(surface_res_info[0])):
    
    tmp_resname = dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]]["resname"]
    
    try:
        tmp_stickiness = dicoSticky[tmp_resname[0:3]]

        if(tmp_stickiness >= sticky_threshold):
            surface_res_info[5,tmp_res_index] = dicoSticky[tmp_resname[0:3]]
        
            index_to_keep.append(tmp_res_index)
        else:
            #Si la valeur de stickiness ne remplie pas la condition de stickiness_threshold, on ne conserve pas le résidu pour la suite du programme et on affecte 0 à sa valeur de bfactor
            for atom in dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]]['atomlist']:
                dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]][atom]["bfactor"] = 0
        
    except:
        #Si on ne trouve pas de valeur de stickiness dans le dicoSticki, on ne conserve pas le résidu pour la suite du programme et on affecte 0 à sa valeur de bfactor
        for atom in dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]]['atomlist']:
            dicoPDB[surface_res_info[0, tmp_res_index]][surface_res_info[1, tmp_res_index]][atom]["bfactor"] = 0

#On ne garde du numpy.array que la sous-partie des résidus pour lesquels ont a pu calculer une valeur de stickiness
surface_res_info = surface_res_info[:,index_to_keep]

#Clustering des régions de résidus de surface à haute stickiness
print('Clustering des résidus proches à haute stickiness ...\n')
nb_cluster = 0

#liste très importante qui contient les indices de résidus du tableau surface_res_info qui n'ont pas encre été intégrés à un cluster
res_index_notclustered = list(range(len(surface_res_info[0])))

while (nb_cluster < nb_max_cluster and len(res_index_notclustered) !=0 ):
    
    nb_cluster += 1
    print('\tCréation du cluster n°',nb_cluster,'\n')
    tmp_nb_res = 1#le centre du cluster
    
    #détermination du centre du nouveau cluster
    max_sticki_nonclustered = max(surface_res_info[5, res_index_notclustered])
    
    tmp_center_index = res_index_notclustered[list(surface_res_info[5, res_index_notclustered]).index(max_sticki_nonclustered)]
    tmp_center_coorlist = list(surface_res_info[2:5,tmp_center_index].astype(float))
    
    #on affecte au centre du cluster le numéro du cluster
    surface_res_info[6, tmp_center_index] = nb_cluster
    
    #on enlève le centre de la liste des résidus non-clusterisés
    res_index_notclustered.remove(tmp_center_index)


    #Calcul de la distance de chaque résidu non cluster-isé au nouveau centre    
    
    #liste des distances au centre des résidus non-clusterisés
    tmp_dist_list = []
    
    for tmp_res_index in res_index_notclustered:

        
        tmp_res_coorlist = list(surface_res_info[2:5, tmp_res_index].astype(float))
        
        tmp_dist_list.append(CV.Norme(CV.Rij(tmp_center_coorlist, tmp_res_coorlist)))
    
    #remplissage du nouveau cluster
    
    #trie de la liste des distances au centre du cluster (on garde les indices de la liste tmp_dist_list triée pour pouvoir remonter aux indices des résidus non-clusterisés)
    sorted_dist_list_index = sorted(range(len(tmp_dist_list)), key=lambda k: tmp_dist_list[k])
    
    i = 0
    
    #indices des résidus que l'on va clusteriser
    res_index_clustered = []
    
    while(tmp_nb_res < nb_max_res and i < len(sorted_dist_list_index)):
        
        #parcours de la liste triée du plus proche au moins proche du centre et test de l'ajout du résidu au cluster
        tmp_dist_index = sorted_dist_list_index[i]
        
        if(tmp_dist_list[tmp_dist_index] <= dist_threshold):
    
            surface_res_info[6, res_index_notclustered[tmp_dist_index]] = nb_cluster
            
            res_index_clustered.append(res_index_notclustered[tmp_dist_index])
            
            tmp_nb_res += 1
            
            if(tmp_nb_res == nb_max_res):
                break
        else:
            pass
        
        i += 1
    
    #on retire tous les résidus clusterisés de la liste des non-clusterisés    
    for tmp_res_index_clustered in res_index_clustered:
        res_index_notclustered.remove(tmp_res_index_clustered)
    

#Reprise du dicoPDB avec le numéro de cluster en bfactor
for row in range(len(surface_res_info[0])):
    
    tmp_chain = surface_res_info[0, row]
    tmp_res = surface_res_info[1, row]
    tmp_cluster = float(surface_res_info[6, row])
    
    for atom in dicoPDB[tmp_chain][tmp_res]["atomlist"]:
        dicoPDB[tmp_chain][tmp_res][atom]['bfactor'] = tmp_cluster

#Ecriture d'un fichier de sauvegarde pour le tableau surface_res_info
print('Ecriture du fichier de sauvegarde \"surface_res_info.txt\" ...\n')

f = open("surface_res_info.txt","w").close()
f = open("surface_res_info.txt","a")

for row in range(len(surface_res_info[0])):
    
    line_to_write = surface_res_info[0, row]
    
    for element in range(len(surface_res_info))[1:]:
 
        line_to_write = line_to_write + '\t' + surface_res_info[element][row]
    
    f.write(line_to_write+'\n')    

f.close()

#Ecriture du fichier pdb output
print('Ecriture du fichier pdb d\'output ...\n')
PDBTools.writePDB(dicoPDB, outfile, True)