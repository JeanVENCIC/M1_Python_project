#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys

###########################
#
#      fonctions
#
###########################

def parseUniprotOrganism(infile, flag_resolution) :
    """
    propos : parser un fichier de de type proteome provenant d'Uniprot et
            extrait pour chaque protéine poosédant un fichier PDB, son Ordered
            Locus Name (OLN), l'ID du fichier PDB ayant la meilleur
            résolution et optionnellement la valeur de résolution.
            
    input : un fichier protéome Uniprot, code réalisé pour le fichier
            "/home/jean/Bureau/UE_Python/Projet/uniprot_yeast_filtered_
            organism3A22Saccharomycescerevisiae28strainATCC204__.txt"
            
    output : Retourne un dico avec l'OLN en clé et l'ID_PDB (résolution) en
            valeur ainsi qu'une clé "OLNlist" associée à la liste des OLN
            compris dans le dico. Créé également un fichier txt contenant les
            informations du dico si demandé en argument.
    """
    #ouverture du fichier ligne par ligne
    f = open(infile, "r")
    lines = f.readlines()
    f.close()
    
    #initialisation du dico output
    dico = {}
    dico["OLNlist"] = []
    
    #initialisation de la variable stockant l'OLN de la prot en cours d'étude
    tmpOLN = ""
    
    #initialisation du dico temporaire stockant les ID_PDB avec leur
    #résolution associée si ils existent pour la prot étudiée
    tmpPDBfiles = {}
    tmpPDBfiles["ID"] = []
            
    for line in lines :
        
        #ID de la ligne contenant l'OLN
        if line[0:2] == "GN" :
            
            #découpage de la ligne
            for info in line.split(";") :
                info = info.strip()
                
                #récupération de l'OLN
                if "OrderedLocusNames=" in info :
                    tmpOLN = info[info.find("=")+1:]
                    
                    #condition spéciale pour supprimer d'éventuelles info
                    #accolées à l'OLN
                    if " " in tmpOLN :
                        tmpOLN = tmpOLN[:tmpOLN.find(" ")]
                    if "\t" in tmpOLN :
                        tmpOLN = tmpOLN[:tmpOLN.find("\t")]
        
        #ID de la ligne contenant les infos du(des) fichier(s) PDB           
        elif line[0:2] == "DR"  :
            if line[5:9] == "PDB;" :
                
                tmpPDB_ID = line.split(";")[1].strip()
                
                #si l'ID PDB n'a pas déjà été retenu pour l'OLN étudiée
                if not tmpPDB_ID in tmpPDBfiles["ID"] :
                    
                    #récupération de sa résolution
                    tmpPDB_res = line.split(";")[3].strip("A").strip()
                    
                    #si la résolution n'est pas indiquée une valeur INT_MAX lui
                    #est associée
                    if tmpPDB_res == "-" :
                        tmpPDB_res = sys.maxsize
                    
                    #on associe l'ID_PDB à sa valeur de résolution
                    tmpPDBfiles["ID"].append(tmpPDB_ID)
                    tmpPDBfiles[tmpPDB_ID] = float(tmpPDB_res)
        
        #ID de la ligne indiquant la fin d'une protéine
        elif line [0:2] == "//" :
            
            #si au moins un fichier PDB à été retenu pour la prot étudiée
            if len(tmpPDBfiles["ID"]) != 0 :
                
                #on initialise l'ID_PDB ayant la plus basse résolution avec le
                #premier PDB du dico tmpPDBfiles
                res_min_ID = tmpPDBfiles["ID"][0]
                res_min_val = tmpPDBfiles[res_min_ID]
                
                #on garde l'ID_PDB ayant la résolution la plus basse ainsi que
                #sa valeur
                for PDB_ID in tmpPDBfiles["ID"] :
                    if tmpPDBfiles[PDB_ID] < res_min_val:
                        res_min_ID = PDB_ID
                        res_min_val = tmpPDBfiles[PDB_ID]
                    
                #on ajoute au dico output l'OLN à OLNlist et associé à l'ID_PDB 
                if not tmpOLN in dico["OLNlist"] :
                    dico["OLNlist"].append(tmpOLN)
                    
                    if(flag_resolution == 1) :
                        dico[tmpOLN] = [res_min_ID, res_min_val]
                        
                    elif(flag_resolution == 0) :
                        dico[tmpOLN] = res_min_ID
            
            #on ré-initialise le dico tmpPDBfiles et l'ONL de la prot étudiée
            tmpOLN = ""
            tmpPDBfiles = {}
            tmpPDBfiles["ID"] = []
    
    
    
    return dico



def parseSauvegarde(infile):
    """
        propos : parser un fichier de de type proteome provenant d'Uniprot
                et extrait pour chaque protéine poosédant un fichier PDB, son
                Ordered Locus Name (OLN) ainsi que l'ID du fichier PDB ayant la
                meilleur résolution.
                
        input : un fichier txt comportant sur chaque couple OLN - ID_PDB - 
                (optionnellement résolution) : fichier créé lors du parsage du
                fichier précédent comme sauvegarde).
                
        outpud : Retourne un dico avec l'OLN en clé et l'ID_PDB (résolution) en
                valeur ainsi qu'une clé "OLNlist" associée à la liste des OLN
                compris dans le dico.
    """
    
    #ouverture du fichier ligne par ligne
    f = open(infile, "r")
    lines = f.readlines()
    f.close()
    
    #initialisation du dico output
    dico = {}
    dico["OLNlist"] = []
    
    #si le fichier contient le couple OLN - ID_PDB
    if (len(lines[0].split("\t")) == 2) :
    
        #ajout de l'ONL à OLN list et association de l'OLN avec son ID_PDB
        for line in lines :
            OLN_PDB = line.split("\t")
            dico["OLNlist"].append(OLN_PDB[0])
            dico[OLN_PDB[0]] = OLN_PDB[1]
    
    #si le fichier contient le couple OLN - [ID_PDB, resolution] 
    elif (len(lines[0].split("\t")) == 3) :
        #ajout de l'ONL à OLN list et association de l'OLN avec son ID_PDB et
        #la résolution associé en angstrom
        for line in lines :
            OLN_PDB = line.split("\t")
            dico["OLNlist"].append(OLN_PDB[0])
            dico[OLN_PDB[0]] = [OLN_PDB[1], OLN_PDB[2]]
            
    return dico

###########################
#
#      MAIN
#
###########################

if ((len(sys.argv) != 3) or ("-h" in sys.argv)) :
    print("\nUsage :", sys.argv[0], "\n\t<infile>\n\t<Y/N> (wether to keep resolution values)\n")
    
else :
    infile = sys.argv[1]
    try:
        f = open(infile, "r")
        f.close()
    except:
        print('Unable to read file :', infile)
        sys.exit()
    
    flag_resolution = sys.argv[2]
    
    if(flag_resolution in ["YES", "yes", 'y', 'Y']) :
        dico = parseUniprotOrganism(infile, 1)
    
        for tmpOLN in dico["OLNlist"]:
            print(tmpOLN, dico[tmpOLN][0], dico[tmpOLN][1],sep = "\t")
            
    elif(flag_resolution in ["NO", "no", "N", "n"]) :
        dico = parseUniprotOrganism(infile, 0)
    
        for tmpOLN in dico["OLNlist"]:
            print(tmpOLN, dico[tmpOLN],sep = "\t")
            
    else :
        print("argument invalide : flag_resolution différent de 'Y' ou 'N'")