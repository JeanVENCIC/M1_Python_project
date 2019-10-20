# -*- coding: utf-8 -*-
"""
Created on Thu May 23 11:00:06 2019

@author: Pierre
"""

"""
fonction permettant de gérer la parti CV
"""

import CV
import os
import sys
import PDBTools as pdbt
import time


slash = '/'
if(os.name=='nt'):
    slash='\\'


"""
fonction de chargement de la liste d'atome d'un pdb
param:
    pdb: un dico pdb
return:
    une liste de de vecteur [x,y,z]
"""
def loadAtomeList(pdb):
    
    #initialisation du tableau d'atome pour CV
    atomList = []
    
    #on parcours les chaines
    for chain in pdb["chains"]:
        
        #on parcours les résidies
        for resi in pdb[chain]["reslist"]:
            
            #on parcours les atomes
            for atome in pdb[chain][resi]["atomlist"]:
            
                #on ajoute tout les atomes à la liste d'atome
                atomList.append([pdb[chain][resi][atome]["x"], pdb[chain][resi][atome]["y"], pdb[chain][resi][atome]["z"]])
                
    return atomList
    
"""
fonction de création d'un fichier .cv pour un fichier pdb donné
param:
    filePDB: le fichier PDB d'entré
    fileCV: le fichier cv de sortie
    method: la methode a appliquer pour le calcul de CV
    rc: la distance que CV doit utiliser
"""
def makeFileCV(filePDB, fileCV, method="CVbarycentre", rc=20):
    
    print("creation du fichier: " + fileCV)
    
    pdb= pdbt.parsePDBMultiChains(filePDB)
    
    out = open(fileCV, "w")
    out.write("#file:" + (filePDB.split(slash)[-1]).split(".")[0] +"\t" "method:" + method + "\trc:" + str(rc) + "\n")
    
    if( method == "CV"):
        
        #initialisation du tableau d'atome pour CV
        atomList = loadAtomeList(pdb)   
            
        #maintenant on a notre de liste d'atome pour CV
        #on réitere sur chaque residu pour calculer le CV
        
        #on parcours les chaines
        for chain in pdb["chains"]:
            
            #on parcours les résidues
            for resi in pdb[chain]["reslist"]:
                
                #initialisation de la structure à faire passé à CV
                atomRes = []
                
                #on parcours les atomes
                for atome in pdb[chain][resi]["atomlist"]:
                    
                    atomRes.append([pdb[chain][resi][atome]["x"], pdb[chain][resi][atome]["y"], pdb[chain][resi][atome]["z"]])
                    
                #calcul de CV
                tmp = CV.CV_res2(atomRes, atomList, rc)
                
                #ecriture du resultat
                out.write(chain + "\t" + resi + "\t" + pdb[chain][resi]["resname"] + "\t" + str(tmp) + "\n")
    
    #cas ou on utilise CV sur le barycentre de l'AA
    elif( method == "CVbarycentre"):
        
        #initialisation du tableau d'atome pour CV
        atomList =loadAtomeList(pdb)   
            
        #maintenant on a notre de liste d'atome pour CV
        #on réitere sur chaque residu pour calculer le CV
        
        #on parcours les chaines
        for chain in pdb["chains"]:
            
            #on parcours les résidues
            for resi in pdb[chain]["reslist"]:
                
                #initialisation de la structure à faire passé à CV
                atomRes = [pdbt.getCentroid(pdb[chain][resi])]
                
                #calcul de CV
                tmp = CV.CV_res2(atomRes, atomList, rc) #deuxieme parametre à null car un parametre inutile dans CV_res

                #ecriture du resultat
                out.write(chain + "\t" + resi + "\t" + pdb[chain][resi]["resname"] + "\t" + str(tmp) + "\n")
                
    out.close()

"""
fonction de chargement d'un fichier .cv
param:
    file: le fichier .cv
return:
    un dico semblable au format pdb mais sans le champ atome, et avec un champ CV suplémentaire
"""
def loadFileCV(file):
    
    #initialisation du dico
    res = {}
    res["chain"] = []
    
    try:    
    
        f = open(file, "r")
        lines = f.readlines()
        f.close()
        
        #on récup les informations
        for line in lines :
            
            if(line[0]!='#'):  
               
               tmp = line.strip("\n")
               tmp = tmp.split("\t")  
               
               if (not(tmp[0] in res["chain"])):
                   res["chain"].append(tmp[0])
                   res[tmp[0]]={}
                   res[tmp[0]]["reslist"]=[]
                   
               res[tmp[0]]["reslist"].append(tmp[1])
               
               residu = {}
               residu["resname"] = tmp[2]
               residu["CV"] = float(tmp[3])
               
               res[tmp[0]][tmp[1]] = residu
                
    except FileNotFoundError:
         print(file + " notFound\n")
         
    return res




"""
fonction servant a créer tout les fichier cv d'un répertoire, ne recréer pas les fichiers déjà existant
param:
    dirIn: le répertoire contennant les .pdb a traiter
    dirOut: le répertoire de sorti des .cv
    method: la methode a appliquer pour le calcul de CV
    rc: la distance que CV doit utiliser
"""
def makeAllFileCV(dirIn, dirOut, method="CVbarycentre", rc=20):
    
    #initialisation de la variable de temp
    start_time = time.time()    
    
    #recupération de la liste des fichier
    filesIn = os.listdir(dirIn)
    
    try:
        os.mkdir(dirOut)
    except OSError:
        pass
        
    
    filesOut = os.listdir(dirOut)
    
    for file in filesIn:
        
        #on applique la correspodnance en .cv de sortie
        fileName = (file.split(slash)[-1]).split(".")[0] + ".cv"
        
        #on verifie s'il n'est pas déjà présent
        if not(fileName in filesOut): #s'il n'est pas dnas le repertoir de sortie
            
            #on créer le .cv
            makeFileCV(dirIn + slash + file, dirOut + slash + fileName, method, rc)
    
    #affichage du temp qu'a pris l'opération
    print("Temps d execution makeAllFileCV: %s secondes --- \n" % (time.time() - start_time))
    
    
"""
fonction permettant de créer un pdb avec le bfactor permettant de déifférencier la surface
param:
        pdb: un dico pdb
        valeurCV: un dico de valeur de CV
        seuil: le seuil a appliquer pour différencier les AA
"""
def makePDB(pdb, valeurCV, pdbOut, seuil):
    
    for chain in pdb["chains"]:
            
            #on parcours les résidues
            for resi in pdb[chain]["reslist"]:
                
                for atom in pdb[chain][resi]["atomlist"]:
                    
                    if valeurCV[chain][resi]["CV"]<seuil:
                        
                        pdb[chain][resi][atom]["bfactor"] = 1
                    
                    else:
                        pdb[chain][resi][atom]["bfactor"] = 0
                        
                
                
    pdbt.writePDB(pdb, pdbOut ,bfactor=True)
    
    
    
    
"""
fonction d'affichage de l'utilisation en ligne de commande
"""
def afficheUsage():
    print("\nUsage:")
    print("\t-d <input .pdb directory> <output .cv directory>")  
    print("\t-f <input .pdb file> <output .cv file>")    
    print("\t-m <method>, default='CVbarycentre'")
    print("\t-dist <distance>, positive integer default=10")
    print("\t-makePDB <inPDBfile><inCVfile><outPDBfile><seuil>, to make a pdb output")
    print("\tif there are -d, -f, and -makePDB, -f and makePDB will be ignored\n")

    
###########################
#
#      MAIN
#
###########################


#verification si le script est lancer en ligne de commande ou non
if (sys.argv[0] == "CVmaker.py") :
    
    if("-h" in sys.argv):
        afficheUsage()
    else:
        
        try :

            try:
                indirPDB = sys.argv[sys.argv.index("-d")+1]
                outdirCV = sys.argv[sys.argv.index("-d")+2]
            except:
                pass
                
            try:
                infilePDB = sys.argv[sys.argv.index("-f")+1]
                outfileCV = sys.argv[sys.argv.index("-f")+2]
            except:
                pass
                
            try:
                inPDB = sys.argv[sys.argv.index("-makePDB")+1]
                inCV = sys.argv[sys.argv.index("-makePDB")+2]
                outPDB = sys.argv[sys.argv.index("-makePDB")+3]
                seuil = float(sys.argv[sys.argv.index("-makePDB")+4])
            except:
                pass
            
            try:
                method = sys.argv[sys.argv.index("-m")+1]
            except:
                method = "CVbarycentre"
                
            try:
                rc = int(sys.argv[sys.argv.index("-dist")+1])
            except:
                rc = 10
                
                
            if ("-d" in sys.argv):
                makeAllFileCV(indirPDB, outdirCV, method=method, rc=rc)
            elif("-f" in sys.argv):
                makeFileCV(infilePDB, outfileCV, method=method, rc=rc)
            elif("-makePDB" in sys.argv):
                
                tmp = loadFileCV(inCV)
                tmp2 = pdbt.parsePDBMultiChains(inPDB)
                
                makePDB(tmp2, tmp, outPDB, seuil)
            else:
                afficheUsage()
                exit()
            
        except:
            afficheUsage()
            exit()
                     
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                