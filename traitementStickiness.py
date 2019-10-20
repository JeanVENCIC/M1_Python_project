# -*- coding: utf-8 -*-
"""
Created on Wed May 15 19:38:36 2019

@author: Pierre
"""

import stickiness as st

import time
import sys
import os

import matplotlib.pyplot as plt

import scipy.stats.mstats as ssm

import numpy as np

"""
fonction principale  gérant la création du fichier de stickiness
param:
        directory: le repertoire contenant les pdb
        dicoNom: un dico de dico cavec comme clé le nom de ficheir de sortie
        et comme champ:
            dico[nomFichier]["method"] #methode appliquer pour le calcul
            dico[nomFichier]["norm"]#methode de normalisation
            dico[nomFichier]["rc"] la distance de seuil pour CV
            dico[nomFichier]["crit"] le critere de séparation exterieur interieur de la proteine
            dico[nomFichier]["valCV"] dossier de valeur de CV
        
        remake: un boolean pour réutiliser le fichier stickiness s'il existe ou le recréer

return:
    un dico correspondants aux résultats des calculs
        
"""
def makeFiles(directory ,dicoNom ,remake=0):
    
    res = {}    
    
    for key in dicoNom:
        
        print("je traite le fichier : "  + key)
        
        start_time = time.time()
        
        try:
            file = open(key, "r");
            
            line = file.readline()
            
            file.close()
            
            line = line.strip("\n")
            line = line.strip("#")
            line = line.replace(":","\t")
            line = line.split("\t")
            
            #on verifie que les informations du fichiers sont conforme aux parametres
            if(len(line)!= 8 or line[1]!=dicoNom[key]["method"] or line[3]!=dicoNom[key]["norm"] or line[5]!=str(dicoNom[key]["crit"]) or line[7]!=str(dicoNom[key]["rc"])):
                print(key + " non conforme")
                print('creation de ' + key + "\n")
                res[key] = st.listeStickiness(directory ,fileOut=key, norm=dicoNom[key]["norm"], method=dicoNom[key]["method"], crit=dicoNom[key]["crit"], rc=dicoNom[key]["rc"], dirCV=dicoNom[key]["valCV"])
            else:
                print("fichier conforme")
                if(remake==1):
                    print('creation de ' + key + "\n")
                    res[key] = st.listeStickiness(directory ,fileOut=key, norm=dicoNom[key]["norm"], method=dicoNom[key]["method"], crit=dicoNom[key]["crit"], rc=dicoNom[key]["rc"], dirCV=dicoNom[key]["valCV"])
                else :
                    res[key] = st.loadFile(key)
            
        except FileNotFoundError:
            print(key + " not found")
            print('creation de ' + key + "\n")
            res[key] = st.listeStickiness(directory ,fileOut=key, norm=dicoNom[key]["norm"], method=dicoNom[key]["method"], crit=dicoNom[key]["crit"], rc=dicoNom[key]["rc"], dirCV=dicoNom[key]["valCV"])
            
            
        print("Temps d execution: %s secondes --- \n" % (time.time() - start_time))
            
    return res



"""
abondance: un dico associant pour un fichier pdb, son nom de gene, et son abondance
stich: un dico de stickiness pour une méthode de calcule
"""
def assoGeneAbonSticky(abondance, stick):    
    
    res = {}
    
    for key in abondance:
        if (key in stick):# and abondance[key]["abondance"]<100):
            res[key]={}
            res[key]["stickiness"] = stick[key]
            res[key]["gene"] = abondance[key]["gene"]
            res[key]["abondance"] = abondance[key]["abondance"]
            
            
    return res #on doit se retrouver avec un dico associant pour un fichier pdb son gene, son abondance, et sa stichiness
           
           
"""
fonction de chargement du fichier de correspondance pdb abondance
param:
    file: un fichier de correspondance: #gene_OLN	PDB_ID	abundance
return:
    un dico avec en clé le PDB_ID et 2 champs de vlaeurs: gene et abondance
"""
def loadAbon(file):

    res={}
    
    try:
        
        f= open(file, "r")
        lines = f.readlines()
        f.close()
        
    except:
        print("erreur : " + file + " not found")
        exit()
        
    for line in lines:
            
            if(line[0]!='#'):            
            
                line = line.split("\t")
                res[line[1]] = {}
                res[line[1]]["gene"] = line[0]
                res[line[1]]["abondance"] = float(line[2])
    
    return res


"""
creation de 4 listes pour faire des plots
param:
    dico: un dico de correspondance abondance/Stickiness
return:
    une liste de liste correspondant à: abonbdance / stickiness surface / stickiness interne / stickiness total
"""
def makeListe(dico):
    res = [[],[],[],[]]
    
    for key in dico:
        res[0].append(dico[key]["abondance"])
        res[1].append(dico[key]["stickiness"][0])
        res[2].append(dico[key]["stickiness"][1])
        res[3].append(dico[key]["stickiness"][2])
    
    return res



"""
fonction d affichage de 3 graph pour les 3 aa enfoui, de surface, et toto
param:
    listeAbondance: une liste correspondant aux valeurs d'abondance
    liste: une liste contennant 3liste de coordonnées x de stickiness
    title: le titre du fichier de sortie
"""
def plot3(listeAbondance, liste, title):
    
    
    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].scatter(liste[0], listeAbondance)
    axarr[0].set_title('Stickiness expose') 
    axarr[0].set_ylabel("Abondance")
    
    axarr[1].scatter(liste[1], listeAbondance)
    axarr[1].set_title('Stickiness enfoui') 
    axarr[1].set_ylabel("Abondance")
    
    axarr[2].scatter(liste[2], listeAbondance)
    axarr[2].set_title('Stickiness totale')
    axarr[2].set_ylabel("Abondance")    
    
    plt.xlabel("Stickiness")
    
    plt.savefig(title)
    
"""
fonction d'affichage d'histogramme de la stickiness
param:
    liste: une liste contennant 3liste de coordonnées x
    title: le titre du fichier de sortie
"""
def hist3(liste, title):
    
    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].hist(liste[0], bins = 20)
    axarr[0].set_title('Stickiness expose')
    axarr[1].hist(liste[1], bins = 20)
    axarr[1].set_title('Stickiness enfoui')
    axarr[2].hist(liste[2], bins = 20)
    axarr[2].set_title('Stickiness totale')
    
    plt.xlabel("Stickiness")
    
    plt.savefig(title)    

"""
fonction d'affichage des informations
param:
        listex, une liste correspondant aux coordonnées sur l'axe x
        listey, une lsite correspondant aux coordonnées sur l'axe y
        title, le titre du graphique
"""
def plotInfo(listex,listey, title):
    
    x=np.array(listex)
    y=np.array(listey)

    #calcul des valeurs
    n = len(x)
    mean = np.mean(x)
    var = np.var(x)
    
    
    f = plt.figure()
    ax = f.add_subplot(111)    
    
    quantiles=ssm.mquantiles(x)

    
    leTexte = (
    "nbPoint: " + '%.0f' % n +"\n"    
    + "mean: " + '%.4f' % mean +"\n"
    + "var: " + '%.4f' % var + "\n"
    + "quantile1: " + '%.4f' % quantiles[0] + "\n"
    + "quantile2: " + '%.4f' % quantiles[2]
    )  
    
    #affichage des infos
    plt.text(0.82,0.80,leTexte,horizontalalignment='center',
     verticalalignment='center', transform = ax.transAxes)

    #affichage du graph
    plt.scatter(x,y, s = 7)
    
    #affichage des quantiles
    plt.axvline(x=quantiles[0], linewidth=3, color='g')    
    plt.axvline(x=quantiles[2], linewidth=3, color='g')  
    
    #labels
    plt.xlabel("Stickiness")
    plt.ylabel("Abondance")
    plt.title(title)
    
    plt.savefig(title)    
    

"""
fonction d'affichage d'un dico d'association abondance stickiness
"""
def afficheAsso(dico):
    
    f= open("debug.txt", "w")    
    
    for key in dico:
        f.write(key +"\t"+ str(dico[key]["abondance"]) + "\t" + str(dico[key]["stickiness"][0])+ "\t" + str(dico[key]["stickiness"][1])+ "\t" + str(dico[key]["stickiness"][2]) +"\n")
    
    f.close()

"""
fonction d'affichage de l'utilisation en ligne de commande
"""
def afficheUsage():
    print("\nUsage:")
    print("necessary argument:")
    print("\t-inPDB <input PDB directory>")
    print("\t-a <input abondanceCorrelation file> (#gene_OLN	PDB_ID	abundance)")
    print("unnecessary argument:")
    print("\t-inCV <input .cv directory>, default='filesCV'")  
    print("\t-oSt <output Stickiness file>, default='StickinessOut.txt'")    
    print("\t-oGr <output graph directory>, default='graphOut'")
    print("\t-d <CV distance>, positive integer default=10")
    print("\t-s <CVcrit>, between 0 and 1 default=0.75")
    print("\t-m <method>, method to calculate CV, choose between CV and CVbarycentre default='CVbarycentre'")
    print("\t-n <norm>, method to normalise Stickiness, choose between Na and nbRes default='Na'")
    print("\t-r <remake>, choose 1 to rebuilt StickinessFile, 0 otherwize default=1\n")
    
###########################
#
#      MAIN
#
###########################


#verification si le script est lancer en ligne de commande ou non
if (sys.argv[0] == "traitementStickiness.py") :
    
    if("-h" in sys.argv):
        afficheUsage()
    else:
        
        try :

            #recuperation des 2 elements minimaux
            try:
                indirPDB = sys.argv[sys.argv.index("-inPDB")+1]
            except:
                afficheUsage()
                sys.exit()
            
            try:
                fileAbon = sys.argv[sys.argv.index("-a")+1]
            except:
                afficheUsage()
                sys.exit()
            
            
            try:
                indirCV = sys.argv[sys.argv.index("-inCV")+1]
            except:
                indirCV = "filesCV"
            
            try:
                nomFichier = sys.argv[sys.argv.index("-oSt")+1]
            except:
                nomFichier = "StickinessOut.txt"
                
            try:
                graphOut = sys.argv[sys.argv.index("-oGr")+1]
            except:
                graphOut = "graphOut"
                
            try:
                distSeuil = int(sys.argv[sys.argv.index("-d")+1])
            except:
                distSeuil=10
                
            try:
                CVseuil = float(sys.argv[sys.argv.index("-s")+1])
            except:
                CVseuil = 0.75
                
            try:
                method = sys.argv[sys.argv.index("-m")+1]
            except:
                method = "CVbarycentre"
                
            try:
                norm = sys.argv[sys.argv.index("-n")+1]
            except:
                norm = "Na"
                
            try:
                rem = int(sys.argv[sys.argv.index("-r")+1])
            except:
                rem = 1
                
            
            #dico des différentes combinaisons
            dico = {} #pour empiler les instructions permet de faire une serie d'instruction à appliquer
            #cle = non du fichier ou on a l'association fichier stickiness
            dico[nomFichier] = {}
            dico[nomFichier]["method"] = method #methode appliquer pour le calcul
            dico[nomFichier]["norm"] = norm #methode de normalisation
            dico[nomFichier]["rc"] = distSeuil
            dico[nomFichier]["crit"] = CVseuil
            dico[nomFichier]["valCV"] = indirCV
            
            #clacul de la stickiness
            laStickiness = makeFiles(indirPDB, dico, remake=rem) #dico de chaque stickyness
            
            #chargement des vlaeurs d'abondance
            dicoValeur = loadAbon(fileAbon)
            
            #association de l'abondance a la stickiness
            dicoValeur = assoGeneAbonSticky(dicoValeur, laStickiness[nomFichier])
            
            #creation d'une lsite correspondant aux résultat
            liste = makeListe(dicoValeur)
            
            
            liste1= liste[0]
            liste2= [liste[1],liste[2],liste[3]]    
            
            
            #creation du dossier de sorti si inexistant
            try:
                os.mkdir(graphOut)
            except OSError:
                #print("dossier deja existant ou erreur")
                pass
            
            #recuperation du path initial
            pathInitial= os.getcwd()
            
            os.chdir(graphOut)
            
            #affichage des graphs
            plot3(liste1, liste2, "AbonEtSticki")
            hist3(liste2, "repartitionSticki")
            plotInfo(liste2[0], liste1, "infoStickiSurface")
            plotInfo(liste2[1], liste1, "infoStickiEnfoui")
            plotInfo(liste2[2], liste1, "infoStickiTotal")
            
            os.chdir(pathInitial)
        
        except:
            
            afficheUsage()




