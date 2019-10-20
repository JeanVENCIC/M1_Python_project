# -*- coding: utf-8 -*-


import PDBTools as pdbt
import CVmaker as CVm


import sys
import os
import shutil



slash = '/'
if(os.name=='nt'):
    slash='\\'


"""
ce fichier doit etre en presence du fichier stickiness.txt
stickiness.txt: fichier de correspondance AA -> stickiness
"""

"""
cette fonction fait un dico associatif entre AA et stickiness à partir d'un fichier
    
param:
    file:: un string correspondant au fichier
    
return:
    un dico associant AA -> valeur de stickiness
"""
def setStickinessAA(file):
    
	#initialisation du dico
    resSticky = {}    
    
    f = open(file) #ouverture du fichier
    lines = f.readlines()       
    f.close   

	#on traite les lignes
    for line in lines: 
        
        tmp = line.strip("a\n")
        tmp = tmp.split("  ") #apparament dans le fichier, la séparation est faite par 3 " "
        resSticky[tmp[0]] = float(tmp[1])
    
    #cas ou le residu est inconnu
    resSticky["UNK"]=0.0

    return resSticky
    
dicoSticky= setStickinessAA("stickiness.txt")



"""
fonction permettant de déterminer les AA exposés. prend file, un string correspondant à un fichier
retourne un dico associant le nom du fichie à deux liste d'AA au même format que le parser du pdb.
    
param:
    file: un string correspondant au fichier pdb
    method: un string correspondant à la méthode à utiliser pour déterminer l'enfouissement des protéines
    crit: un float servant de critere de séparation entre exposé/ non exposé
    rc: distance seuil pour determiner les atomes environnants un atome i
    
return:
    un dico composé de deux liste, "expo" correspondant à la liste des residues exposés, et "nonexpo" pour le contraire
"""
def sepExpose(file, method="CV", crit=0.5, rc=20, dirOut="valeurCV"): #pour l'instant qu'une méthode, la deuxieme sera freesasa, crit correspond à un parametre de seuil, utile en particulier pour CV
    
    #initialisation du dico à retourner
    res={}
    res["expo"] = []
    res["nonexpo"] = []
        
    tmp=None
    
    fileName = (file.split(slash)[-1]).split(".")[0] + ".cv"
    path= dirOut + slash + fileName


    if( method == "CV" or method == "CVbarycentre"):

        #verification que le fichier existe
        if(os.path.isfile(path)):
                
            tmp= CVm.loadFileCV(path)
            
        #sinon on le créer sur le tas
        else:
            CVm.makeFileCV(file, path, method, rc)
            tmp= CVm.loadFileCV(path)
            
            
        for chain in tmp["chain"]:
            
            for resi in tmp[chain]["reslist"]:
            
                if (tmp[chain][resi]["CV"] < crit):
                    res["expo"].append(tmp[chain][resi])
                else:
                    res["nonexpo"].append(tmp[chain][resi])
    
    else:
        print("methode inconnue")
        sys.exit()
                    
    return res


"""
cette fonction permet de calculer la stickiness d'une proteine, elle prend en entré une liste d'AA, et le fichier de correspondance AA ->stickiness
    
param:
    liste: une liste de residus format dicodPDB
	
return:
    la valeur de la stickiness de la liste de residus
"""
def stickiness(liste, norm="Na"): #Na = pas de normallisation
    
    res=0
    tmpNorm=1; #initialisation aux cas pas de normalisation
    
    if(norm == "nbRes"):#normalisation par le nombre de résidu
        tmpNorm = len(liste)
  
    
    #on itere sur les residus
    for resi in liste:      
        
        #on incremente la valeur de stickiness
        if (resi["resname"] in dicoSticky):        
            res = res + dicoSticky[resi["resname"]]    
        else:
            print("AA inconnue " + resi["resname"])
    
    try:
        res = res / tmpNorm
    except:
        print("erreur div par zero\n")
        pass
        
        
    return res
        
    

"""
la fonction qui va permetre de lister les valeurs de stickiness de chaque fichier pdb
prend une liste de String correspondant aux nom des fichiers pdb
    
param:
    listFile: une liste de string correspondant aux pdb à traiter
    fileOut: parametre optionnel, servant à la création d'un fichier de sortie
    norm: un string correpondant au critere de normalisation
    method: un string correspondant à la méthode à utliser pour estimer l'enfouissement
    crit: un float servant de seuil pour determiner l'enfouissement
    rc: un entier servant de seuil de distance pour déerminer l'environnement d'un atome
        
return:
    un dico associant le nom du fichier à une liste composer de la stickiness est AA exposés, des non exposés, et de l'ensemble des 2
    
	"""
def listeStickiness(directory, fileOut=None, norm="Na", method="CV", crit=0.5, rc=20, dirCV="valeurCV"):
    
    
    #pathInitial = os.getcwd()   
    
    dicoSticki = {};    
    dicoSticki["listFile"] = []
    
	#declaration du fichier de sortie
    out=None    
    
    if (fileOut!=None): #petite amélioration, si deuxieme parametre, ecriture dans un ficheir de sortie, (peut-etre utile sur une couche superieur)
        out = open(fileOut, "w")
        out.write("#method:" + method + "\tnorm:" + norm + "\tcrit:" + str(crit) + "\trc:" + str(rc) + "\n")
       
    
    #on recupere les nons de fichier
    listFile = os.listdir(directory)
    
        
    try:
        os.mkdir(dirCV)
    except OSError:
        pass
        
    
    CVm.makeAllFileCV(directory, dirCV, method, rc)  
    
    #on parcours tout les fichiers
    for file in listFile: 
        
        print("je traite: " + file)
            
        #on recupere la liste des residues enfouies
        exposition = sepExpose(directory + slash + file, method, crit, rc, dirCV)        
        
        #on calcule la stickiness
        valExpo = stickiness(exposition["expo"], norm) # exposition doit etre un dico respectant le schema 
        valNonExpo = stickiness(exposition["nonexpo"], norm)
        valTot = stickiness(exposition["expo"] + exposition["nonexpo"], norm)
        
        #stockage dans une liste
        valeurSticky = [valExpo, valNonExpo, valTot]
            
        nameTmp=file.split(".")[0]
        
        #on ajoute les valeurs aux dico de sortie
        dicoSticki["listFile"].append(nameTmp)          
        dicoSticki[nameTmp] = valeurSticky
        
		#s'il existe un deuxieme parametre, on écrit dans un fichier de sortie
        if (fileOut!=None): 
            out.write(nameTmp + "\t" + str(valExpo) +"\t" + str(valNonExpo) +"\t" + str(valTot) + "\n")   
    
    
    #on ferme le fichier d'écriture
    if (fileOut!=None):
        out.close()
        
    
    
    return dicoSticki


"""
fonction servant à charger des valeurs de stickiness pour des protéines
param:
    file: un string correspondant à un fichier de correspondance entre proteine et sa stickiness
return:
    un dictionnaire de stickiness
"""
def loadFile(file):
    
    dicoSticki = {};    
    dicoSticki["listFile"] = []    
    
    try:    
    
        f = open(file, "r")
        lines = f.readlines()
        f.close()
        
        for line in lines :
            
            if(line[0]!='#'):        
            
                tmp = line.strip("\n")
                tmp = tmp.split("\t")        
                
                dicoSticki["listFile"].append(tmp[0])          
                dicoSticki[tmp[0]] = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
                
    except FileNotFoundError:
         print(file + " notFound\n")
        
    
    return dicoSticki







#fonction utilitaire

"""
fonction permettant d'afficher un dico de stickiness dans un terminal
param:
    dicoSticki: un dico de stickiness
"""
def afficheDicoSticky(dicoSticki):

    for a in dicoSticki["listFile"]:
        
        print("{}\t{}\t{}\t{}" .format(a,dicoSticki[a][0],dicoSticki[a][1],dicoSticki[a][2]))


"""
fonction permettan de dire si un fichier contient des elements non prévu par le dico de stickiness
param:
    file: un fichier pdb
return:
    un booleen correspondant au resultat
"""
def tcheckIfUnkAA(file):

        pdb = pdbt.parsePDBMultiChains(file)
        
        for a in pdb["chains"]:
            for b in pdb[a]["reslist"]:
                if not(pdb[a][b]["resname"] in dicoSticky):
                    print("probleme sur: " + file)
                    return 1
                
        return 0
        
        
"""
fonction permettant de recupérer les fichier comportant uniquement des leements reconnue par le dico
param:
    directory: le répertoire contenant les pdb
"""
def filtre(directory):
    
    print("lancement du filtre\n")
    
    pathInitial = os.getcwd()
    listTmp= []
    listFile=os.listdir(directory)
    os.chdir(directory)
    
    for file in listFile:
        
        if(tcheckIfUnkAA(file)==0):
            listTmp.append((directory + slash + file)) #on ajoute le chemin du dossierou est le fichier et le fichier
    
    print("\ncopie des fichier\n")
    os.chdir(pathInitial)    
    try:
        os.mkdir("filtreOut")
    except OSError:
        print("dossier deja existant ou erreur")
        pass
    
    for file in listTmp:
        
        try:
            print("copie de " + file)
            shutil.copy(file, "filtreOut")
        except :
                pass
        
    print("\nfin de la copie\n")




"""
fonction d'affichage de l'utilisation en ligne de commande
"""
def afficheUsage():
    print("\nUsage:")
    print("necessary argument:")
    print("\t-d <input PDB directory>")
    print("or")
    print("\t-f <input pdb file>")
    print("unnecessary argument:")
    print("\t-cv <input .cv directory>, default='filesCV'")  
    print("\t-o <output Stickiness file>, default='StickinessOut.txt'")    
    print("\t-m <method>, method to calculate CV, choose between CV and CVbarycentre default='CVbarycentre'")
    print("\t-n <norm>, method to normalise Stickiness, choose between Na and nbRes default='Na'")
    print("\t-dist <CV distance>, positive integer default=10")
    print("\t-s <CVcrit>, between 0 and 1 default=0.75\n")
    

###########################
#
#      MAIN
#
###########################
    
    
#verification si le script est lancer en ligne de commande ou non
if (sys.argv[0] == "stickiness.py") :
    
    if("-h" in sys.argv):
        afficheUsage()
    else:
        
        try :

             #recuperation des 2 elements minimaux
            try:
                dirPDB = sys.argv[sys.argv.index("-d")+1]
            except:
                pass
            
            try:
                fileIn = sys.argv[sys.argv.index("-f")+1]
            except:
                pass
            
            try:
                dirCV = sys.argv[sys.argv.index("-cv")+1]
            except:
                dirCV = "filesCV"
            
            try:
                fileOut = sys.argv[sys.argv.index("-o")+1]
            except:
                fileOut = "StickinessOut.txt"
                
            try:
                method = sys.argv[sys.argv.index("-m")+1]
            except:
                method = "CVbarycentre"
                
            try:
                norm = sys.argv[sys.argv.index("-n")+1]
            except:
                norm = "Na"
                
            try:
                crit = float(sys.argv[sys.argv.index("-c")+1])
            except:
                crit = 0.75
    
            try:
                rc = int(sys.argv[sys.argv.index("-dist")+1])
            except:
                rc = 10
                
                
            if ("-f" in sys.argv):
                    #on recupere la liste des residues enfouies
                    
                    exposition = sepExpose(fileIn, method, crit, rc, dirCV)        

                    #on calcule la stickiness
                    valExpo = stickiness(exposition["expo"], norm) # exposition doit etre un dico respectant le schema 
                    valNonExpo = stickiness(exposition["nonexpo"], norm)
                    valTot = stickiness(exposition["expo"] + exposition["nonexpo"], norm)
                    
                    #stockage dans une liste
                    print("stickiness expo = "  + str(valExpo))
                    print("stickiness nonexpo = "  + str(valNonExpo))
                    print("stickiness Tot = "  + str(valTot))
                    
            elif ("-d" in sys.argv):
                listeStickiness(dirPDB, fileOut=fileOut, norm=norm, method=method, rc=rc, dirCV=dirCV)
                
            else:
                afficheUsage()
                exit()
                
        except:
            afficheUsage()
            exit()
    
    