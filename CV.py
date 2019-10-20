# -*- coding: utf-8 -*-



import math, numpy
from statistics import mean

def Rij(list1, list2):
    """
    Fonction qui calcule et retourne les coordonnées du vecteur Rij,
    à partir des coordinnées des atomes i et j
    
    param:
        list1: une liste correspondant aux coordonées x, y et z de l'atome i
        list2: une liste correspondant aux coordonées x, y et z de l'atome j
    
    return:
        une liste correspondant au vecteur(i,j)
    """
    
    #on calcul les valeurs du vecteur
    x = (list1[0]-list2[0])
    y = (list1[1]-list2[1])
    z = (list1[2]-list2[2])
    return([x,y,z])
    
def Norme(list1):
    """
    fonction qui retourne la norme d'un vecteur à partir de ses coordonnées.
   
    param:
        list1: une liste correspondant aux coordonées du vecteurs
    
    return:
        la norme du vecteur
    """
    return(math.sqrt(list1[0]*list1[0]+list1[1]*list1[1]+list1[2]*list1[2]))
    
    
#################################################################################################################
#methode correspondant au format du tp    
#################################################################################################################

def Env_i(i, rc, listXYZ):
    """
    fonction qui calcule les coordonnées atomiques des atomes à moins de rc A de l'atome i
    
    param:
        i: une liste correspondant aux coordonnées de l'atome i
        rc: un entier correspondant à la distance seuil
        listXYZ: une liste composé d'une liste correspondants aux x des atomes , d'une lsite correspondant aux y des atomes, et d'une liste correspondant au z des atomes.
    
    return:
        une liste d'atome sous le même format que listXYZ correspondant aux atome à moin de rc A de i
    """
    
    listRes=[[],[],[]] 
    
    #on itere sur toutes les coordonées
    for j in range(len(listXYZ[0])):
        
        if (i != [listXYZ[0][j], listXYZ[1][j], listXYZ[2][j]]): #traite le cas ou on a comparer le même atome
        
            #on récupère la valeurs de la distance
            dist = Norme(Rij(i, [listXYZ[0][j], listXYZ[1][j], listXYZ[2][j]] ))
            
            #si la distance ainsi calculé est inférieur au seuil, on l'ajoute au vecteur résultat
            if((dist!=0) and (dist<=rc)):
                
                listRes[0].append(listXYZ[0][j])
                listRes[1].append(listXYZ[1][j])
                listRes[2].append(listXYZ[2][j])
    
    return (listRes)



def CVi(i, listXYZ):
    """
    fonction calculant la variance circulaire d'un atome i
    
    param:
        i: une liste correspondant aux coordonées de l'atome i
        listXYZ: une liste correspondant aux atomes de la protéines composé d'une liste correspondants aux x des atomes , d'une liste correspondant aux y des atomes, et d'une liste correspondant au z des atomes.
    
    return:
       la variance circulaire de l'atome i 
    """
    #initialisation des valeurs
    res=0
    sommeTmp=[0,0,0]
    
    #on parcours la liste des atomes
    for j in range(len(listXYZ[0])):

        #on calcule un element de la somme
        
        if (i != [listXYZ[0][j], listXYZ[1][j], listXYZ[2][j]]): #traite le cas ou on a comparer le même atome
        
            #on calcule la valeur du vecteur (i,j)
            vectTmp = Rij(i, [listXYZ[0][j], listXYZ[1][j], listXYZ[2][j]] )
            
            #norme du vecteur
            norm=Norme(vectTmp)
            
            #incrementation de la somme
            sommeTmp[0]+=vectTmp[0]/norm
            sommeTmp[1]+=vectTmp[1]/norm
            sommeTmp[2]+=vectTmp[2]/norm
        
    
    #estimation du nombre d'atome proche de i
    nbAtomeProche = len(listXYZ[0])
    
    res=1.0-((1.0/(nbAtomeProche))*Norme(sommeTmp))
    
    return res


def CV_res(listXYZres, listXYZ, seuil = 20):
    """
    fonction qui calcule et retourne la variance circulaire d'un résidu
    
    param:
        listXYZres, une liste correspondants aux atomes du résidu
        listXYZ: une liste contenant les coordonées de tous les atomes d'une proteine
        seuil: le seuil pour déterminer si un atome est proche d'un autre
    
    return:
        la variance circulaire du résidu
    """
    
    listTmp=[]
    
    for i in range(len(listXYZres[0])):
        
        #estimation de l'environnement de l'atome i
        atomeProche = Env_i([listXYZres[0][i],listXYZres[1][i],listXYZres[2][i]], seuil, listXYZ)
        
        #on calcule pour chaque atome du residu sa variance circulaire
        listTmp.append(CVi([listXYZres[0][i],listXYZres[1][i],listXYZres[2][i]], atomeProche))
        
    return mean(listTmp)
    
    
#################################################################################################################
#nouvelle version
#################################################################################################################    
def Env_i2(i, rc, listXYZ):
    """
    fonction qui calcule les coordonnées atomiques des atomes à moins de rc A de l'atome i
    n'ajoute pasl'element i
    
    param:
        i: une liste correspondant aux coordonnées de l'atome i
        rc: un entier correspondant à la distance seuil
        listXYZ: une liste de [x,y,z] correspondant au coordonnées atomique d'un atome
    
    return:
        une liste d'atome sous le même format que listXYZ correspondant aux atome à moin de rc A de i
    """
    
    listRes=[] 
    
    #on itere sur toutes les coordonées
    for j in range(len(listXYZ)):
        
        if (i != listXYZ[j]): #traite le cas ou on a comparer le même atome
        
            #on récupère la valeurs de la distance
            dist = Norme(Rij(i, listXYZ[j] ))
            
            #si la distance ainsi calculé est inférieur au seuil, on l'ajoute au vecteur résultat
            if((dist!=0) and (dist<=rc)):
                
                listRes.append(listXYZ[j])
    
    return (listRes)     
    

def CVi2(i, listXYZ):
    """
    fonction calculant la variance circulaire d'un atome i
    
    param:
        i: une liste correspondant aux coordonées de l'atome i
        listXYZ: une liste de [x,y,z] correspondant aux coordonnées atomique d'un atome
    
    return:
       la variance circulaire de l'atome i 
    """
    #initialisation des valeurs
    res=0
    sommeTmp=[0,0,0]
    
    #on parcours la liste des atomes
    for j in range(len(listXYZ)):

        #on calcule un element de la somme
        
        if (i != listXYZ[j]): #traite le cas ou on a comparer le même atome
        
            #on calcule la valeur du vecteur (i,j)
            vectTmp = Rij(i, listXYZ[j] )
            
            #norme du vecteur
            norm=Norme(vectTmp)
            
            sommeTmp[0]+=vectTmp[0]/norm
            sommeTmp[1]+=vectTmp[1]/norm
            sommeTmp[2]+=vectTmp[2]/norm
        
    
    #estimation du nombre d'atome proche de i
    nbAtomeProche = len(listXYZ)
    
    res=1.0-((1.0/(nbAtomeProche))*Norme(sommeTmp))
    
    return res
    
    
def CV_res2(listXYZres, listXYZ, seuil = 20): #ne traite pas la lsite d'environnement
    """
    fonction qui calcule et retourne la variance circulaire d'un résidu
    
    param:
        listXYZres, une liste correspondants aux coordonnées à tester 
        listXYZ: une liste contenant les coordonées de tous les atomes d'une proteine
        seuil: le seuil pour déterminer si un atome est proche d'un autre
    
    return:
        la variance circulaire du résidu
    """
    
    listTmp=[]
    
    for i in range(len(listXYZres)):
        
        #estimation de l'environnement de l'atome i
        atomeProche = Env_i2(listXYZres[i], seuil, listXYZ)
          
        #on calcule pour chaque atome du residu sa variance circulaire
        listTmp.append(CVi2(listXYZres[i], atomeProche))
        
    return mean(listTmp)    
    

def allatoms_coor_list(dicoPDB):
    coor_list = [[],[],[]]
    
    for chain in dicoPDB["chains"]:
        
        for residue in dicoPDB[chain]["reslist"]:
            
            for atom in dicoPDB[chain][residue]["atomlist"]:
                
                coor_list = numpy.c_[coor_list, [dicoPDB[chain][residue][atom]["x"],
                                                 dicoPDB[chain][residue][atom]["y"],
                                                 dicoPDB[chain][residue][atom]["z"]]]
    return coor_list
    
    

def allbary_coor_list(dicoPDB):
    bary_coor_list = [[],[],[]]
    
    for chain in dicoPDB["chains"]:
        
        for residue in dicoPDB[chain]["reslist"]:
            
            bary_coor_list = numpy.c_[bary_coor_list, bary_res(dicoPDB[chain][residue])]
            
    return bary_coor_list
    
def bary_res(res_dico):

    sum_x = 0
    sum_y = 0
    sum_z = 0
    nb_atom = len(res_dico["atomlist"])

    for atom in res_dico["atomlist"]:
        sum_x += res_dico[atom]["x"]
        sum_y += res_dico[atom]["y"]
        sum_z += res_dico[atom]["z"]

    xbary = sum_x/nb_atom
    ybary = sum_y/nb_atom
    zbary = sum_z/nb_atom

    return([xbary, ybary, zbary])
