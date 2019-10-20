#!/usr/bin/env python3

import math


######################################
#
#       Parsers & Writers
#
########################################

def parsePDBMultiChains(infile) :
    """ propos : parser un fichier PDB
        input : un fichier PDB
        outpud : dico de type dPDB"""


    # lecture du fichier PDB
    #-----------------------
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    # var ini
    #---------
    dPDB = {}
    #dPDB["reslist"] = []
    dPDB["chains"] = []

    
    # parcoure le PDB
    #-----------------
    for line in lines :
        if line[0:4] == "ATOM" :

            # on recupere l'info de la chaine
            chain = line[21]

            # si la chaine n'existe pas, on cree la cle correspondante et on ajoute la chaine a la liste des chaines
            if not chain in dPDB["chains"] :
                dPDB["chains"].append(chain) # ajout de "chain" a la liste des chaines
                dPDB[chain] = {} # creation du sous-dico pour la chaine
                # on prepare la structure de donnees pour cette chaine
                dPDB[chain]["reslist"] = []

            # on recupere l'info du residu    
            curres = "%s"%(line[22:26]).strip()

            # si le residu pour cette chaine "chain" n'existe pas, on cree la cle correspondante et on ajoute le res a la liste des res 
            if not curres in dPDB[chain]["reslist"] :
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                # on prepare la structure de donnees pour ce residu
                dPDB[chain][curres]["atomlist"] = []
                # on recupere l'info du residu
                dPDB[chain][curres]["resname"] = line[17:20].strip()

            # on recupere les info pour l'atome de ce res de cette chaine (type atomique + coords x, y, z)    
            atomtype = line[12:16].strip()
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {} 
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()
            dPDB[chain][curres][atomtype]["bfactor"] = float(line[60:67].strip())

            
    return dPDB

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

#################################################
#           WRITING TOOLS
#################################################


def writePDB(dPDB, filout = "out.pdb", bfactor = False) :
    """purpose: according to the coordinates in dPDB, writes the corresponding PDB file.
       If bfactor = True, writes also the information corresponding to the key bfactor
       of each atom in dPDB.
       input: a dico with the dPDB format
       output: PDB file.
    """

    fout = open(filout, "w")

    for chain in dPDB["chains"]:
        for res in dPDB[chain]["reslist"] :
            for atom in dPDB[chain][res]["atomlist"] :
                if bfactor :
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00%7.3f X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"],dPDB[chain][res][atom]["bfactor"] ))
                else:
                    fout.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["resname"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"] ))
                    
    fout.close()


######################################
#
#           extract info from PDB
#
########################################



# ceci est une possibilite, vous pouviez aussi renvoyer la sortie avec une liste de listes
# de type [[chainname, nb_aa1],[chainname, nb_aa2],etc] --> [["A", 23],["B", 29],["C", 11]]
def getNbAA(dPDB, aa) :
    """ Purpose: counts the number of aa per chain of a protein
        Input: dPDB (a dico PDB) and the type of aa (string) in 3letter code (i.e. ALA, CYS...)
        Output: a dico containing for each chain, the nb of residues corresponding to "aa" 
    """

    # init var
    d_aaPerChain = {}
    
    for chaini in dPDB["chains"] :
        # initiation de la cle chaini. La cle chaini est pointe vers un compteur initialise a zero
        d_aaPerChain[chaini] = 0

        for resi in dPDB[chaini]["reslist"] :
            # compte nb aa pr la chaine i
            if dPDB[chaini][resi]["resname"] == aa :
                d_aaPerChain[chaini]+= 1

    return d_aaPerChain
                
######################################
#
#           3D manipulations
#
########################################


def distancePoints(L1,L2):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates 
       output: distance"""

    #print(x1, x2)

    x = L1[0]-L2[0]
    y = L1[1]-L2[1]
    z = L1[2]-L2[2]
    return math.sqrt(x*x+y*y+z*z)


def getCentroid(d_res):
    """Purpose: Calculates the center of mass of a residue 
       Input: a dico residue
       Output: coords x, y, z of the centroid (tuple format)
    """
        
    x = y = z = 0.0
    
    # loop over all atoms
    for atom in d_res["atomlist"] :
        x +=d_res[atom]["x"]
        y +=d_res[atom]["y"]
        z +=d_res[atom]["z"]
        
    Xcen = float(x)/len(d_res["atomlist"]) 
    Ycen = float(y)/len(d_res["atomlist"])
    Zcen = float(z)/len(d_res["atomlist"])

    return (Xcen, Ycen, Zcen)    




def compDistance(d_res1, d_res2, mode = "minvalue") :
    """Purpose: computes the distance between 2 residues (res1 and res2)
                Distance can be calculated either as the min distance between all the pairwise distances (atom-atom)
                or between the 2 centroids of the residues
 
       Inputs:  d_res1, d_res2 which are dico corresponding to res 1 and res 2 respectively
                mode: either "minvalue" or "centroid" (string)
                
       Output:  distance (float) """
     
    # mode minvalue: min distance atom-atom
    #--------------------------------------
    if mode == "minvalue" : 
        distance = 1000000

        #*** computes all pairwise distances between atoms of res1, res2

        # loop over all atoms from res1
        for atom1 in d_res1["atomlist"] :
            coord1 = [d_res1[atom1]["x"], d_res1[atom1]["y"], d_res1[atom1]["z"]]
            # lloop over all atoms from res2
            for atom2 in d_res2["atomlist"] :
                coord2 = [d_res2[atom2]["x"], d_res2[atom2]["y"], d_res2[atom2]["z"]]
                # computes distance between atom1 from res1 and atom2 from res2
                dist_tmp = distancePoints((coord1[0], coord1[1], coord1[2]),(coord2[0],coord2[1], coord2[2]))
                if distance > dist_tmp :
                    distance = dist_tmp

    # mode centroid: distance between the centroids of the 2 given residues
    #----------------------------------------------------------------------
    elif mode == "centroid" : 
        cent1 = getCentroid(d_res1)
        cent2 = getCentroid(d_res2)
        #print(cent1, cent2)

        distance = distancePoints(cent1, cent2)
        
    return distance
