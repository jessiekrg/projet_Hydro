#Importation des modules nécessaires de Biopython
from Bio.PDB import PDBParser, PPBuilder
import matplotlib.pyplot as plt


#Étape 1 : Initialiser le parser pour lire le fichier PDB
parser = PDBParser(QUIET=True)  # QUIET=True pour éviter les avertissements lors du parsing

#Étape 2 : Charger la structure depuis le fichier PDB
structure = parser.get_structure("structure_proteine", "3d9s.pdb")

#Étape 3 : Initialiser le constructeur polypeptides
ppb = PPBuilder()
seq = []

#Étape 4 : Parcourir les polypeptides construits et extraire les séquences
for chaine in ppb.build_peptides(structure):
    sequence = chaine.get_sequence()  #Obtient la séquence d'acides aminés d'une chaine
    seq.append(sequence) # Ajoute une chaine à la liste séquence (4 : A,B,C,D)

#Étape 5 : Echelle d'hydrophobicité
Indice_H = {
    'I':4.5, 
}

#Étape x : Calculer l'indice d'hydrophobicité de chaque aa
def Calcul_Profil_H(Chaine):
    return [Indice_H[aa] for aa in Chaine] # Retourne une liste de la valeur de chaque aa

#Étape x : Maplotlib