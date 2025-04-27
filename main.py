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
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
    'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
    'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

#Étape 6 : Calculer l'indice d'hydrophobicité de chaque aa
def Calcul_Profil_H(Chaine):
    return [Indice_H[aa] for aa in Chaine] # Retourne une liste de la valeur de chaque aa

indice = []
for i in seq:
    indice += Calcul_Profil_H(i)

#Etape 7 : Calcul du profil d’hydrophobicité avec une moyenne glissante.

fenetre = 20
liste_profil = []
for i in range (len(indice)):
    f = fenetre//2
    apres = indice[i+1:min(len(indice), i + 1 + f)]
    avant = indice[max(0, i - f): i]
    nbr_valeur = len(avant) + len(apres) + 1
    profil = (sum(apres) + sum(avant) + indice[i])/nbr_valeur

    liste_profil.append(profil)

#Étape x : Maplotlib 
plt.plot(range(len(liste_profil)), liste_profil, color='black')
plt.xlim(0, len(liste_profil))
plt.xticks(range(0, len(liste_profil), 50 ))
plt.xticks(rotation = 45)
plt.ylim(-4.5, 4.5)
plt.title("Profil d'hydrophobicité lissé")
plt.xlabel("Position dans la séquence")
plt.ylabel("Hydrophobicité moyenne")
plt.axhline(0, color='gray', linestyle='--')  # Ligne de séparation
plt.grid(True)
plt.show()


