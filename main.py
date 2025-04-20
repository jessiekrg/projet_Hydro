# Importation des modules nécessaires de Biopython
from Bio.PDB import PDBParser, PPBuilder

# Étape 1 : Initialiser le parser pour lire le fichier PDB
parser = PDBParser(QUIET=True)  # QUIET=True pour éviter les avertissements lors du parsing

# Étape 2 : Charger la structure depuis le fichier PDB
# Remplacez "7cjs.pdb" par le nom de votre fichier PDB
structure = parser.get_structure("structure_proteine", "7cjs.pdb")

# Étape 3 : Initialiser le constructeur de polypeptides
ppb = PPBuilder()

# Étape 4 : Parcourir les polypeptides construits et extraire les séquences
for polypeptide in ppb.build_peptides(structure):
    sequence = polypeptide.get_sequence()  # Obtient la séquence d'acides aminés
    print(sequence)  # Affiche la séquence

