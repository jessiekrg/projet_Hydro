# Profil d'hydrophobicité d'une protéine à partir d'un fichier PDB

Projet réalisé par Lina et Jessica.

Ce projet trace le **profil d'hydrophobicité** d'une protéine à partir de sa structure 3D (fichier PDB), en utilisant l'échelle de **Kyte & Doolittle** et une **moyenne glissante**. L'exemple traité est l'**aquaporine 5 (AQP5)**, une protéine membranaire (PDB : `3d9s`).

## Contexte

L'hydrophobicité des acides aminés influence fortement le repliement des protéines. Pour une protéine membranaire, les segments hydrophobes correspondent souvent aux domaines transmembranaires (hélices α). Le profil d'hydrophobicité permet donc de repérer ces segments le long de la séquence.

L'échelle de Kyte & Doolittle (1982) attribue à chaque acide aminé une valeur :

- valeur **positive** → acide aminé hydrophobe
- valeur **négative** → acide aminé hydrophile

## Contenu du dossier

| Fichier | Description |
|---|---|
| `LinaJessica.py` | Script Python principal |
| `3d9s.pdb` | Structure de l'AQP5 (tétramère, chaînes A, B, C, D) téléchargée sur [RCSB PDB](https://www.rcsb.org/) |
| `Rapportfinal.pdf` | Rapport complet (introduction, résultats, discussion, bibliographie) |
| `Figure_1.png` | Profil obtenu avec le script (fenêtre de 20) |
| `Echelle_Hydrophobicite.png` | Tableau comparatif de trois échelles d'hydrophobicité (Kyte-Doolittle, Wimley-White, Hessa) |

## Fonctionnement du script

1. **Lecture du PDB** : `PDBParser` de Biopython charge `3d9s.pdb`.
2. **Extraction des séquences** : `PPBuilder` construit les polypeptides et en récupère la séquence d'acides aminés (4 chaînes : A, B, C, D).
3. **Échelle d'hydrophobicité** : un dictionnaire `Indice_H` associe à chaque acide aminé (code 1 lettre) sa valeur de Kyte-Doolittle.
4. **Conversion** : la fonction `Calcul_Profil_H` transforme chaque séquence en liste de valeurs ; les 4 chaînes sont **concaténées** en une seule liste.
5. **Lissage** : moyenne glissante centrée sur chaque position. Avec `fenetre = 20`, on prend jusqu'à 10 résidus avant et 10 après ; aux extrémités, la moyenne est calculée sur les valeurs disponibles.
6. **Tracé** avec Matplotlib : hydrophobicité moyenne en fonction de la position dans la séquence, avec une ligne de référence à 0.

## Installation et utilisation

Prérequis : Python 3 et les bibliothèques suivantes.

```bash
pip install biopython matplotlib
```

Lancer le script depuis le dossier du projet (le fichier `3d9s.pdb` doit être dans le même dossier) :

```bash
python LinaJessica.py
```

Une fenêtre s'ouvre avec le graphique du profil d'hydrophobicité lissé.

### Paramètres modifiables

- **Taille de la fenêtre** : variable `fenetre` (étape 7). Le rapport présente une fenêtre de 9 (Figure 1 du rapport) et une fenêtre de 20 (Figure 2 du rapport). La version actuelle du script utilise 20.
- **Protéine** : remplacer le fichier PDB dans `parser.get_structure(...)`.
- **Échelle** : modifier le dictionnaire `Indice_H` (voir `Echelle_Hydrophobicite.png` pour d'autres échelles).

## Résultats

Le profil de l'AQP5 (tétramère, environ 975 résidus) montre une alternance de **pics hydrophobes** et de **creux hydrophiles**. Les pics sont compatibles avec les 6 hélices transmembranaires (TM1 à TM6) de chaque monomère, ce qui soutient le caractère de protéine membranaire intégrale. Les creux correspondent plutôt à des régions exposées à l'eau ou aux boucles impliquées dans le passage sélectif de l'eau. Le détail de l'analyse est dans `Rapportfinal.pdf`.

## Limites et pistes d'amélioration

- Les chaînes A à D sont mises bout à bout : la moyenne glissante « déborde » donc d'une chaîne sur l'autre aux jonctions.
- Les zones transmembranaires ne sont pas détectées automatiquement ni surlignées sur le graphique (étapes 6 et 7 de la stratégie du rapport) ; seule la lecture visuelle du profil est utilisée.
- `PPBuilder` ignore les résidus non standard ou les ruptures de chaîne, et le dictionnaire ne gère que les 20 acides aminés standards.
- Le chemin du fichier PDB est écrit en dur dans le script.

## Références

- Kyte, J. & Doolittle, R. F. (1982). *A simple method for displaying the hydropathic character of a protein*. Journal of Molecular Biology, 157(1), 105-132. https://doi.org/10.1016/0022-2836(82)90515-0
- RCSB Protein Data Bank : https://www.rcsb.org/
- ExPASy ProtScale : https://web.expasy.org/protscale/

La bibliographie complète figure dans le rapport.