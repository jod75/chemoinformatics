#!/usr/bin/env python

###############################################################################
# Tanimoto Similarity test using RDKit and Zinc15 database
# Joseph D'Emanuele
#

import urllib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw


###############################################################################
# This runs through the molecules in a SMILES file and returns
# a list of Molecules.
def process_smiles_file(filename):
    smiles = {}
    with open(filename, "r") as infile:
        infile.readline()  # skip header
        for line in infile:
            parts = line.split()
            m = Chem.MolFromSmiles(parts[0])
            if m is None:
                continue
            smiles[parts[1]] = m
    return smiles

###############################################################################
# similarity test

# download smiles file
urllib.urlretrieve("http://files.docking.org/2D/AA/AAAA.smi", "../data/AAAA.smi")

# process file and create a list of Molecules
molecules = process_smiles_file("../data/AAAA.smi")

# use first molecule as query fingerprint
fp_query = AllChem.GetMorganFingerprintAsBitVect(molecules[molecules.keys()[0]], 2)

# dictionary to keep similarity index
similarities = {}

# compute Tanimoto similarity for all molecules in our file
for moleculeKey in molecules.keys():
    fp2 = AllChem.GetMorganFingerprintAsBitVect(molecules[moleculeKey], 2)
    similarity = DataStructs.FingerprintSimilarity(fp_query, fp2)
    similarities[moleculeKey] = similarity

# get top 20 similar molecules
top20 = sorted(similarities, key=similarities.get, reverse=True)[:20]
top20.insert(0, molecules.keys()[0])  # this is the query molecule

# get bottom 20 similar molecules
bottom20 = sorted(similarities, key=similarities.get, reverse=False)[:20]
bottom20.insert(0, molecules.keys()[0])  # this is the query molecule

# draw top20 similar molecules
img = Draw.MolsToGridImage([molecules[x] for x in top20], molsPerRow=2, subImgSize=(400, 400),
                           legends=["%s - %f" % (x, similarities[x]) for x in top20])
img.save("../out/similarities_top20.png")

# draw bottom20 similar molecules
img = Draw.MolsToGridImage([molecules[x] for x in bottom20], molsPerRow=2, subImgSize=(400, 400),
                           legends=["%s - %f" % (x, similarities[x]) for x in bottom20])
img.save("../out/similarities_bottom20.png")