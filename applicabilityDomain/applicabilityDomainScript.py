import argparse
import csv
import os
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardizer
from rdkit.Chem import rdFingerprintGenerator

def argParse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", nargs=1, help="Location of the csv containing the SMILES you wish to test")
    parser.add_argument("-p", "--proteinFile", nargs=1, help="Location of the csv containing the compounds related to your protein")
    parser.add_argument("-o", "--outputFile", nargs=1, help="File location for the output csv")

    args = parser.parse_args()

    if args.inputFile == None or args.inputFile[0] == None:
        parser.error("Input file needed to run applicability domain tool")
    if args.outputFile == None or args.outputFile[0] == None:
        parser.error("Output file needed to run applicability domain tool")
    if args.proteinFile == None or args.proteinFile[0] == None:
        parser.error("Protein file needed to run applicability domain tool")
    return args.inputFile[0], args.proteinFile[0], args.outputFile[0]

class Compound:
    def __init__(self, name, original):
        self.name = name
        self.original = original
        self.smiles = None
        self.fp = None


def standardize(compounds):
    for compound in compounds:
        try:
            compound.smiles = standardize_one_smiles(compound.original) 
        except:
            print(f"Failed to standardize SMILES: {compound.original}")

def standardize_one_smiles(smiles):
    largest_Fragment = rdMolStandardize.LargestFragmentChooser()
    if smiles is not None:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            largest_mol = largest_Fragment.choose(molecule)
            m_no_salt = standardizer.get_parent_mol(largest_mol)
            molecule_std = standardizer.standardize_mol(m_no_salt[0])
            smiles_std = Chem.MolToSmiles(molecule_std)
            return smiles_std
    return None

def fingerprint(compound):
    if compound.smiles != None:

        mol = Chem.MolFromSmiles(compound.smiles)

        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        compound.fp = mfpgen.GetFingerprint(mol)

def getCompounds(path):
    with open(path, 'r') as file:
        csvreader = csv.reader(file)
        compounds = []
        for row in csvreader:
            compounds.append(Compound(row[0], row[1]))
        compounds.pop(0)
    standardize(compounds)
    for compound in compounds:
        fingerprint(compound)
    return compounds

def applicabilityDomain(testCompounds, proteinCompounds):
    for compound in testCompounds:
        if compound.fp != None:
            sdc = {}
            maxTanimotoSimilarity = 0
            for proteinCompound in proteinCompounds:
                maxTanimotoSimilarity = max(maxTanimotoSimilarity, DataStructs.FingerprintSimilarity(compound.fp, proteinCompound.fp, metric=DataStructs.TanimotoSimilarity))
            sdc = maxTanimotoSimilarity
            compound.sdc = sdc

def write_csv(compounds, outputChannel):
    writer = csv.writer(outputChannel, lineterminator='\n')

    header = ["Name", "SMILES (original)", "SMILES (standardized)", "SDC"]
    _ = writer.writerow(header)

    for compound in compounds:
        row = [compound.name, compound.original, compound.smiles, compound.sdc]
        _ = writer.writerow(row)

if __name__ == "__main__":

    testCompoundsPath, proteinCompoundsPath, outputPath = argParse()
    output = open(outputPath, 'w')

    testCompounds = getCompounds(testCompoundsPath)
    proteinCompounds = getCompounds(proteinCompoundsPath)

    applicabilityDomain(testCompounds, proteinCompounds)
    write_csv(testCompounds, output)

    output.close()
    