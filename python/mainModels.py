import csv
import os
import chemprop
from compound import Compound
from preprocessing import *
from constants import *
from rdkit import DataStructs
import math

dir_path = os.path.dirname(os.path.realpath(__file__))

def getGCNNCompounds(channels):# Returns a map where keys are the GCNN proteins and the values are lists of Compounds (Compounds must have a fingerprint)
    proteins = proteinsGCNN
    allCompounds = {}
    for protein in proteins:
        # csv_path = f"/server/plugin/monstrous/gcnn/compounds/{protein}.csv"
        csv_path = f"{dir_path}/../gcnn/compounds/{protein}.csv"
        with open(csv_path, 'r') as file:
            csvreader = csv.reader(file)
            compounds = []
            for row in csvreader:
                compounds.append(Compound("", row[1]))
            compounds.pop(0) #remove header row from compounds
        standardize(compounds, channels)
        for compound in compounds:
            fingerprint(compound)
        allCompounds[protein] = compounds
    return allCompounds

def getSEACompounds(channels):# Returns a map where keys are the SEA proteins and the values are lists of Compounds (Compounds must have a fingerprint) NOTE: these are the app domain compounds
    print("Getting Similarity Approach compounds...")
    proteins = proteinsSEA
    allCompounds = {}
    for protein in proteins:
        print("\tGetting compounds for " + protein)
        # csv_path = f"/server/plugin/monstrous/similarityApproach/compounds/{protein}.csv"
        csv_path = f"{dir_path}/../similarityApproach/compounds/{protein}.csv"
        with open(csv_path, 'r') as file:
            csvreader = csv.reader(file)
            compounds = []

            smilesCol = 1
            if (protein.split("_")[1] == "inhibitor"): smilesCol = 2
            print(protein + "'s smiles col: " + str(smilesCol))

            for row in csvreader:
                compounds.append(Compound("", row[smilesCol]))
            compounds.pop(0) #remove header row from compounds
        print("\tFound " + str(len(compounds)) + " compounds for " + protein)
        standardize(compounds, channels)
        for compound in compounds:
            fingerprint(compound)
        allCompounds[protein] = compounds
    return allCompounds

def applicabilityDomain(compounds, channels):
    gcnnCompounds = getGCNNCompounds(channels)
    seaCompounds = getSEACompounds(channels)

    for compound in compounds:
        if compound.fp != None:
            outsideApplicabilityDomainGCNN = {}
            for gcnnProtein in proteinsGCNN:
                try:
                    maxTanimotoSimilarity = 0
                    for gcnnCompound in gcnnCompounds[gcnnProtein]:
                        maxTanimotoSimilarity = max(maxTanimotoSimilarity, DataStructs.FingerprintSimilarity(compound.fp, gcnnCompound.fp, metric=DataStructs.TanimotoSimilarity))
                    outsideApplicabilityDomainGCNN[gcnnProtein] = maxTanimotoSimilarity > APPLICABILITY_DOMAIN_THRESHOLD
                except:
                    channels.writeErrorMsg(f"Failed to generate applicability domain for: {gcnnProtein}, {compound.original}")
            compound.outsideApplicabilityDomainGCNN = outsideApplicabilityDomainGCNN
            outsideApplicabilityDomainSEA = {}
            for seaProtein in proteinsSEA:
                try:
                    maxTanimotoSimilarity = 0
                    for seaCompound in seaCompounds[seaProtein]:
                        maxTanimotoSimilarity = max(maxTanimotoSimilarity, DataStructs.FingerprintSimilarity(compound.fp, seaCompound.fp, metric=DataStructs.TanimotoSimilarity))
                    outsideApplicabilityDomainSEA[seaProtein] = maxTanimotoSimilarity > APPLICABILITY_DOMAIN_THRESHOLD
                except:
                    channels.writeErrorMsg(f"Failed to generate applicability domain for: {seaProtein}, {compound.original}")
            compound.outsideApplicabilityDomainSEA = outsideApplicabilityDomainSEA

def run_sea(compounds, channels):
    seaCompounds = getSEACompounds(channels)
    for compound in compounds:
        try:
            resultsSEA = {}
            if compound.fp != None:
                for protein in proteinsSEA:
                    rs = 0
                    for seaCompound in seaCompounds[protein]:
                        ts = DataStructs.FingerprintSimilarity(compound.fp, seaCompound.fp, metric=DataStructs.TanimotoSimilarity)
                        if ts >= 0.57:
                            rs += ts
                    resultValue = calcP(rs, len(seaCompounds[protein]))
                    resultsSEA[protein] = (resultValue <= SEA_THRESHOLD)
                compound.resultsSEA = resultsSEA
        except:
            channels.writeErrorMsg(f"Failed to run Similarity Approach predictions for SMILES: {compound.original}")

def calcP(sumOfRelevantSimilarities, numMols):
    z = (sumOfRelevantSimilarities - (0.000424 * numMols)) / (0.00449 * (numMols ** 0.665))
    xz = -(math.e ** ((-1 * z * math.pi) / (math.sqrt(6) - 0.577215665)))
    if z <= 28:
        return (1 - (math.e ** xz))
    else:
        return ((-1 * xz) - ((xz ** 2) / 2) - ((xz ** 3) / 6))

def divideCompounds(compounds, batchSize):
    for i in range(0, len(compounds), batchSize):
        yield compounds[i:i + batchSize]

def GCNN(compounds, batchSize, channels):
    compound_batch_list = list(divideCompounds(compounds, batchSize))

    proteins = proteinsGCNN
    for i in range(len(proteins)):
        channels.writeStateMsg(f"Running GCNN protein {i}/{len(proteins)}...")
        batchNum = 0
        print(f"Running GCNN for {proteins[i]} ({i + 1}/{len(proteins)})")
        for compound_batch in compound_batch_list:
            batchNum += 1
            filteredBatch = [compound for compound in compound_batch if not compound.filtered]

            smiles = list(map(lambda n:[n.smiles], filteredBatch))
            try:  
                arguments = [
                '--test_path', '/dev/null',
                '--preds_path', '/dev/null',
                '--checkpoint_dir', f"{dir_path}/../gcnn/models/{proteins[i]}"
                ]
                args = chemprop.args.PredictArgs().parse_args(arguments)
                #preds = [list(x) for x in zip(*chemprop.train.make_predictions(args=args, smiles=smiles))][0]
                resultsLists = chemprop.train.make_predictions(args=args, smiles=smiles) # runs gcnn, returns each result in its own list, all within one final list
                resultsListTransposed = zip(*resultsLists) # transposes the resultsLists into an iterator of length 1, containing the results as a tuple
                preds = [list(x) for x in resultsListTransposed][0] # converts resultsListTransposed from an iterator of tuples to a list of lists
                for j in range(len(filteredBatch)):
                    filteredBatch[j].resultsGCNN[proteins[i]] = preds[j] >= 0.5
            except:
                channels.writeErrorMsg(f"Failed to run GCNN predictions for protein: {proteins[i]}, batch: {batchNum} (Compounds {(batchNum-1)*batchSize + 1}-{min(len(compounds), (batchNum)*batchSize)}) ")

def runMainModels(compounds, batch_size, channels):
    channels.writeStateMsg("Running GCNN models...")
    GCNN(compounds, batch_size, channels)
    channels.writeStateMsg("Running Similarity Approach models...")
    run_sea(compounds, channels)
    channels.writeStateMsg("Finding applicability domain...")
    applicabilityDomain(compounds, channels)
