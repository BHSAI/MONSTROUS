import csv
import json
from constants import *

def write_csv(compounds, outputChannel, proteinGroup, applyApplicabilityDomainFilter):
    writer = csv.writer(outputChannel, lineterminator='\n')

    header = ["Name", "SMILES (original)", "SMILES (standardized)"]
    if proteinGroup == "inhibitor": header += inhibitorsGCNN + inhibitorsSEA
    if proteinGroup == "substrate": header += substratesGCNN + substratesSEA
    _ = writer.writerow(header)

    for compound in compounds:
        row = [compound.name, compound.original, compound.smiles]
        if proteinGroup == "inhibitor" and applyApplicabilityDomainFilter:
            for protein in inhibitorsGCNN: 
                if compound.outsideApplicabilityDomainGCNN == None or compound.outsideApplicabilityDomainGCNN[protein]: row.append("-")
                else: 
                    if protein in compound.resultsGCNN: row.append(compound.resultsGCNN[protein])
                    else: row.append("-")
            for protein in inhibitorsSEA:
                if compound.outsideApplicabilityDomainSEA == None or compound.outsideApplicabilityDomainSEA[protein]: row.append("-")
                else: 
                    if protein in compound.resultsSEA: row.append(compound.resultsSEA[protein])
                    else: row.append("-")
        elif proteinGroup == "substrate" and applyApplicabilityDomainFilter:
            for protein in substratesGCNN: 
                if compound.outsideApplicabilityDomainGCNN == None or compound.outsideApplicabilityDomainGCNN[protein]: row.append("-")
                else: 
                    if protein in compound.resultsGCNN: row.append(compound.resultsGCNN[protein])
                    else: row.append("-")
            for protein in substratesSEA:
                if compound.outsideApplicabilityDomainSEA == None or compound.outsideApplicabilityDomainSEA[protein]: row.append("-")
                else: 
                    if protein in compound.resultsSEA: row.append(compound.resultsSEA[protein])
                    else: row.append("-")
        elif proteinGroup == "inhibitor":
            for protein in inhibitorsGCNN: 
                if protein in compound.resultsGCNN: row.append(compound.resultsGCNN[protein])
                else: row.append("-")
            for protein in inhibitorsSEA:
                if protein in compound.resultsSEA: row.append(compound.resultsSEA[protein])
                else: row.append("-")
        elif proteinGroup == "substrate":
            for protein in substratesGCNN: 
                if protein in compound.resultsGCNN: row.append(compound.resultsGCNN[protein])
                else: row.append("-")
            for protein in substratesSEA:
                if protein in compound.resultsSEA: row.append(compound.resultsSEA[protein])
                else: row.append("-")
        _ = writer.writerow(row)


def write_json(compounds, outputChannel):
    jsonList = []
    for compound in compounds:
        jsonList.append({"name" : compound.name,
                         "original_smiles" : compound.original,
                         "standardized_smiles" : compound.smiles,
                         "results_gcnn" : compound.resultsGCNN,
                         "results_sea" : compound.resultsSEA,
                         "outside_applicability_domain_gcnn" : compound.outsideApplicabilityDomainGCNN,
                         "outside_applicability_domain_sea" : compound.outsideApplicabilityDomainSEA,
                         "filtered" : compound.filtered
                         })
    res_dict = {"compounds": jsonList}
    print(json.dumps(res_dict), file=outputChannel)
    

def write_database(compounds, outputChannel):
    writer = csv.writer(outputChannel, lineterminator='\n')

    header = ["name", "smiles_original", "smiles_std", "inhibitor/substrate", "protein", "result", "outside_applicability_domain"]
    _ = writer.writerow(header)

    for compound in compounds:
        val = "-"
        appDomain = "-"
        for protein in inhibitorsGCNN:
            if protein in compound.resultsGCNN: val = compound.resultsGCNN[protein]
            else: val = "-"
            if compound.outsideApplicabilityDomainGCNN != None and protein in compound.outsideApplicabilityDomainGCNN: appDomain = compound.outsideApplicabilityDomainGCNN[protein]
            else: appDomain = "-"
            row = [compound.name, compound.original, compound.smiles, "inhibitor", protein, val, appDomain]
            _ = writer.writerow(row)
        for protein in inhibitorsSEA:
            if protein in compound.resultsSEA: val = compound.resultsSEA[protein]
            else: val = "-"
            if compound.outsideApplicabilityDomainSEA != None and protein in compound.outsideApplicabilityDomainSEA: appDomain = compound.outsideApplicabilityDomainSEA[protein]
            else: appDomain = "-"
            row = [compound.name, compound.original, compound.smiles, "inhibitor", protein, val, appDomain]
            _ = writer.writerow(row)
        for protein in substratesGCNN:
            if protein in compound.resultsGCNN: val = compound.resultsGCNN[protein]
            else: val = "-"
            if compound.outsideApplicabilityDomainGCNN != None and protein in compound.outsideApplicabilityDomainGCNN: appDomain = compound.outsideApplicabilityDomainGCNN[protein]
            else: appDomain = "-"
            row = [compound.name, compound.original, compound.smiles, "substrate", protein, val, appDomain]
            _ = writer.writerow(row)
        for protein in substratesSEA:
            if protein in compound.resultsSEA: val = compound.resultsSEA[protein]
            else: val = "-"
            if compound.outsideApplicabilityDomainSEA != None and protein in compound.outsideApplicabilityDomainSEA: appDomain = compound.outsideApplicabilityDomainSEA[protein]
            else: appDomain = "-"
            row = [compound.name, compound.original, compound.smiles, "substrate", protein, val, appDomain]
            _ = writer.writerow(row)

    

def write_output(fileFormat, compounds, outputChannels):
    if fileFormat == 'csv':
        write_csv(compounds, outputChannels[0], "inhibitor", applyApplicabilityDomainFilter=True)
        write_csv(compounds, outputChannels[1], "substrate", applyApplicabilityDomainFilter=True)
        write_csv(compounds, outputChannels[2], "inhibitor", applyApplicabilityDomainFilter=False)
        write_csv(compounds, outputChannels[3], "substrate", applyApplicabilityDomainFilter=False)
    elif fileFormat == 'json':
        write_json(compounds, outputChannels[0])
    else: # format is database
        write_database(compounds, outputChannels[0])