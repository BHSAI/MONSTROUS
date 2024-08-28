class Compound:
    def __init__(self, name, original):
        self.name = name
        self.original = original
        self.smiles = None
        self.fp = None
        self.resultsGCNN = {}
        self.resultsSA = {}
        self.outsideApplicabilityDomainGCNN = None
        self.outsideApplicabilityDomainSA = None
        self.filtered = False # Marks compounds that should be filtered out of the GCNN due to having too small a molecular weight or having inorganic atoms
    
def getCompounds(names, smiles):
    compounds = []
    for i in range(len(smiles)):
        compounds.append(Compound(names[i], smiles[i]))
    return compounds