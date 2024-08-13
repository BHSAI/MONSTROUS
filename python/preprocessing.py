from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardizer
from rdkit.Chem import rdFingerprintGenerator

def standardize(compounds, channels):
    for compound in compounds:
        try:
            compound.smiles = standardize_one_smiles(compound.original) 
        except:
            channels.writeErrorMsg(f"Failed to standardize SMILES: {compound.original}")

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

def organicMWFilter(compound):
    if compound.smiles == None:
        compound.filtered = True
    else:
        mol = Chem.MolFromSmiles(compound.smiles)
        mw = Descriptors.MolWt(mol)
        inorganic = False
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in ["H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]:
                inorganic = True
                break
        compound.filtered = (inorganic or (mw < 50))


def processCompounds(compounds, channels):
    channels.writeStateMsg("Standardizing compounds...")
    standardize(compounds, channels)
    channels.writeStateMsg("Fingerprinting compounds...")
    for compound in compounds:
        try:
            fingerprint(compound)
        except:
            channels.writeErrorMsg(f"Failed to generate fingerprint for SMILES: {compound.original}")
        organicMWFilter(compound)