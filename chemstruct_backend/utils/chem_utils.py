

from rdkit import Chem
from rdkit.Chem import Descriptors

def get_basic_mol_info(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return {
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds()
    }
