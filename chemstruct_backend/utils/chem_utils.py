from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import requests


def resolve_to_smiles(input_str):
    """
    Resolve a common name or input (e.g., 'ethanol') to a valid SMILES string using PubChem.
    If input is already valid SMILES, return it directly.
    """
    try:
        # Case 1: Already a valid SMILES
        mol = Chem.MolFromSmiles(input_str)
        if mol is not None:
            return input_str

        # Case 2: Resolve using PubChem name-to-SMILES API
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{input_str}/property/CanonicalSMILES,ConnectivitySMILES/JSON"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            props = response.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            return props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
        else:
            print(f"[resolve_to_smiles] PubChem error code: {response.status_code}")
            return None

    except Exception as e:
        print(f"[resolve_to_smiles ERROR] Failed to resolve: {input_str}. Error: {e}")
        return None


def get_basic_mol_info(input_str):
    """
    Resolve input to SMILES and return basic molecular metadata (formula, weight, atoms, bonds).
    """
    smiles = resolve_to_smiles(input_str)
    if not smiles:
        raise ValueError("Unable to resolve input to valid SMILES.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES is invalid or could not be parsed by RDKit.")

    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 3),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds()
    }
