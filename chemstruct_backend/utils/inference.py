import torch
from rdkit import Chem
from src.data_utils import mol_to_graph
from utils.model_loader import load_model

def predict_bond(smiles: str):
    model = load_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        logits = model(data.x, data.edge_index, data.edge_attr)
        probs = torch.sigmoid(logits)

        bond_conf_map = {}
        for idx in range(data.edge_index.shape[1]):
            atom1 = data.edge_index[0][idx].item()
            atom2 = data.edge_index[1][idx].item()
            bond = tuple(sorted((atom1, atom2)))  # ensure undirected

            if bond not in bond_conf_map:
                bond_conf_map[bond] = []
            bond_conf_map[bond].append(probs[idx].item())

        if not bond_conf_map:
            return {
                "predicted_bond": None,
                "predicted_bond_atoms": None,
                "confidence": 0.0
            }

        # Average probability of both directions
        bond_confidences = [(a1, a2, sum(confs) / len(confs)) for (a1, a2), confs in bond_conf_map.items()]
        atom1, atom2, confidence = max(bond_confidences, key=lambda x: x[2])

        # Safely convert predicted atom indices to symbols
        try:
            mol = Chem.MolFromSmiles(smiles)
            if atom1 >= mol.GetNumAtoms() or atom2 >= mol.GetNumAtoms():
                raise IndexError("Predicted atom index out of bounds.")
            sym1 = mol.GetAtomWithIdx(atom1).GetSymbol()
            sym2 = mol.GetAtomWithIdx(atom2).GetSymbol()
            bond_str = f"{sym1}({atom1}) - {sym2}({atom2})"
        except Exception as e:
            bond_str = f"{atom1} - {atom2} (atom symbols unavailable)"
            print(f"[WARNING] RDKit indexing issue: {e}")

        return {
            "predicted_bond": [atom1, atom2],
            "predicted_bond_atoms": bond_str,
            "confidence": round(confidence, 3)
        }
