# inference.py
import torch
from src.data_utils import mol_to_graph
from utils.model_loader import load_model

def predict_bond(smiles: str):
    model = load_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        logits = model(data.x, data.edge_index, data.edge_attr)
        probs = torch.sigmoid(logits)

        seen = set()
        bond_confidences = []

        for i in range(data.edge_index.shape[1]):
            atom1 = data.edge_index[0][i].item()
            atom2 = data.edge_index[1][i].item()

            # Ensure undirected uniqueness
            bond = tuple(sorted((atom1, atom2)))
            if bond not in seen:
                seen.add(bond)
                bond_confidences.append((atom1, atom2, probs[i].item()))

        if not bond_confidences:
            return {"predicted_bond": None, "confidence": 0.0}

        # Select bond with highest reactivity
        atom1, atom2, confidence = max(bond_confidences, key=lambda x: x[2])

        return {
            "predicted_bond": [atom1, atom2],
            "confidence": round(confidence, 3)
        }
