import torch
from src.data_utils import mol_to_graph
from utils.model_loader import load_model

def predict_bond(smiles: str):
    model = load_model()
    data = mol_to_graph(smiles)

    with torch.no_grad():
        logits = model(data.x, data.edge_index, data.edge_attr)
        probs = torch.sigmoid(logits)
        top_idx = torch.argmax(probs).item()

        atom1 = data.edge_index[0][top_idx].item()
        atom2 = data.edge_index[1][top_idx].item()
        confidence = probs[top_idx].item()

        return {
            "predicted_bond": [atom1, atom2],
            "confidence": round(confidence, 3)
        }
