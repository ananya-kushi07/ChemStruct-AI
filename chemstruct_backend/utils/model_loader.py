import torch
from src.model import BondBreakGNN
from config.config import Config

_model = None


def load_model():
    global _model
    if _model is None:
        _model = BondBreakGNN(node_feat_dim=4)
        _model.load_state_dict(torch.load(Config.MODEL_PATH, map_location='cpu'))
        _model.eval()
    return _model
