from flask import Blueprint, request, jsonify
from utils.inference import predict_bond
from utils.chem_utils import get_basic_mol_info

predict_bp = Blueprint('predict_bp', __name__)

@predict_bp.route('/predict_bond_reactivity', methods=['POST'])
def bond_prediction():
    data = request.get_json()
    smiles = data.get("smiles", "")

    if not smiles:
        return jsonify({"error": "SMILES input required."}), 400

    try:
        pred = predict_bond(smiles)
        meta = get_basic_mol_info(smiles)
        return jsonify({
            "prediction": pred,
            "molecule_info": meta
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500
