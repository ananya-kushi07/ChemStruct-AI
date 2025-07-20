from flask import Blueprint, request, jsonify
from utils.inference import predict_bond
from utils.chem_utils import get_basic_mol_info, resolve_to_smiles

predict_bp = Blueprint('predict_bp', __name__)

@predict_bp.route('/predict_bond_reactivity', methods=['POST'])
def bond_prediction():
    data = request.get_json()
    input_str = data.get("smiles", "").strip()

    if not input_str:
        return jsonify({"error": "Molecule input is required (SMILES or name)."}), 400

    try:
        # Try resolving input to canonical SMILES (name → SMILES, formula → SMILES, or pass through if valid)
        smiles = resolve_to_smiles(input_str)
        if not smiles:
            return jsonify({"error": "Could not resolve input to a valid molecule (SMILES)."}), 400

        pred = predict_bond(smiles)
        meta = get_basic_mol_info(smiles)
        return jsonify({
            "input_used": input_str,
            "smiles_used": smiles,
            "prediction": pred,
            "molecule_info": meta
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500
