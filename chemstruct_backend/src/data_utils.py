from rdkit import Chem
import torch
from torch_geometric.data import Data

def mol_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    node_feats = []
    for atom in mol.GetAtoms():
        node_feats.append([
            atom.GetAtomicNum(),
            int(atom.GetIsAromatic()),
            atom.GetTotalNumHs(includeNeighbors=True),
            atom.GetFormalCharge()
        ])
    x = torch.tensor(node_feats, dtype=torch.float)

    edge_index = []
    edge_feats = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        edge_index += [[i, j], [j, i]]
        bond_type = bond.GetBondTypeAsDouble()
        edge_feats += [[bond_type], [bond_type]]

    edge_index = torch.tensor(edge_index).t().contiguous()
    edge_attr = torch.tensor(edge_feats, dtype=torch.float)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
