from torch_geometric.nn import GCNConv
import torch.nn.functional as F
import torch.nn as nn
import torch


class BondBreakGNN(nn.Module):
    def __init__(self, node_feat_dim, hidden_dim=64):
        super().__init__()
        self.conv1 = GCNConv(node_feat_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.lin = nn.Sequential(
            nn.Linear(2 * hidden_dim + 1, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )

    def forward(self, x, edge_index, edge_attr):
        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))

        row, col = edge_index
        edge_inputs = torch.cat([x[row], x[col], edge_attr], dim=1)

        return self.lin(edge_inputs).squeeze(-1)
