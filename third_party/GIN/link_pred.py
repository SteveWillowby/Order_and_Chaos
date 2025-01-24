import math
import os.path as osp
import random
import sys

import torch
from sklearn.metrics import roc_auc_score

import torch_geometric.transforms as T
# from torch_geometric.datasets import Planetoid
from torch_geometric.nn import GCNConv, GINConv, GIN
from torch_geometric.utils import negative_sampling

if torch.cuda.is_available():
    device = torch.device('cuda')
elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
    device = torch.device('mps')
else:
    device = torch.device('cpu')

"""
transform = T.Compose([
    T.NormalizeFeatures(),
    T.ToDevice(device),
    T.RandomLinkSplit(num_val=0.05, num_test=0.1, is_undirected=True,
                      add_negative_train_samples=False),
])
path = osp.join(osp.dirname(osp.realpath(__file__)), '..', 'data', 'Planetoid')
dataset = Planetoid(path, name='Cora', transform=transform)
# After applying the `RandomLinkSplit` transform, the data is transformed from
# a data object to a list of tuples (train_data, val_data, test_data), with
# each element representing the corresponding split.
train_data, val_data, test_data = dataset[0]
"""

assert len(sys.argv) == 3
edgelist_filename = sys.argv[1]
output_filename = sys.argv[2]

class LocalUndirectedGraph:

    def __init__(self, edgelist_filename, add_neg_samples=False, seed=540):
        f = open(edgelist_filename, "r")
        lines = f.readlines()
        f.close()

        lines = [l.strip().split(" ") for l in lines]
        edges = [(int(l[0]), int(l[1])) for l in lines]
        edges = set([(min(a, b), max(a, b)) for (a, b) in edges])
        nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])
        N = max(nodes) + 1
        M = len(edges)
        
        self.num_nodes = N
        self.num_edges = M
        self.edges = edges

        input_mode = ["1_hot", "random"][0]

        if input_mode == "1_hot":
            self.x = torch.tensor([[float(int(i == j)) for j in range(0, N)] for i in range(0, N)]).to(device)
            self.num_features = N

        else:
            random.seed(seed)
            self.num_features = (int(math.ceil(math.log2(N))) + 1) * ((M * 2) // N)
            self.x = torch.tensor([[random.randint(0, 1) for _ in range(0, self.num_features)] for i in range(0, N)]).to(device)

        edges_both_ways = list(edges) + [(b, a) for (a, b) in edges]
        self.edge_index = torch.tensor([[edges_both_ways[j][i] for j in range(0, 2 * M)] for i in range(0, 2)]).to(device)

        if add_neg_samples:
            self.edge_label_index = torch.tensor([[i // N for i in range(0, N * N)], [i % N for i in range(0, N * N)]]).to(device)
            self.edge_label = torch.tensor([float(int((min(i // N, i % N), max(i // N, i % N)) in edges)) for i in range(0, N * N)]).to(device)
        else:
            self.edge_label = torch.tensor([float(1) for _ in range(0, 2 * M)]).to(device)
            self.edge_label_index = torch.tensor([[edges_both_ways[j][i] for j in range(0, 2 * M)] for i in range(0, 2)]).to(device)


train_data = LocalUndirectedGraph(edgelist_filename)
val_data   = LocalUndirectedGraph(edgelist_filename, add_neg_samples=True)
# test_data  = LocalUndirectedGraph(edgelist_filename)
dataset    = LocalUndirectedGraph(edgelist_filename)

"""
print(train_data.x.shape)
print(train_data.edge_index.shape)
row_1 = [int(v.item()) for v in list(train_data.edge_index[0,:])]
row_2 = [int(v.item()) for v in list(train_data.edge_index[1,:])]
# print(min(row_1))
# print(max(row_1))
edges = set([(min(row_1[i], row_2[i]), max(row_1[i], row_2[i])) for \
            i in range(0, len(row_1))])
print(len(edges))
print(len(row_1))
print(train_data.y.shape)
print(train_data.edge_label.shape)
print("Number of edge labels = %d" % (len(set([int(v.item()) for v in list(train_data.edge_label)]))))
"""

class Net(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)

    def encode(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv2(x, edge_index)

    def decode(self, z, edge_label_index):
        return (z[edge_label_index[0]] * z[edge_label_index[1]]).sum(dim=-1)

    def decode_all(self, z):
        prob_adj = z @ z.t()
        return (prob_adj > 0).nonzero(as_tuple=False).t()

    def undirected_edges_by_rank(self, z, N):
        prob_adj = z @ z.t()
        all_edges = []
        for i in range(0, N):
            for j in range(i + 1, N):
                edge = (i, j)
                rank = prob_adj[i,j].item()
                all_edges.append((-rank, edge))
        all_edges.sort()
        return [(-rank, edge) for (rank, edge) in all_edges]

class Net2(GIN):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__(in_channels, hidden_channels, 1, out_channels=out_channels, \
                         dropout=0.0)

    def encode(self, x, edge_index):
        v = self.forward(x, edge_index)
        return v

    def decode(self, z, edge_label_index):
        return (z[edge_label_index[0]] * z[edge_label_index[1]]).sum(dim=-1)

    def decode_all(self, z):
        prob_adj = z @ z.t()
        return (prob_adj > 0).nonzero(as_tuple=False).t()

    def undirected_edges_by_rank(self, z, N):
        prob_adj = z @ z.t()
        all_edges = []
        for i in range(0, N):
            for j in range(i + 1, N):
                edge = (i, j)
                rank = prob_adj[i,j].item()
                all_edges.append((-rank, edge))
        all_edges.sort()
        return [(-rank, edge) for (rank, edge) in all_edges]

model = Net2(dataset.num_features, 128, 64).to(device)
optimizer = torch.optim.Adam(params=model.parameters(), lr=0.01)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=50, gamma=0.5)
criterion = torch.nn.BCEWithLogitsLoss()


def train():
    model.train()
    optimizer.zero_grad()
    z = model.encode(train_data.x, train_data.edge_index)

    partial_objective = True

    if partial_objective:
        # We perform a new round of negative sampling for every training epoch:
        neg_edge_index = negative_sampling(
            edge_index=train_data.edge_index, num_nodes=train_data.num_nodes,
            num_neg_samples=train_data.edge_label_index.size(1), method='sparse')

        edge_label_index = torch.cat(
            [train_data.edge_label_index, neg_edge_index],
            dim=-1,
        )
        edge_label = torch.cat([
            train_data.edge_label,
            train_data.edge_label.new_zeros(neg_edge_index.size(1))
        ], dim=0)

    else:
        edge_label       = val_data.edge_label
        edge_label_index = val_data.edge_label_index

    out = model.decode(z, edge_label_index).view(-1)
    loss = criterion(out, edge_label)
    loss.backward()
    optimizer.step()
    return loss


@torch.no_grad()
def test(data):
    model.eval()
    z = model.encode(data.x, data.edge_index)
    out = model.decode(z, data.edge_label_index).view(-1).sigmoid()
    return roc_auc_score(data.edge_label.cpu().numpy(), out.cpu().numpy())

best_val_auc = final_test_auc = 0
for epoch in range(1, 351):
    loss = train()
    val_auc = test(val_data)
    if val_auc > best_val_auc:
        best_val_auc = val_auc
    scheduler.step()
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Val: {val_auc:.4f}')
    sys.stdout.flush()
print("Best val_auc: %f" % best_val_auc)
sys.stdout.flush()

z = model.encode(train_data.x, train_data.edge_index)
# final_edge_index = model.decode_all(z)
edges_by_rank = model.undirected_edges_by_rank(z, train_data.num_nodes)

"""
edges = []
i = 0
while len(edges) < train_data.num_edges:
    (a, b) = edges_by_rank[i]
    i += 1
    if a != b:
        edges.append((a, b))
edges.sort()
"""

best_F1 = None
best_F1_endpoint = None
best_F1_precision = None
best_F1_recall = None
i = 0
n = 0
correct = 0
same_size_endpoint = None
best_size = None
while i < len(edges_by_rank):
    (rank, (a, b)) = edges_by_rank[i]
    if a != b:
        n += 1
        if (a, b) in train_data.edges:
            correct += 1
    if i == (len(edges_by_rank) - 1) or rank != edges_by_rank[i + 1][0]:
        # We have a new F1 evaluation point
        precision = correct / n
        recall = correct / len(train_data.edges)
        F1 = precision * recall
        if best_F1 is None or F1 > best_F1:
            best_F1 = F1
            best_F1_endpoint = i + 1
            best_F1_precision = precision
            best_F1_recall = recall

        diff = abs(n - len(train_data.edges))
        if best_size is None or diff < best_size:
            same_size_endpoint = i + 1
            best_size = diff
    i += 1

print("Best F1 of %.2f with a precision of %.2f and a recall of %.2f" % (best_F1, best_F1_precision, best_F1_recall))
print("Best F1 Endpoint: %d" % best_F1_endpoint)
print("Original number of edges: %d" % len(train_data.edges))
sys.stdout.flush()

# edges = [edge for (rank, edge) in edges_by_rank]
# edges = edges[:best_F1_endpoint]
# edges = edges[:same_size_endpoint]

new_edges_by_rank = []
for (rank, (a, b)) in edges_by_rank:
    if a == b:
        continue
    new_edges_by_rank.append((rank, (a, b)))
edges_by_rank = new_edges_by_rank

f = open(output_filename, "w")
for i in range(0, len(edges_by_rank)):
    (rank, (a, b)) = edges_by_rank[i]
    f.write("%d %d %f" % (a, b, rank))
    if i < len(edges_by_rank) - 1:
        f.write("\n")
f.close()
