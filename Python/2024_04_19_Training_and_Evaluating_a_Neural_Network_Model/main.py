# %% Python 3.10.4 | Pytorch 1.12.1 | CUDA 12.4.1
from itertools import chain
from pathlib import Path

import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import roc_auc_score
from torch.utils.data import DataLoader
from torch.utils.data.dataset import TensorDataset

from src.pytorch_nn_demo import NeuralNetwork, device, model_train_eval

# %% Set up paths
ROOT_DIR = Path.cwd()
TRAIN_METADATA_PATH = ROOT_DIR / "input/train_metadata.tsv"
TEST_METADATA_PATH = ROOT_DIR / "input/test_metadata_filtered.tsv"
TRAIN_COUNT_MATRIX_PATH = ROOT_DIR / "input/train_counts_transformed_scaled.tsv.gz"
TEST_COUNT_MATRIX_PATH = ROOT_DIR / "input/test_counts_transformed_scaled.tsv.gz"
FINAL_GENE_LIST_PATH = ROOT_DIR / "input/final_gene_list.tsv"
OUTPUT_DIR = ROOT_DIR / "output"

# %% Import refined gene list
gene_list = pd.read_csv(FINAL_GENE_LIST_PATH, sep="\t")
gene_list = gene_list["index"].to_list()

# %% Import train dataset
train_count_matrix = pd.read_csv(TRAIN_COUNT_MATRIX_PATH, sep="\t", index_col=0)
train_count_matrix = train_count_matrix[train_count_matrix.index.isin(gene_list)]
train_count_matrix = train_count_matrix.sort_index()
train_count_matrix = train_count_matrix.transpose()

train_metadata = pd.read_csv(TRAIN_METADATA_PATH, sep="\t")
train_metadata["label"] = train_metadata["label"].astype(float)

# %% Import test dataset
test_count_matrix = pd.read_csv(TEST_COUNT_MATRIX_PATH, sep="\t", index_col=0)
test_count_matrix = test_count_matrix[test_count_matrix.index.isin(gene_list)]
test_count_matrix = test_count_matrix.sort_index()
test_count_matrix = test_count_matrix.transpose()

test_metadata = pd.read_csv(TEST_METADATA_PATH, sep="\t")
test_metadata["label"] = test_metadata["label"].astype(float)

# %% Convert datasets to PyTorch FloatTensors
train_tensor_dataset = TensorDataset(
    torch.tensor(train_count_matrix.values),
    torch.tensor(train_metadata["label"].values),
)

test_tensor_dataset = TensorDataset(
    torch.tensor(test_count_matrix.values),
    torch.tensor(test_metadata["label"].values),
)
# %% Prepare datasets for batching
train_dataloader = DataLoader(train_tensor_dataset, batch_size=512, shuffle=True)

test_dataloader = DataLoader(test_tensor_dataset, batch_size=512, shuffle=False)

# %% Configure neural network
n_features = train_count_matrix.shape[1]
n_output = 1
n_neurons = int((n_features + n_output) * (2 / 3))  # rule of thumb
dropout_rate = 0.2
activation_function = nn.ReLU()

model = NeuralNetwork(
    n_features=n_features,
    n_neurons=n_neurons,
    n_output=n_output,
    activation_function=activation_function,
    dropout=dropout_rate,
).to(device)

# %% Run training/evaluation iterations
optimizer = optim.Adam(model.parameters(), lr=0.001)
loss_function = nn.BCELoss()

model_train_eval(
    model=model,
    n_epochs=250,
    train_dataloader=train_dataloader,
    test_dataloader=test_dataloader,
    loss_function=loss_function,
    optimizer=optimizer,
    output_path="output"
)

# %% Load best weights (and biases) after iterations
model = model
model.load_state_dict(torch.load("output/best_weights.pth"))
model.eval()
out = model(torch.tensor(test_count_matrix.values).to(device))

# %% Evaluate model loaded with best weights (and biases)
predictions_list = []

# %% Convert PyTorch Tensor to Python list
y_preds = torch.squeeze(out)
predictions_list.append(y_preds.tolist())

## Use Complementary probabilities
predictions_flat = [1 - el for el in list(chain.from_iterable(predictions_list))]

# %% Calculate ROC AUC score
eval_roc_score = roc_auc_score(y_true=test_metadata["label"], y_score=predictions_flat)

# %% Check evaluation using complementary probabilities
model_train_eval(
    model=model,
    n_epochs=250,
    train_dataloader=train_dataloader,
    test_dataloader=test_dataloader,
    loss_function=loss_function,
    optimizer=optimizer,
    output_path="output",
    complementary_prob=True
)
