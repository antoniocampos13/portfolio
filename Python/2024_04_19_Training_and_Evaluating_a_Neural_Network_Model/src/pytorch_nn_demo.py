# %% Python 3.10.4 | Pytorch 1.12.1 | CUDA 12.4.1
import copy
from collections.abc import Callable
from itertools import chain
from pathlib import Path
from typing import Type

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score
from torch.utils.data import DataLoader

# %% Set computation device (GPU or CPU) depending on available hardware
if torch.cuda.is_available():
    device = torch.device("cuda:0")
else:
    device = torch.device("cpu")


class NeuralNetwork(nn.Module):
    """A class to perform the forward pass calculations of a neural network. Inherits from the Module class from PyTorch."""

    def __init__(
        self,
        n_features: int,
        n_neurons: int,
        n_output: int,
        activation_function: Callable[..., torch.FloatTensor],
        dropout: Callable[[torch.FloatTensor, float, bool, bool], torch.FloatTensor],
    ):
        """Constructs all the necessary attributes for the neural network object forward pass calculations. Intended for binary classification problems.

        Args:
            n_features (int): The number of features of the dataset (columns, independent variables).
            n_neurons (int): The number of neurons from the hidden layers.
            n_output (int): The number of classes of the output variable.
            activation_function (Callable[..., torch.FloatTensor]): The desired activation function of the hidden layers.
            dropout (Callable[[torch.FloatTensor, float, bool, bool], torch.FloatTensor]): The desire dropout rate (probability) of neurons.
        """
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        self.n_features = n_features
        self.n_neurons = n_neurons
        self.n_output = n_output
        self.activation_function = activation_function
        self.sigmoid = nn.Sigmoid()

        self.layer1 = nn.Linear(n_features, n_neurons)
        self.layer2 = nn.Linear(n_neurons, n_neurons)
        self.layer3 = nn.Linear(n_neurons, n_output)
        self.double()

    def forward(self, input: torch.FloatTensor) -> torch.FloatTensor:
        """_summary_

        Args:
            input (torch.FloatTensor): The dataset in PyTorch Tensor format.

        Returns:
            torch.FloatTensor: A PyTorch Tensor of the gradients of each element of the input Tensor.
        """
        fwd = self.activation_function(self.layer1(input))
        fwd = self.dropout(fwd)
        fwd = self.activation_function(self.layer2(fwd))
        fwd = self.dropout(fwd)
        fwd = self.sigmoid(self.layer3(fwd))
        return fwd


def train(
    model: Type[NeuralNetwork],
    train_dataloader: Type[DataLoader],
    loss_function: Callable[..., float],
    optimizer: Callable[..., torch.FloatTensor],
) -> tuple[float, torch.FloatTensor]:
    """A train function for a neural network intended for binary classification problems. It adds weights (and biases) optimizer and loss function to act upon the forward pass calculations.

    Args:
        model (Type[NeuralNetwork]): The neural network model.
        train_dataloader (Type[DataLoader]): A PyTorch DataLoader object containing a TensorDataset of the training data.
        loss_function (Callable[..., float]): The desired loss function to be applied on the results of the forward pass.
        optimizer (Callable[..., torch.FloatTensor]): The desired weight (and biases) optimizer.

    Returns:
        tuple[float, torch.FloatTensor]: A tuple containing the loss of a single forward pass and the calculated weights (and biases).
    """

    model.train()
    running_loss = 0.0
    for batch in train_dataloader:
        features, labels = batch
        features = features.to(device)
        labels = labels.to(device)
        optimizer.zero_grad()
        outputs = model(features)
        loss = loss_function(torch.squeeze(outputs), labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item() * labels.size(dim=0)

    train_loss = running_loss / len(train_dataloader.dataset)

    model_weights = copy.deepcopy(model.state_dict())

    return train_loss, model_weights


def test(
    model: Type[NeuralNetwork],
    model_weights: torch.FloatTensor,
    test_dataloader: Type[DataLoader],
    loss_function: Callable[..., float],
    complementary_prob: bool = False,
) -> tuple[float, float]:
    """A test function to evaluate a trained neural network intended for binary classification problems. Loads model weights and calculates ROC score by comparing actual output labels with predicted values (y_pred). The complementary probabilities (1 - y_pred) of the predicted labels are used if desired (useful if the general trend of ROC scores are below 0.50, meaning the model is actually predicting the negative class instead of the positive class).

    Args:
        model (Type[NeuralNetwork]): The neural network model.
        model_weights (torch.FloatTensor): The weights (and biases) calculated and optimized during the backward pass.
        test_dataloader (Type[DataLoader]):  A PyTorch DataLoader object containing a TensorDataset of the testing data.
        loss_function (Callable[..., float]): The desired loss function to be applied on the results of the forward pass.
        complementary_prob (bool, optional): Use complementary probabilities of the predicted labels instead of the default probabilities. Defaults to False.

    Returns:
        tuple[float, float]: A tuple containing the loss of a single forward pass and the calculated ROC score.
    """

    model.load_state_dict(model_weights)

    model.eval()
    running_loss = 0.0
    labels_list = []
    predictions_list = []
    with torch.no_grad():
        for batch in test_dataloader:
            features, labels = batch
            labels_list.append(labels.tolist())
            features = features.to(device)
            labels = labels.to(device)
            outputs = model(features)

            y_preds = torch.squeeze(outputs)
            predictions_list.append(y_preds.tolist())

            loss = loss_function(torch.squeeze(outputs), labels)
            running_loss += loss.item() * labels.size(dim=0)

        test_loss = running_loss / len(test_dataloader.dataset)

    labels_flat = list(chain.from_iterable(labels_list))

    if complementary_prob:
        predictions_flat = [
            1 - el for el in list(chain.from_iterable(predictions_list))
        ]
    else:
        predictions_flat = list(chain.from_iterable(predictions_list))

    roc_auc = roc_auc_score(y_true=labels_flat, y_score=predictions_flat)

    return test_loss, roc_auc


def model_train_eval(
    model: Type[NeuralNetwork],
    n_epochs: int,
    train_dataloader: Type[DataLoader],
    test_dataloader: Type[DataLoader],
    loss_function: Callable[..., float],
    optimizer: Callable[..., torch.FloatTensor],
    output_path: str,
    complementary_prob: bool = False,
) -> None:
    """A function for training and evaluating neural network models trough several epochs.

    Args:
        model (Type[NeuralNetwork]): The neural network model.
        n_epochs (int): The desired number of epochs for the training and evaluation rounds.
        train_dataloader (Type[DataLoader]): A PyTorch DataLoader object containing a TensorDataset of the training data.
        test_dataloader (Type[DataLoader]): A PyTorch DataLoader object containing a TensorDataset of the testing data.
        loss_function (Callable[..., float]): The desired loss function to be applied on the results of the forward pass.
        optimizer (Callable[..., torch.FloatTensor]): The desired weight (and biases) optimizer.
        output_path (str): A local path where the best weights (and biases) will be saved with the PyTorch's model state dictionary format.
        complementary_prob (bool, optional): Use complementary probabilities of the predicted labels instead of the default probabilities. Defaults to False.
    """

    path = Path(output_path)

    if not path.exists():
        path.mkdir()

    best_roc_auc = -np.inf
    best_weights = None

    for epoch in range(n_epochs):
        train_loss, model_weights = train(
            model=model,
            train_dataloader=train_dataloader,
            loss_function=loss_function,
            optimizer=optimizer,
        )

        test_loss, roc_auc = test(
            model=model,
            model_weights=model_weights,
            test_dataloader=test_dataloader,
            loss_function=loss_function,
            complementary_prob=complementary_prob,
        )

        print(
            f"Epoch {epoch + 1}/{n_epochs}. Train loss: {train_loss}. Test loss: {test_loss}. ROC AUC: {roc_auc}"
        )

        if roc_auc > best_roc_auc:
            best_roc_auc = roc_auc
            best_weights = model_weights

            torch.save(best_weights, str(path / "best_weights.pth"))
            print(f"Current best epoch: {epoch + 1}")
