# Import modules
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix, roc_curve, precision_recall_curve
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV, KFold, cross_val_score
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

# Set up constants
RANDOM_SEED = 42
SPLITS = 10

# Set up file paths
project_folder = Path().resolve().parent.parent
features_folder = project_folder / "data" / "processed"

# Import features
X_train = pd.read_pickle(features_folder / "X_train")
X_test = pd.read_pickle(features_folder / "X_test")
y_train = pd.read_pickle(features_folder / "y_train")
y_test = pd.read_pickle(features_folder / "y_test")

# Classification model options
models = [
    DummyClassifier,
    LogisticRegression,
    DecisionTreeClassifier,
    KNeighborsClassifier,
    GaussianNB,
    SVC,
    RandomForestClassifier,
    xgboost.XGBClassifier,
]

for model in models:
    cls = model()

    kfold = KFold(n_splits=SPLITS, random_state=RANDOM_SEED)

    s = cross_val_score(cls, X_train, y_train, scoring="f1", cv=kfold)

    print(f"{model.__name__:22}  F1 Score: " f"{s.mean():.3f} STD: {s.std():.2f}")

# Let's use logistic regression
estimator = LogisticRegression()

# Logistic regression hyperparameters grid
grid = {"penalty":["l1", "l2", "elasticnet","none"],
"dual":[False, True],
"C":np.logspace(-3,3,7)
}

# Fit model while searching best parameters
logreg_cv=GridSearchCV(estimator,grid,cv=SPLITS)
logreg_cv.fit(X_train,y_train)

# Summary of the best parameters
print("tuned hyperparameters :(best parameters) ",logreg_cv.best_params_)

# Fit model with the best parameters
logreg = LogisticRegression(**logreg_cv.best_params_)
logreg.fit(X_train, y_train)

# Predict labels of samples in the testing dataset
predictions = logreg.predict(X_test)
predicted_proba = logreg.predict_proba(X_test)

# Print confusion matrix
conf_matrix = pd.DataFrame(data=confusion_matrix(y_test, predictions).T,index=["Prediction: controls", "Prediction: cases"], columns=["Actual: controls", "Actual: cases"])
print(conf_matrix)

# Print model metrics
print(f"Accuracy: {accuracy_score(y_test, predictions):.3f}")
print(f"Precision: {precision_score(y_test, predictions):.3f}")
print(f"Recall: {recall_score(y_test, predictions):.3f}")
print(f"F1 Score: {f1_score(y_test, predictions):.3f}")

# Plot ROC curve with AUC
false_pos_rate, true_pos_rate, proba = roc_curve(y_test, predicted_proba[:, -1])
plt.figure()
plt.plot([0,1], [0,1], linestyle="--") # plot random curve
plt.plot(false_pos_rate, true_pos_rate, marker=".", label=f"AUC = {roc_auc_score(y_test, predicted_proba[:, -1]):.2f}")
plt.title("ROC Curve")
plt.ylabel("True Positive Rate")
plt.xlabel("False Positive Rate")
plt.legend(loc="lower right")

# Calculate range of F1 scores with corresponding thresholds
precision, recall, thresholds = precision_recall_curve(y_test, predicted_proba[:, -1])

f1_scores = 2 * recall * precision / (recall + precision)

# Get optimal threshold (the one associated with maximum F1 score) and make new predictions
optimal_proba_cutoff = thresholds[np.argmax(f1_scores)]

f1_predictions = [1 if i >= optimal_proba_cutoff else 0 for i in predicted_proba[:, -1]]

# Print updated confusion matrix
conf_matrix_th = pd.DataFrame(data=confusion_matrix(y_test, f1_predictions).T, index=["Prediction: controls", "Prediction: cases"], columns=["Actual: controls", "Actual: cases"])
print(conf_matrix_th)

# Print updated metrics
print(f"Accuracy Before: {accuracy_score(y_test, predictions):.3f} --> Now: {accuracy_score(y_test, f1_predictions):.3f}")
print(f"Precision Before: {precision_score(y_test, predictions):.3f} --> Now: {precision_score(y_test, f1_predictions):.3f}")
print(f"Recall Before: {recall_score(y_test, predictions):.3f} --> Now: {recall_score(y_test, f1_predictions):.3f}")
print(f"F1 Score Before: {f1_score(y_test, predictions):.3f} --> Now: {f1_score(y_test, f1_predictions):.3f}")

# Save (Serialize) model and optimal cutoff
pickle.dump(logreg, open(project_folder / "models" / "logreg", "wb"))
pickle.dump(optimal_proba_cutoff, open(project_folder / "models" / "logreg", "wb"))
