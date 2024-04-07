# %%
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from hyperopt import STATUS_OK, fmin, hp, tpe
from hyperopt.pyll import scope
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier, plot_importance

# %% Set up paths
ROOT_DIR = Path.cwd()
TRAIN_METADATA_PATH = ROOT_DIR / "intermediate/train_metadata.tsv"
TEST_METADATA_PATH = ROOT_DIR / "intermediate/test_metadata_filtered.tsv"
TRAIN_COUNT_MATRIX_PATH = (
    ROOT_DIR / "intermediate/train_counts_transformed_scaled.tsv.gz"
)
TEST_COUNT_MATRIX_PATH = ROOT_DIR / "intermediate/test_counts_transformed_scaled.tsv.gz"
FINAL_GENE_LIST_PATH = ROOT_DIR / "intermediate/final_gene_list.tsv"
XGBOOST_DIR = ROOT_DIR / "output" / "xgboost"

# %% Import refined gene list
gene_list = pd.read_csv(FINAL_GENE_LIST_PATH, sep="\t")
gene_list = gene_list["index"].to_list()

# %% Import train dataset
train_count_matrix = pd.read_csv(TRAIN_COUNT_MATRIX_PATH, sep="\t", index_col=0)
train_count_matrix = train_count_matrix[train_count_matrix.index.isin(gene_list)]
train_count_matrix = train_count_matrix.sort_index()
train_count_matrix = train_count_matrix.transpose()

train_metadata = pd.read_csv(TRAIN_METADATA_PATH, sep="\t")

# %% Import test dataset
test_count_matrix = pd.read_csv(TEST_COUNT_MATRIX_PATH, sep="\t", index_col=0)
test_count_matrix = test_count_matrix[test_count_matrix.index.isin(gene_list)]
test_count_matrix = test_count_matrix.sort_index()
test_count_matrix = test_count_matrix.transpose()

test_metadata = pd.read_csv(TEST_METADATA_PATH, sep="\t")

# %% Hyperparameter space for Bayesian optimization. Based on Rith Pansangas' tutorial available at https://medium.com/@rithpansanga/optimizing-xgboost-a-guide-to-hyperparameter-tuning-77b6e48e289d
space = {
    "max_depth": scope.int(hp.quniform("max_depth", 3, 7, 1)),
    "min_child_weight": scope.int(hp.quniform("min_child_weight", 1, 10, 1)),
    "subsample": hp.uniform("subsample", 0.5, 1),
    "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1),
    "learning_rate": hp.loguniform("learning_rate", -5, -2),
    "gamma": hp.uniform("gamma", 0.5, 5),
    "objective": "binary:logistic",
    "n_estimators": 1000,
    "early_stopping_rounds": 50,
    "eval_metric": "auc",
}


# %% Define function for score minimization. Here I set 1 - ROC_AUC so hyperopt can provide the hyperparameters associated with the maximum ROC AUC
def objective(params):
    xgb_model = XGBClassifier(**params)

    xgb_model.fit(
        train_count_matrix,
        train_metadata["label"],
        eval_set=[
            (train_count_matrix, train_metadata["label"]),
            (test_count_matrix, test_metadata["label"]),
        ],
    )

    y_pred = xgb_model.predict_proba(test_count_matrix)

    score = 1 - roc_auc_score(y_true=test_metadata["label"], y_score=y_pred[:, 1])

    return {"loss": score, "status": STATUS_OK}


# %% Perform Bayesian optimization
best_params = fmin(objective, space, algo=tpe.suggest, max_evals=100)

# %% Save best params dictonary as a JSON file
with open(str(XGBOOST_DIR / "best_params.json"), "w") as file:
    json.dump(best_params, file)

# %% Repeat model fit, this time with the best hyperparameters
xgb_model = XGBClassifier(
    **best_params, n_estimators=1000, early_stopping_rounds=50, eval_metric=["auc"]
)

xgb_model.fit(
    train_count_matrix,
    train_metadata["label"],
    eval_set=[
        (train_count_matrix, train_metadata["label"]),
        (test_count_matrix, test_metadata["label"]),
    ],
)

y_pred = xgb_model.predict_proba(test_count_matrix)

roc_auc = roc_auc_score(y_true=test_metadata["label"], y_score=y_pred[:, 1])

with open(XGBOOST_DIR / "roc_auc.txt", "w") as file:
    file.write(f"{str(roc_auc)}\n")

# %% Collect importance metrics into a dataframe
feature_importance_df = (
    pd.DataFrame(
        [
            xgb_model.get_booster().get_score(importance_type="gain"),
            xgb_model.get_booster().get_fscore(),
        ],
        index=["gain", "weight"],
    )
    .transpose()
    .reset_index()
    .rename(columns={"index": "gene"})
    .sort_values("weight", ascending=False)
)

feature_importance_df.to_csv(
    XGBOOST_DIR / "feature_importance.tsv", sep="\t", index=False
)

# %% Save importance (weight metric) plot with top 15 genes
plot_importance(xgb_model, max_num_features=15)
plt.savefig(str(XGBOOST_DIR / "top15_feature_importance_weight.png"))

# %% Save model as JSON file
xgb_model.save_model(str(XGBOOST_DIR / "xgb_model.json"))

# %%
