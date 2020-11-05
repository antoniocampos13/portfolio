# Import modules
from pathlib import Path

import janitor as jn
import pandas as pd
from sklearn import preprocessing, model_selection

# Set up constants
RANDOM_SEED=42
TEST_SIZE=0.30
INDEX=676
CORR_THRESHOLD=0.80

# Set up paths
project_folder = Path().resolve().parent.parent
CSV_FILE = project_folder / "data" / "interim" / "prostate_cancer_dataset_processed.csv"

# Load CSV into pandas.DataFrame
df = pd.read_csv(CSV_FILE)
orig_df = df

# One hot encoding for categorical variables
df = pd.get_dummies(df, drop_first=True)

# Split train and test datasets
X, y = jn.get_features_targets(df, target_columns="status")

X_train, X_test, y_train, y_test = model_selection.train_test_split(
    X, y, test_size=TEST_SIZE, stratify=y, random_state=RANDOM_SEED
)

# Summary of labels
print(f"Training Set has {sum(y_train)} Positive Labels (cases) and {len(y_train) - sum(y_train)} Negative Labels (controls)")
print(f"Test Set has {sum(y_test)} Positive Labels (cases) and {len(y_test) - sum(y_test)} Negative Labels (controls)")

# Get numeric columns
num_cols = list(X.columns)[:INDEX]

# Scale numeric columns
sca = preprocessing.StandardScaler()
X_train[num_cols] = sca.fit_transform(X_train[num_cols])
X_test[num_cols] = sca.transform(X_test[num_cols])

# Identify correlated features
correlated_features = set()
correlation_matrix = X_train.corr()

for i in range(len(correlation_matrix.columns)):
    for j in range(i):
        if abs(correlation_matrix.iloc[i, j]) > CORR_THRESHOLD:
            colname = correlation_matrix.columns[i]
            correlated_features.add(colname)

# Check the number correlated feature pairs
print(f"Number of correlated feature pairs: {len(correlated_features)}")

# Remove correlated features
X_train = X_train.drop(list(correlated_features),axis=1)
X_test = X_test.drop(list(correlated_features),axis=1)

# Serialize features
X_train.to_pickle(project_folder/"data"/"processed"/"X_train")
X_test.to_pickle(project_folder/"data"/"processed"/"X_test")
y_train.to_pickle(project_folder/"data"/"processed"/"y_train")
y_test.to_pickle(project_folder/"data"/"processed"/"y_test")
