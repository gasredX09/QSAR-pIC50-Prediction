import pandas as pd
import numpy as np
import os
import pickle
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error


class ActivityModel:
    """
    Trains a regression model to predict pIC50 from molecular descriptors.
    Handles:
    - train/test split
    - preprocessing (drop non-features, remove low-variance columns)
    - model fitting
    - evaluation
    - saving/loading model
    """

    def __init__(self):
        self.model = None
        self.featureNames = None

    # -------------------------------
    # LOAD DATA
    # -------------------------------
    def loadDataset(self, path):
        df = pd.read_csv(path)
        return df

    # -------------------------------
    # CLEAN FEATURES
    # -------------------------------
    def prepareXY(self, df):
        """
        Cleans the descriptor dataset:
        - drops non-numeric columns
        - removes inf and -inf
        - removes rows with NaN
        - removes zero-variance columns
        """

        df = df.copy()

        # Drop ID column if present
        if "molecule_chembl_id" in df.columns:
            df = df.drop(columns=["molecule_chembl_id"])

        # Ensure pIC50 exists
        if "pIC50" not in df.columns:
            raise ValueError("pIC50 column not found in dataset")

        # Separate target
        y = df["pIC50"]
        X = df.drop(columns=["pIC50"])

        # Convert all descriptor columns to numeric (force errors to NaN)
        X = X.apply(pd.to_numeric, errors="coerce")

        # Replace inf with NaN
        X = X.replace([np.inf, -np.inf], np.nan)

        # Drop any rows containing NaN
        X = X.dropna(axis=0)
        y = y.loc[X.index]  # align y with cleaned X

        # Remove zero-variance columns
        nunique = X.nunique()
        zeroVarCols = nunique[nunique == 1].index
        X = X.drop(columns=zeroVarCols)

        self.featureNames = X.columns.tolist()
        return X, y

    # -------------------------------
    # TRAIN MODEL
    # -------------------------------
    def train(self, X, y, testSize=0.2, randomState=42):
        Xtrain, Xtest, ytrain, ytest = train_test_split(
            X, y, test_size=testSize, random_state=randomState
        )

        # Baseline Random Forest
        self.model = RandomForestRegressor(
            n_estimators=500, max_depth=None, random_state=42, n_jobs=-1
        )

        self.model.fit(Xtrain, ytrain)

        # predictions
        preds = self.model.predict(Xtest)

        # metrics
        r2 = r2_score(ytest, preds)
        mse = mean_squared_error(ytest, preds)  # no squared kwarg
        rmse = mse**0.5  # manual sqrt
        mae = mean_absolute_error(ytest, preds)

        print("Model Performance:")
        print(f"R²:   {r2:.4f}")
        print(f"RMSE: {rmse:.4f}")
        print(f"MAE:  {mae:.4f}")

        return r2, rmse, mae

    # -------------------------------
    # SAVE MODEL
    # -------------------------------
    def saveModel(self, path="model.pkl"):
        # Create directory if needed
        folder = os.path.dirname(path)
        if folder != "" and not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)

        # Save model + features
        with open(path, "wb") as f:
            pickle.dump((self.model, self.featureNames), f)

        print(f"Saved model to {path}")

    # -------------------------------
    # LOAD MODEL
    # -------------------------------
    def loadModel(self, path="model.pkl"):
        with open(path, "rb") as f:
            self.model, self.featureNames = pickle.load(f)
        print(f"Loaded model from {path}")

    # -------------------------------
    # PREDICT
    # -------------------------------
    def predict(self, Xrow):
        """
        Xrow should be a pandas Series or a 1-row DataFrame
        containing the same columns as featureNames.
        """
        Xrow = Xrow[self.featureNames].values.reshape(1, -1)
        return self.model.predict(Xrow)[0]
