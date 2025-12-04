import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

from sklearn.ensemble import (
    RandomForestRegressor,
    ExtraTreesRegressor,
    GradientBoostingRegressor,
)
from sklearn.tree import DecisionTreeRegressor


class ModelComparator:
    """
    Compares a small set of strong, tree-based regression models.
    Outputs:
      - performance table
      - CSV summary
      - RMSE plot
    """

    def __init__(self, datasetPath, resultsFolder="../results/comparison"):
        self.datasetPath = datasetPath
        self.resultsFolder = resultsFolder
        self.df = None
        self.X = None
        self.y = None
        self.results = None

    # -------------------------------
    # LOAD + CLEAN
    # -------------------------------
    def loadData(self):
        self.df = pd.read_csv(self.datasetPath)
        self.df.columns = self.df.columns.str.strip()

        print(f"Loaded {self.df.shape[0]} rows, {self.df.shape[1]} columns.")

        # Drop ID column
        if "molecule_chembl_id" in self.df.columns:
            self.df = self.df.drop(columns=["molecule_chembl_id"])

        # Ensure pIC50
        if "pIC50" not in self.df.columns:
            raise KeyError("pIC50 column missing!")

        # Remove zero-variance
        self.df = self.df.loc[:, self.df.nunique() > 1]

        # Replace inf values
        self.df = self.df.replace([float("inf"), float("-inf")], pd.NA)

        # Drop NaN rows
        self.df = self.df.dropna(axis=0).reset_index(drop=True)

        # Clip extreme values
        self.df = self.df.clip(lower=-1e6, upper=1e6)

        print("After cleaning:", self.df.shape)

        # Extract X, y
        self.y = self.df["pIC50"]
        self.X = self.df.drop(columns=["pIC50"])

        # Remove near-constant features (std < 0.01)
        lowVarMask = self.X.std() < 0.01
        self.X = self.X.loc[:, ~lowVarMask]

        print("After filtering:", self.X.shape, "features")

    # -------------------------------
    # COMPARE TREE MODELS
    # -------------------------------
    def runComparison(self, testSize=0.2, randomState=42):
        Xtrain, Xtest, ytrain, ytest = train_test_split(
            self.X, self.y, test_size=testSize, random_state=randomState
        )

        models = {
            "RandomForest": RandomForestRegressor(
                n_estimators=400, random_state=42, n_jobs=-1
            ),
            "ExtraTrees": ExtraTreesRegressor(
                n_estimators=400, random_state=42, n_jobs=-1
            ),
            "GradientBoosting": GradientBoostingRegressor(random_state=42),
            "DecisionTree": DecisionTreeRegressor(random_state=42),
        }

        results = []
        print("\nRunning model comparison...\n")

        for name, model in models.items():
            print(f"Training {name}...")
            model.fit(Xtrain, ytrain)
            preds = model.predict(Xtest)

            r2 = r2_score(ytest, preds)
            rmse = mean_squared_error(ytest, preds) ** 0.5
            mae = mean_absolute_error(ytest, preds)

            results.append([name, r2, rmse, mae])

        dfResults = pd.DataFrame(results, columns=["Model", "R2", "RMSE", "MAE"])
        dfResults = dfResults.sort_values("RMSE")

        self.results = dfResults
        return dfResults

    # -------------------------------
    # SAVE CSV
    # -------------------------------
    def saveResults(self):
        os.makedirs(self.resultsFolder, exist_ok=True)
        outPath = f"{self.resultsFolder}/regressor_comparison.csv"
        self.results.to_csv(outPath, index=False)
        print(f"\nSaved comparison table to {outPath}")
        return outPath

    # -------------------------------
    # PLOT
    # -------------------------------
    def plotTop(self):
        os.makedirs(self.resultsFolder, exist_ok=True)

        df = self.results.copy()

        plt.figure(figsize=(10, 5))
        plt.barh(df["Model"], df["RMSE"], color="skyblue")
        plt.xlabel("RMSE (lower is better)")
        plt.title("QSAR Model Comparison")
        plt.gca().invert_yaxis()

        outPath = f"{self.resultsFolder}/rmse_plot.png"
        plt.savefig(outPath, dpi=150)
        plt.close()

        print(f"Saved RMSE plot to {outPath}")
        return outPath

    # -------------------------------
    # RUN FULL PIPELINE
    # -------------------------------
    def run(self):
        self.loadData()
        results = self.runComparison()
        self.saveResults()
        self.plotTop()
        return results
