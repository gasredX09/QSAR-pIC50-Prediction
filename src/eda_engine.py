import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors
import os


class EDAEngine:
    """
    Handles exploratory data analysis (EDA) on the cleaned dataset.
    Computes Lipinski descriptors, summary stats, and generates plots.
    """

    def __init__(self, csvPath, resultsFolder="../results/eda"):
        self.csvPath = csvPath
        self.resultsFolder = resultsFolder
        self.df = None

    # -------------------------------
    # LOAD CLEAN DATA
    # -------------------------------
    def loadData(self):
        self.df = pd.read_csv(self.csvPath)
        print(f"Loaded {len(self.df)} rows from: {self.csvPath}")

    def ensureFolder(self):
        if not os.path.exists(self.resultsFolder):
            os.makedirs(self.resultsFolder)
            print(f"Created folder: {self.resultsFolder}")

    # -------------------------------
    # COMPUTE BASIC DESCRIPTORS
    # -------------------------------
    def computeLipinski(self):
        """
        Adds Lipinski descriptors:
        - MW (Molecular Weight)
        - LogP
        - HBD (H-bond donors)
        - HBA (H-bond acceptors)
        """

        mwList = []
        logpList = []
        hbdList = []
        hbaList = []

        for smi in self.df["canonical_smiles"]:
            mol = Chem.MolFromSmiles(smi)

            if mol is None:
                mwList.append(None)
                logpList.append(None)
                hbdList.append(None)
                hbaList.append(None)
                continue

            mwList.append(Descriptors.MolWt(mol))
            logpList.append(Descriptors.MolLogP(mol))
            hbdList.append(Descriptors.NumHDonors(mol))
            hbaList.append(Descriptors.NumHAcceptors(mol))

        self.df["MW"] = mwList
        self.df["LogP"] = logpList
        self.df["HBD"] = hbdList
        self.df["HBA"] = hbaList

        print("Computed Lipinski descriptors.")

    # -------------------------------
    # BASIC SUMMARY STATISTICS
    # -------------------------------
    def computeStats(self):
        stats = self.df.describe()
        stats.to_csv(f"{self.resultsFolder}/summary_stats.csv")
        print("Saved summary statistics.")
        return stats

    # -------------------------------
    # BOXPLOTS
    # -------------------------------
    def boxPlot(self, column):
        plt.figure(figsize=(6, 5))
        sns.boxplot(data=self.df, y=column)
        plt.title(f"{column} Distribution")
        plt.savefig(f"{self.resultsFolder}/{column}_boxplot.png", dpi=150)
        plt.close()
        print(f"Saved: {column}_boxplot.png")

    # -------------------------------
    # SCATTERPLOTS vs pIC50
    # -------------------------------
    def scatterPlot(self, x):
        plt.figure(figsize=(6, 5))
        sns.scatterplot(data=self.df, x=x, y="pIC50")
        plt.title(f"{x} vs pIC50")
        plt.savefig(f"{self.resultsFolder}/{x}_vs_pIC50.png", dpi=150)
        plt.close()
        print(f"Saved: {x}_vs_pIC50.png")

    # -------------------------------
    # RUN FULL EDA
    # -------------------------------
    def run(self):
        self.loadData()
        self.ensureFolder()
        self.computeLipinski()
        self.computeStats()

        for col in ["MW", "LogP", "HBD", "HBA"]:
            self.boxPlot(col)
            self.scatterPlot(col)

        return self.df
