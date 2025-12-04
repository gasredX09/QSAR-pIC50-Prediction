import os
import pandas as pd
import subprocess


class DescriptorEngine:
    """
    Generates molecular descriptors using PaDEL-Descriptor.
    - Writes SMILES file (.smi)
    - Runs PaDEL via subprocess
    - Merges descriptors with pIC50 labels
    """

    def __init__(self, cleanCSV, padelJarPath, outputFolder="../data/descriptors"):
        self.cleanCSV = cleanCSV
        self.outputFolder = outputFolder
        self.padelJar = padelJarPath

        self.smilesPath = f"{self.outputFolder}/molecules.smi"
        self.descriptorCSV = f"{self.outputFolder}/descriptors.csv"
        self.finalCSV = f"{self.outputFolder}/final_dataset.csv"

        os.makedirs(self.outputFolder, exist_ok=True)
        self.df = None
        self.descDF = None

    # -------------------------------
    # LOAD CLEAN DATA
    # -------------------------------
    def loadData(self):
        self.df = pd.read_csv(self.cleanCSV)
        print(f"Loaded {len(self.df)} rows from {self.cleanCSV}")

    # -------------------------------
    # WRITE SMILES FILE (.smi)
    # -------------------------------
    def writeSmilesFile(self):
        """
        Writes SMILES and molecule IDs in the format:
            SMILES    ID
        """
        with open(self.smilesPath, "w") as f:
            for _, row in self.df.iterrows():
                smi = row["canonical_smiles"]
                cid = row["molecule_chembl_id"]
                f.write(f"{smi}\t{cid}\n")

        print(f"Wrote SMILES file: {self.smilesPath}")

    # -------------------------------
    # RUN PADEL DESCRIPTOR GENERATION
    # -------------------------------
    def runPadel(self):
        print("Running PaDEL...")

        command = [
            "java",
            "-jar",
            self.padelJar,
            "-dir",
            self.outputFolder,
            "-file",
            self.descriptorCSV,
            "-2d",
            "-fingerprints",
            "-threads",
            "4",
            "-descriptortypes",
            "../app/PaDEL-Descriptor/descriptors.xml",
        ]

        subprocess.run(command, check=True)
        print("PaDEL completed.")

    # -------------------------------
    # LOAD DESCRIPTORS
    # -------------------------------
    def loadDescriptors(self):
        self.descDF = pd.read_csv(self.descriptorCSV)
        print(f"Loaded descriptor matrix with {self.descDF.shape[1]} features.")

    # -------------------------------
    # MERGE DESCRIPTORS + pIC50
    # -------------------------------
    def mergeDescriptors(self):
        """
        Merge descriptor matrix (rows indexed by molecule ID)
        with the pIC50 labels from clean CSV.
        """
        merged = self.descDF.merge(
            self.df[["molecule_chembl_id", "pIC50"]],
            left_on="Name",
            right_on="molecule_chembl_id",
        )
        merged = merged.drop(columns=["Name"])
        merged.to_csv(self.finalCSV, index=False)

        print(f"Final ML-ready dataset saved to: {self.finalCSV}")
        return merged

    # -------------------------------
    # FULL PIPELINE
    # -------------------------------
    def run(self):
        self.loadData()
        self.writeSmilesFile()
        self.runPadel()
        self.loadDescriptors()
        return self.mergeDescriptors()
