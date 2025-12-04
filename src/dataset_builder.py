import pandas as pd
from chembl_webresource_client.new_client import new_client
from molecule import Molecule


class DatasetBuilder:
    """
    Handles Part 1:
    - Find target by name (correct search method)
    - Filter to human + single protein
    - Fetch activity data with LIMIT to avoid hangs
    - Clean dataset
    - Build Molecule objects
    - Save final CSV
    """

    def __init__(self, targetName, limit=200):
        self.targetName = targetName
        self.targetID = None
        self.limit = limit
        self.rawDF = None
        self.cleanDF = None
        self.molecules = []

    # -----------------------------------------------------------
    # STEP 1 – Find target by NAME using .search() (correct)
    # -----------------------------------------------------------
    def fetchTargetID(self):
        target = new_client.target
        query = target.search(self.targetName)
        df = pd.DataFrame.from_dict(query)

        if df.empty:
            raise ValueError(f"No targets found for '{self.targetName}'")

        # Filter Human + Single Protein
        df = df[df["organism"] == "Homo sapiens"]
        df = df[df["target_type"] == "SINGLE PROTEIN"]

        if df.empty:
            raise ValueError("No human single-protein targets found")

        self.targetID = df.iloc[0]["target_chembl_id"]
        print("Selected Target ID:", self.targetID)

    # -----------------------------------------------------------
    # STEP 2 – Fetch activity data SAFELY with LIMIT
    # -----------------------------------------------------------
    def fetchFromChEMBL(self):
        if self.targetID is None:
            self.fetchTargetID()

        activity = new_client.activity

        fields = [
            "molecule_chembl_id",
            "canonical_smiles",
            "standard_type",
            "standard_value",
            "standard_units",
        ]

        print("Fetching activity data...")

        # Apply both filters BEFORE slicing
        result = (
            activity.filter(target_chembl_id=self.targetID)
            .filter(standard_type="IC50")
            .only(fields)[: self.limit]
        )

        self.rawDF = pd.DataFrame.from_dict(list(result))
        print("Fetched:", len(self.rawDF), "rows")

    # -----------------------------------------------------------
    # STEP 3 – Clean data
    # -----------------------------------------------------------
    def cleanData(self):
        df = self.rawDF.copy()

        df = df[df["standard_type"] == "IC50"]
        df = df[df["standard_value"].notna()]

        df = df[
            [
                "molecule_chembl_id",
                "canonical_smiles",
                "standard_value",
                "standard_units",
            ]
        ]

        df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")

        df = df.drop_duplicates(subset="molecule_chembl_id")
        df = df.reset_index(drop=True)

        self.cleanDF = df

    # -----------------------------------------------------------
    # STEP 4 – Build Molecule objects
    # -----------------------------------------------------------
    def buildMoleculeObjects(self):
        for _, row in self.cleanDF.iterrows():
            m = Molecule(
                chemblID=row["molecule_chembl_id"],
                smiles=row["canonical_smiles"],
                ic50=row["standard_value"],
            )
            m.computePIC50()
            if m.pic50 is not None:
                self.molecules.append(m)

    def classifyBioactivity(self):
        """
        Adds a 'bioactivity_class' column using thresholds:
        - IC50 <= 1000 nM → active
        - IC50 >= 10000 nM → inactive
        - otherwise → intermediate
        """
        classes = []
        for val in self.cleanDF["standard_value"]:
            v = float(val)
            if v >= 10000:
                classes.append("inactive")
            elif v <= 1000:
                classes.append("active")
            else:
                classes.append("intermediate")

        self.cleanDF["bioactivity_class"] = classes

    # -----------------------------------------------------------
    # STEP 5 – Save CSV
    # -----------------------------------------------------------
    def saveCleanCSV(self, path):
        # Convert molecules to dict-based dataframe
        df = pd.DataFrame([m.toDict() for m in self.molecules])

        # Attach bioactivity_class column from cleanDF
        df["bioactivity_class"] = self.cleanDF["bioactivity_class"].values

        df.to_csv(path, index=False)
        return df

    # -----------------------------------------------------------
    # STEP 6 – Run full pipeline
    # -----------------------------------------------------------
    def runPipeline(self, outputCSV):
        self.fetchFromChEMBL()
        self.cleanData()
        self.classifyBioactivity()
        self.buildMoleculeObjects()
        return self.saveCleanCSV(outputCSV)
