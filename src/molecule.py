class Molecule:
    """
    Represents a single molecule with its metadata, IC50 value,
    computed pIC50, and (later) molecular descriptors.
    """

    def __init__(self, chemblID, smiles, ic50):
        self.chemblID = chemblID
        self.smiles = smiles
        self.ic50 = ic50
        self.pic50 = None
        self.descriptors = None

    def computePIC50(self):
        """
        Converts IC50 (nM) to pIC50 using:
            pIC50 = -log10(IC50 in molar)
        """

        # Ensure IC50 is a float
        try:
            ic50_float = float(self.ic50)
        except:
            self.pic50 = None
            return

        capped = min(ic50_float, 1e8)  # cap at 100 million nM
        molar = capped / 1e9  # nM → molar

        if molar <= 0:
            self.pic50 = None
        else:
            import math

            self.pic50 = -math.log10(molar)

    def setDescriptors(self, descDict):
        """
        Stores descriptor dictionary for the molecule.
        """
        self.descriptors = descDict

    def toDict(self):
        """
        Returns a clean dictionary representation for CSV export.
        """
        return {
            "molecule_chembl_id": self.chemblID,
            "canonical_smiles": self.smiles,
            "standard_value_nM": self.ic50,
            "pIC50": self.pic50,
        }
