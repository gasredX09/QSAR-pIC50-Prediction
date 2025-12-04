# QSAR pIC50 Prediction for Acetylcholinesterase

A comprehensive computational drug discovery pipeline that predicts molecular bioactivity (pIC50) for **Acetylcholinesterase**, a key drug target for Alzheimer's disease and related neurological conditions.

## Overview

This project implements a complete **QSAR (Quantitative Structure-Activity Relationship)** workflow using:
- **Data**: Bioactivity data from ChEMBL database
- **Descriptors**: Molecular fingerprints and descriptors via PaDEL-Descriptor
- **ML Model**: Random Forest regression with comparative analysis of alternative regressors
- **Web App**: Streamlit-based interactive prediction interface

## Features

✅ **Automated Data Pipeline**
- Fetch real bioactivity data from ChEMBL API
- Clean and curate datasets with IC50 values
- Classify bioactivity (active/inactive/intermediate)

✅ **Molecular Descriptor Generation**
- Integrates PaDEL-Descriptor for comprehensive fingerprinting
- Generates 1800+ molecular descriptors per molecule
- Handles SMILES string preprocessing

✅ **Machine Learning**
- Random Forest regression model with hyperparameter tuning
- Model comparison: RF, Gradient Boosting, SVR, and KNN
- Cross-validation and performance metrics (R², RMSE, MAE)

✅ **Interactive Web App**
- Streamlit UI for real-time predictions
- Batch processing of multiple molecules
- CSV/SMILES file uploads

✅ **Exploratory Data Analysis**
- Statistical summaries and distributions
- Feature importance ranking
- Data quality assessments

## Project Structure

```
QSAR-pIC50-Prediction/
├── project_main.py              # Main orchestration script
├── requirements.txt             # Python dependencies
├── README.md                    # This file
│
├── app/                         # Streamlit web application
│   ├── app.py                   # Main Streamlit app
│   ├── logo.png                 # App branding
│   ├── temp/                    # Temporary files (SMILES, descriptors)
│   └── PaDEL-Descriptor/        # PaDEL descriptor tool
│       ├── PaDEL-Descriptor.jar
│       ├── descriptors.xml      # Descriptor configuration
│       ├── lib/                 # Dependencies
│       └── license/
│
├── src/                         # Core OOP modules
│   ├── activity_model.py        # RandomForest model & training
│   ├── dataset_builder.py       # ChEMBL API client & data cleaning
│   ├── descriptor_engine.py     # PaDEL orchestration
│   ├── eda_engine.py            # Statistical analysis
│   ├── model_comparator.py      # Multi-model benchmark
│   └── molecule.py              # Molecule abstraction layer
│
├── data/
│   ├── raw/                     # Original bioactivity CSV from ChEMBL
│   ├── processed/               # Cleaned datasets
│   │   └── acetylcholinesterase_bioactivity.csv
│   └── descriptors/
│       ├── molecules.smi        # SMILES strings
│       ├── descriptors.csv      # Computed descriptors
│       └── final_dataset.csv    # Merged dataset (bioactivity + descriptors)
│
├── notebooks/                   # Jupyter analysis & prototyping
│   ├── 01_data_preprocessing.ipynb
│   ├── 02_exploratory_data_analysis.ipynb
│   ├── 03_descriptor_generation.ipynb
│   ├── 04_model_training.ipynb
│   └── 05_compare_regressors.ipynb
│
├── model/                       # Trained model artifacts
│   └── final_model.pkl          # Pickled RandomForest + feature names
│
└── results/                     # Output artifacts
    ├── eda/
    │   └── summary_stats.csv
    └── comparison/
        └── regressor_comparison.csv
```

## Installation

### Prerequisites
- Python 3.8+
- Java (required for PaDEL-Descriptor)
- Git

### Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/gasredX09/QSAR-pIC50-Prediction.git
   cd QSAR-pIC50-Prediction
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On macOS/Linux
   # OR
   venv\Scripts\activate  # On Windows
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Verify Java installation:**
   ```bash
   java -version
   ```

## Usage

### Option 1: Full Pipeline (Recommended for First Run)

Run the complete data → descriptors → model → comparison workflow:

```bash
python project_main.py
```

This will:
1. Fetch bioactivity data for Acetylcholinesterase from ChEMBL
2. Generate molecular descriptors using PaDEL
3. Train a Random Forest model
4. Compare with other regression algorithms
5. Save the trained model to `model/final_model.pkl`

### Option 2: Interactive Web Application

Launch the Streamlit app for real-time predictions:

```bash
streamlit run app/app.py
```

Then:
1. Upload a SMILES file (`.smi`, `.txt`, or `.csv`)
2. Click "Predict"
3. View predicted pIC50 values and download results

**Input formats:**
- `.smi` or `.txt`: One SMILES string per line
- `.csv`: Spreadsheet with a `canonical_smiles` or first column containing SMILES

### Option 3: Jupyter Notebooks

Explore the analysis step-by-step:

```bash
jupyter notebook notebooks/
```

Run notebooks in order:
1. `01_data_preprocessing.ipynb` — ChEMBL data fetch & cleaning
2. `02_exploratory_data_analysis.ipynb` — Statistical analysis
3. `03_descriptor_generation.ipynb` — PaDEL integration
4. `04_model_training.ipynb` — Model training & evaluation
5. `05_compare_regressors.ipynb` — Algorithm benchmarking

## Key Modules

### `src/dataset_builder.py`
- Fetches targets from ChEMBL by name
- Retrieves IC50 bioactivity data
- Cleans and filters for human single-protein targets
- Builds Molecule objects with pIC50 calculations

### `src/activity_model.py`
- Trains Random Forest regression models
- Handles feature selection and preprocessing
- Saves/loads pickled models with feature metadata
- Predicts pIC50 for new molecules

### `src/descriptor_engine.py`
- Orchestrates PaDEL-Descriptor execution
- Converts SMILES to molecular descriptors
- Standardizes and validates descriptor outputs

### `src/eda_engine.py`
- Computes statistical summaries
- Generates visualization-ready data
- Exports summary statistics to CSV

### `src/model_comparator.py`
- Benchmarks Random Forest, Gradient Boosting, SVR, KNN
- Cross-validates and reports metrics
- Exports comparison results

### `src/molecule.py`
- Core Molecule class with SMILES/ChEMBL metadata
- Computes pIC50 from IC50 values
- Serializes to dictionary for CSV export

## Data Sources

- **ChEMBL Database**: Real bioactivity data for Acetylcholinesterase
  - Target: CHEMBL220 (Human Acetylcholinesterase)
  - Bioactivity type: IC50 (nM)
  - Filtering: Human + Single Protein targets only

## Model Performance

The Random Forest model is trained on 200+ molecules with IC50 bioactivity data:

- **R² Score**: ~0.85 (typical)
- **RMSE**: ~0.3–0.4 log units
- **MAE**: ~0.2–0.3 log units

(Exact metrics vary based on dataset size and descriptor selection.)

## Requirements

Key dependencies are managed in `requirements.txt`:

- **pandas** — Data manipulation
- **numpy** — Numerical computing
- **scikit-learn** — ML algorithms & metrics
- **matplotlib, seaborn** — Visualization
- **chembl_webresource_client** — ChEMBL API access
- **rdkit** — Cheminformatics utilities
- **streamlit** — Web app framework
- **Pillow** — Image handling
- **jupyter** — Interactive notebooks (optional)

For a full list and pinned versions, see `requirements.txt`.

## Contributing

Contributions are welcome! To extend the project:

1. Add new regression algorithms to `src/model_comparator.py`
2. Implement additional descriptor preprocessing in `src/descriptor_engine.py`
3. Create new analysis notebooks in `notebooks/`
4. Improve the Streamlit UI in `app/app.py`

## Troubleshooting

### "PaDEL JAR not found"
Ensure `app/PaDEL-Descriptor/PaDEL-Descriptor.jar` exists or update the `PADEL_JAR` path in `app/app.py`.

### "Java not found"
Install Java:
```bash
# macOS
brew install openjdk

# Ubuntu/Debian
sudo apt-get install openjdk-11-jdk

# Windows
# Download from https://www.oracle.com/java/technologies/downloads/
```

### ChEMBL API timeout
Increase the `limit` parameter in `DatasetBuilder` or retry the pipeline. The API is rate-limited.

### Memory issues during descriptor generation
Reduce batch size or process molecules in smaller chunks. Edit `DescriptorEngine` to split SMILES files.

## License

This project is provided as-is for educational and research purposes. See `LICENSE` for details.

## Citation & Credits

- **Original concept**: Data Professor (YouTube QSAR series)
- **Implementation & OOP refactoring**: Aryan
- **Descriptor generation**: PaDEL-Descriptor (Yap et al.)
- **Bioactivity data**: ChEMBL (European Bioinformatics Institute)

## References

1. ChEMBL Database: https://www.ebi.ac.uk/chembl/
2. PaDEL-Descriptor: https://padel.nus.edu.sg/software/padeldescriptor/
3. RDKit: https://www.rdkit.org/
4. Scikit-learn: https://scikit-learn.org/

## Contact & Support

For questions or issues, please open a GitHub Issue or contact the repository maintainer.

---

**Last Updated**: December 2024  
**Status**: Active  
**Python Version**: 3.8+
