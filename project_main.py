"""
project_main.py

Entry point for the QSAR Computational Drug Discovery pipeline
(Acetylcholinesterase - pIC50 regression).

This script orchestrates the full pipeline:

1. Fetch and clean bioactivity data from ChEMBL
2. Run exploratory data analysis (EDA)
3. Generate molecular descriptors using PaDEL
4. Train a Random Forest regression model and save it
5. Compare multiple regressors and save comparison results
"""

import os
import sys
import pandas as pd


# ---------------------------------------------------------
# Resolve project root and add src/ to import path
# ---------------------------------------------------------
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(ROOT_DIR, "src")
sys.path.append(SRC_DIR)

from dataset_builder import DatasetBuilder
from eda_engine import EDAEngine
from descriptor_engine import DescriptorEngine
from activity_model import ActivityModel
from model_comparator import ModelComparator


# ---------------------------------------------------------
# Paths (relative to project root)
# ---------------------------------------------------------
DATA_PROCESSED = os.path.join("data", "processed")
DATA_DESCRIPTORS = os.path.join("data", "descriptors")
RESULTS_DIR = os.path.join("results")
RESULTS_EDA = os.path.join(RESULTS_DIR, "eda")
RESULTS_COMPARISON = os.path.join(RESULTS_DIR, "comparison")
MODEL_DIR = os.path.join("model")

os.makedirs(DATA_PROCESSED, exist_ok=True)
os.makedirs(DATA_DESCRIPTORS, exist_ok=True)
os.makedirs(RESULTS_EDA, exist_ok=True)
os.makedirs(RESULTS_COMPARISON, exist_ok=True)
os.makedirs(MODEL_DIR, exist_ok=True)

# filenames
TARGET_NAME = "acetylcholinesterase"
CLEAN_CSV = os.path.join(DATA_PROCESSED, f"{TARGET_NAME}_bioactivity.csv")
FINAL_DATASET = os.path.join(DATA_DESCRIPTORS, "final_dataset.csv")
MODEL_PATH = os.path.join(MODEL_DIR, "final_model.pkl")

# PaDEL jar inside app/PaDEL-Descriptor
PADEL_JAR = os.path.join("app", "PaDEL-Descriptor", "PaDEL-Descriptor.jar")
PADEL_OUTPUT_FOLDER = DATA_DESCRIPTORS


# ---------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------


def step1_build_dataset():
    print("\n[STEP 1] Fetching and cleaning ChEMBL bioactivity data...")
    builder = DatasetBuilder(TARGET_NAME)
    clean_df = builder.runPipeline(CLEAN_CSV)
    print(f"  -> Saved clean bioactivity CSV to {CLEAN_CSV}")
    print(f"  -> Shape: {clean_df.shape}")
    return clean_df


def step2_run_eda():
    print("\n[STEP 2] Running exploratory data analysis (EDA)...")
    try:
        eda = EDAEngine(cleanCSV=CLEAN_CSV, outputFolder=RESULTS_EDA)
        eda.run()
        print(f"  -> EDA figures and summaries saved to {RESULTS_EDA}")
    except Exception as e:
        print(f"  !! EDA step failed (continuing pipeline): {e}")


def step3_generate_descriptors():
    print("\n[STEP 3] Generating molecular descriptors with PaDEL...")
    engine = DescriptorEngine(
        cleanCSV=CLEAN_CSV,
        padelJarPath=PADEL_JAR,
        outputFolder=PADEL_OUTPUT_FOLDER,
    )
    final_df = engine.run()
    print(f"  -> Final descriptor dataset saved to {FINAL_DATASET}")
    print(f"  -> Shape: {final_df.shape}")
    return final_df


def step4_train_model():
    print("\n[STEP 4] Training Random Forest regression model...")
    df = pd.read_csv(FINAL_DATASET)
    am = ActivityModel()
    X, y = am.prepareXY(df)
    am.train(X, y)
    am.saveModel(MODEL_PATH)
    print(f"  -> Trained model saved to {MODEL_PATH}")


def step5_compare_regressors():
    print("\n[STEP 5] Comparing regressors on the same dataset...")
    mc = ModelComparator(FINAL_DATASET)
    results = mc.run()
    print("  -> Top models (head):")
    print(results.head())
    print(f"  -> Full comparison CSV saved under {RESULTS_COMPARISON}")


def main():
    print("========== QSAR Pipeline: Acetylcholinesterase ==========")
    print(f"Project root: {ROOT_DIR}")

    # 1. dataset
    step1_build_dataset()

    # 2. eda
    step2_run_eda()

    # 3. descriptors
    step3_generate_descriptors()

    # 4. model training
    step4_train_model()

    # 5. model comparison
    step5_compare_regressors()

    print("\n✅ Pipeline complete.")
    print("You can now run the Streamlit app with:")
    print("    streamlit run app/app.py")


if __name__ == "__main__":
    main()
