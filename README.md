# QSAR-DrugDiscovery

This repository contains a **Computational Drug Discovery project** focused on predicting the inhibitory activity of compounds using **Quantitative Structure–Activity Relationship (QSAR)** modeling. The workflow demonstrates how to collect biological activity data, compute molecular descriptors, and build predictive machine learning models for drug discovery.

---

## 📌 Project Overview

1. **Data Collection**
   - Downloaded biological activity data from the **ChEMBL database**.
   - Selected **Acetylcholinesterase inhibitors** as the target protein.

2. **Feature Engineering**
   - Converted compound SMILES strings into molecular descriptors.
   - Computed **Lipinski’s Rule of Five** descriptors (MW, LogP, H-bond donors/acceptors).
   - Generated extended descriptors using **PADEL-Descriptor**.

3. **Exploratory Data Analysis (EDA)**
   - Visualized active vs. inactive compounds using boxplots and scatterplots.
   - Explored trends in physicochemical properties.

4. **Model Building**
   - Prepared datasets (X = descriptors, Y = pIC50 values).
   - Built regression models using **Scikit-learn**.
   - Benchmarked multiple algorithms with **LazyPredict** to identify top-performing models.

---

## 🛠️ Tools and Technologies

- **Programming:** Python (Pandas, NumPy, Scikit-learn, Matplotlib, Seaborn)
- **Databases:** ChEMBL  
- **Descriptor Calculation:** PADEL-Descriptor  
- **ML Utilities:** LazyPredict

---

## 📂 Repository Structure

```
QSAR-DrugDiscovery/
│── data/                 # Raw and processed datasets
│── descriptors/          # Molecular descriptors generated
│── notebooks/            # Jupyter notebooks for each step
│── models/               # Trained regression models
│── results/              # Plots, benchmark tables, evaluation metrics
│── README.md             # Project documentation
```

---

## 🚀 How to Run

1. Clone this repository:  
   ```bash
   git clone https://github.com/<your-username>/QSAR-DrugDiscovery.git
   cd QSAR-DrugDiscovery
   ```

2. Install dependencies:  
   ```bash
   pip install -r requirements.txt
   ```

3. Run Jupyter notebooks step by step from the `notebooks/` folder.

---

## 📊 Results

- Successfully built QSAR models for predicting **pIC50 values of Acetylcholinesterase inhibitors**.  
- Identified key molecular descriptors influencing compound activity.  
- Compared regression models to highlight best-performing algorithms.

---

## 📖 References

- [ChEMBL Database](https://www.ebi.ac.uk/chembl/)  
- [PADEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/)  
- Relevant QSAR and computational drug discovery literature.

---

## 👤 Author

**Aryan Sharan Reddy Guda**  
Master of Science – Quantitative Biology & Bioinformatics, Carnegie Mellon University  
[LinkedIn](https://www.linkedin.com/in/gasredx09/) | [GitHub](https://github.com/gasredX09)
