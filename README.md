# LUAD Structural Consequence Pipeline

## 🧬 Project Overview
This project integrates **Genomic Differential Expression** with **Structural Biology** to identify druggable targets in Lung Adenocarcinoma (LUAD). 

Standard RNA-seq analysis identifies *which* genes are upregulated, but not *if* they are viable drug targets. This pipeline bridges that gap by fetching AlphaFold structures for top candidates and analyzing their biophysical stability.

## 🎯 Biological Goal
To filter "statistically significant" cancer genes based on protein structural quality, distinguishing between **Disordered/Fibrous candidates** (poor targets) and **Stable/Globular candidates** (ideal targets).

## 🛠 Tech Stack
* **R (DESeq2):** Differential expression analysis of TCGA-LUAD patient data.
* **Python (Pandas, Requests):** API integration with MyGene.info and AlphaFold EBI.
* **BioPython:** Structural analysis (pLDDT extraction, Instability Index calculation).
* **Seaborn/Matplotlib:** Visualization of the "Structural Volcano Plot."

## 📊 Key Findings
The pipeline analyzed ~60,000 genes and filtered down to 100 top candidates.

* **Top Candidate:** `PSAT1` (Log2FC: 3.5, Instability: 28.7, Confidence: 97%).
    * *Verdict:* High-priority metabolic target.
* **False Positive Identified:** `COL10A1` (Log2FC: 4.6, Confidence: 58%).
    * *Verdict:* Disordered matrix protein; likely poor drug target despite high expression.

## 📂 Repository Structure
* `scripts/`:
    * `01_differential_expression.R`: Normalizes and filters raw counts.
    * `02_fetch_structures.py`: Downloads PDB files from AlphaFold Database.
    * `03_analyze_structures.py`: Calculates Instability Index and pLDDT.
    * `04_visualize_results.py`: Merges data and plots the Structural Map.
* `data/`: Raw counts and downloaded PDB structures.
* `results/`: CSVs containing expression and structural metrics.
* `plots/`: Final visualizations.
