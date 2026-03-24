# 🧬 Protein Redox Analysis Pipeline (BioDC-based)

> ⚡ Designed for scalable redox analysis across protein variants in computational protein engineering workflows.

A computational pipeline for large-scale analysis of redox properties in protein systems, extending BioDC with automated batch processing and structured result extraction.

This project enables systematic evaluation of redox potentials, electron transfer energetics, and coupling across multiple protein variants — supporting protein engineering, bioelectronics, and computational biophysics applications.

---

## 🚀 Key Features

- ⚡ Automated execution of BioDC across multiple protein inputs  
- 🔄 Batch processing pipeline for scalable analysis  
- 📊 Structured extraction of redox potential, reorganisation energy, and electronic coupling  
- 📁 Organised output directories and reproducible workflows  
- 🧠 Designed for HPC and large dataset processing  

---

## 🧬 My Contribution

- Adapted and extended parts of the BioDC workflow to better support protein structure inputs, with a focus on haem-containing (heme) proteins and their structural preparation
- Built a pipeline (`run_pipeline.py`) to automate BioDC execution across multiple protein systems  
- Developed analysis scripts (`analysis.py`) to extract key redox and electron transfer metrics  
- Designed a reproducible workflow for structured computational experiments  
- Enabled comparative analysis across protein variants  

---
## 📁 Project Structure
protein-redox-pipeline/
├── src/
│ ├── run_pipeline.py
│ ├── analysis.py
├── biodc/ # Updated BioDC code (MIT Licensed)
├── configs/
├── examples/
├── results/
├── README.md
├── LICENSE
└── .gitignore

---

## ⚙️ Installation

Clone the repository:

```bash
git clone https://github.com/your-username/protein-redox-pipeline.git
cd protein-redox-pipeline

Install dependencies:
pip install -r requirements.txt

---

## ⚙️ How to Run

### Step 1: Run pipeline

python src/run_pipeline.py \
  --biodc-script biodc/V2.2/BioDCv2.py \
  --input-dir examples \
  --template-input configs/input.txt \
  --output-root results

### Step 2: Analyse results
python src/analysis.py \
  --results-root results/runs_YYYYMMDD_HHMMSS \
  --output-csv results/aggregated_results.csv

### 📊 Example Output
sample	redox_potential	reorganisation_energy	electronic_coupling
protein1	-0.12	0.45	0.003
protein2	-0.08	0.39	0.005

---

## 🧠 Applications
- Protein engineering and redox tuning
- Electron transfer modelling in multi-heme systems
- Bioelectronic material design
- Computational biophysics workflows
- AI-driven drug discovery and protein design pipelines

---

## ⚠️ License & Attribution
This project builds upon BioDC (Version 2.2), developed by Matthew J. Guberman-Pfeffer.

Copyright (c) 2024 Matthew J. Guberman-Pfeffer

BioDC is licensed under the MIT License.
A copy of the original license is included in this repository.

Modifications, pipeline development, and analysis scripts in this project are my own.

---

## 📚 Original BioDC Documentation

BioDC is a Python program that automates and accelerates the computation of redox potentials, cooperativities, and conductivities in (polymeric) multi-heme cytochromes.

### Core Capabilities
- Preparation of proteins containing one or more heme groups for molecular dynamics with the AMBER forcefield
- Estimation of reaction and reorganisation free energies and electronic couplings for electron transfer
- Computation of redox cooperativity across heme groups
- Modelling charge transport using diffusion and steady-state flux models

### References
Derrida, B. J. Stat. Phys. (1983)
Guberman-Pfeffer, M. J. J. Phys. Chem. B (2022)
Jansson, F. (2011)
Nenashev et al. Phys. Rev. B (2010)
Breuer et al. PNAS (2014)
Jiang et al. JACS (2017)
Jiang et al. J. Phys. Chem. Lett. (2020)




