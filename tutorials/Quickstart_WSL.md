# Quickstart Guide: SPAC Visualization on Windows with WSL
This guide provides detailed instructions for setting up and launching the Spacy Visualization dashboard on a Windows machine using WSL (Windows Subsystem for Linux). It supports local testing and development workflows including environment setup, dependency resolution, and launching the app for analysis.

## Prerequisites
Before starting, ensure you’ve installed the following:
- Git for Windows
- Miniconda (Python=3.9.13)
- Ubuntu on WSL
- Visual Studio Code (WSL extension recommended)
Optional: GitHub CLI or access tokens for authentication.

## Windows + WSL Setup Instructions
### 1. Setup your Workspace 

```bash
mkdir -p ~/SPAC_Workplace
cd ~/SPAC_Workplace
```

### 2. Cloning the SPAC Shiny repository

After installing Git for Windows, clone the main branch to ensure the latest version of the app is available. 

```bash
git clone https://github.com/Summer2025-SPAC/SPAC_Shiny 
cd SPAC_Shiny 
git checkout main
```

### 3. Creating the Conda environment

Install the dependencies from the environment.yml file. Ensure python=3.9.13

```bash
conda env create -f environment.yml
conda activate shiny
```

This will install all dependencies listed in the file and create the environment with the name specified inside. This will also activate the environment.

### 4. Open VS Code in WSL

```bash
code ./
```

### 5. Run the app.
Use this command to run the SPAC app in the local browser.

```bash
shiny run --reload app.py
```

Visit http://127.0.0.1:8000 in your browser. Use this sample [dataset](https://zenodo.org/records/6376767/files/healthy_lung_adata.h5ad?download=1) set to test out the app. 

