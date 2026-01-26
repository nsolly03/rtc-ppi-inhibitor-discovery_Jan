## Receptor Preparation - Step A Complete

### Script: 02b_prepare_receptors_simple.py

**Purpose:** Convert protein structures to PDBQT format for docking

**Usage:**
```bash
python scripts/02b_prepare_receptors_simple.py
```

**What it does:**
1. Extracts relevant protein chains
2. Converts PDB to PDBQT format
3. Assigns basic AutoDock atom types
4. Calculates center of mass for grid box placement

**Output files:**
- `7DFG_receptor.pdbqt` - Main target (RdRp complex)
- `6W4H_receptor.pdbqt` - NSP10-NSP16 interface
- `7EDI_receptor.pdbqt` - NSP10-NSP14 interface

**Note:** For production docking on CECI HPC, VirtualFlow will optimize these structures (add hydrogens, refine charges).

---# Complete Setup Guide: RTC-PPI Inhibitor Discovery Pipeline

**Author:** Olivier Nsekuye  
**Date:** January 2025  
**Lab:** GIGA-VIN, University of Liège  
**Project:** FRIA Doctoral Research - Pan-coronavirus inhibitor discovery

---

## Table of Contents

1. [Overview](#overview)
2. [System Requirements](#system-requirements)
3. [Software Installation](#software-installation)
4. [Environment Setup](#environment-setup)
5. [VS Code Configuration](#vs-code-configuration)
6. [Project Structure](#project-structure)
7. [GitHub Setup](#github-setup)
8. [Testing the Environment](#testing-the-environment)
9. [Troubleshooting](#troubleshooting)
10. [Daily Workflow](#daily-workflow)

---

## Overview

This guide documents the complete setup process for the RTC-PPI inhibitor discovery computational pipeline. The setup enables:

- **Virtual screening** of billions of compounds
- **Ligand preparation** and filtering
- **Receptor preparation** for docking
- **Data analysis** and visualization
- **Version control** with Git/GitHub
- **Reproducible research** environment

**Important:** Local docking is not set up on Mac (Apple Silicon limitations). Ultra-large-scale docking will be performed on CECI HPC cluster where VirtualFlow and AutoDock Vina are pre-installed.

---

## System Requirements

### Hardware
- **Computer:** MacBook Air (Apple Silicon M-series chip)
- **RAM:** 8 GB minimum (16 GB recommended)
- **Storage:** 50 GB free space minimum
- **Internet:** Required for package downloads and GitHub

### Operating System
- **OS:** macOS (Apple Silicon - ARM64 architecture)
- **Version:** macOS 11.0 or later

---

## Software Installation

### 1. Install Miniconda

**What is Miniconda?**  
Miniconda is a minimal Python distribution manager that allows you to create isolated environments with specific package versions. This prevents conflicts between different projects.

**Why Miniconda over Anaconda?**
- Lighter (~500 MB vs ~5 GB)
- Faster environment creation
- Only installs what you need
- Standard in computational chemistry
- Better for HPC deployment

**Installation Steps:**

1. Download Miniconda for Apple Silicon:
   - Go to: https://docs.conda.io/en/latest/miniconda.html
   - Download: **Miniconda3 macOS Apple M1 64-bit pkg**

2. Install:
   - Double-click the downloaded `.pkg` file
   - Follow installation prompts
   - Choose "Install for me only"

3. Initialize conda:
```bash
   /opt/miniconda3/bin/conda init zsh
```

4. Close and reopen Terminal

5. Verify installation:
```bash
   conda --version
```
   Expected output: `conda 23.x.x` or higher

6. Configure conda:
```bash
   # Add conda-forge channel (best for scientific packages)
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   
   # Disable auto-activation of base environment
   conda config --set auto_activate_base false
```

---

### 2. Install Git

**What is Git?**  
Git is a version control system that tracks changes to your code over time. It enables you to:
- Save snapshots (commits) of your work
- Revert to previous versions if needed
- Collaborate with others
- Backup code to cloud (GitHub)

**Why Git?**
- Essential for reproducible research
- Required for thesis/publication code sharing
- Industry standard for code management
- Integrates with GitHub for backup

**Installation Steps:**

1. Check if Git is already installed:
```bash
   git --version
```

2. If not installed, install via Homebrew:
```bash
   # Install Homebrew first (if needed)
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   
   # Install Git
   brew install git
```

3. Configure Git with your identity:
```bash
   git config --global user.name "Olivier Nsekuye"
   git config --global user.email "olivier.nsekuye@uliege.be"
   git config --global core.editor "nano"
```

4. Verify:
```bash
   git config --list
```

---

### 3. Install VS Code

**What is VS Code?**  
Visual Studio Code is a free, professional code editor that provides:
- Syntax highlighting for Python
- Integrated terminal
- Git integration
- Jupyter notebook support
- Extensions for scientific computing

**Why VS Code?**
- Industry standard for Python development
- Excellent Jupyter notebook integration
- Built-in Git support
- Free and open-source
- Active community and extensions

**Installation Steps:**

1. Download VS Code:
   - Go to: https://code.visualstudio.com/
   - Download for macOS (Universal or Apple Silicon)

2. Install:
   - Open the downloaded `.dmg` file
   - Drag VS Code to Applications folder

3. Add `code` command to Terminal:
   - Open VS Code
   - Press `Cmd+Shift+P`
   - Type: `Shell Command: Install 'code' command in PATH`
   - Press Enter

4. Verify:
```bash
   code --version
```

---

## Environment Setup

### 1. Create Conda Environment

**What is a Conda Environment?**  
An isolated Python environment with specific package versions. Different projects can have different environments without conflicts.

**Why Separate Environments?**
- Prevents package version conflicts
- Reproducible setup for collaborators
- Easy to recreate on different machines
- Can be exported and shared

**Create the rtc-chem environment:**
```bash
# Create environment with Python 3.10 and RDKit
conda create -n rtc-chem python=3.10 rdkit=2023.09 openbabel -c conda-forge

# Wait for package resolution and confirmation
# Type 'y' when prompted

# Activate the environment
conda activate rtc-chem
```

You should see `(rtc-chem)` at the start of your terminal prompt.

---

### 2. Install Python Packages

**Why These Packages?**

Each package serves a specific purpose in your pipeline:

| Package | Purpose | Why Essential |
|---------|---------|---------------|
| **meeko** | Ligand preparation for docking | Converts SMILES to PDBQT format |
| **biopython** | Protein structure handling | Read/write PDB files, extract chains |
| **pandas** | Data analysis | Handle screening results, IC50 data |
| **numpy** | Numerical computing | Array operations, calculations |
| **scipy** | Scientific computing | Statistical analysis, curve fitting |
| **scikit-learn** | Machine learning | Clustering, PCA, chemical space analysis |
| **matplotlib** | Plotting | Generate publication-quality figures |
| **seaborn** | Statistical visualization | Heatmaps, distribution plots |
| **py3Dmol** | 3D molecular visualization | View proteins and ligands in notebooks |
| **jupyterlab** | Notebook interface | Interactive analysis and visualization |
| **psutil** | System monitoring | Check RAM/CPU usage |
| **requests** | HTTP requests | Download PDB files, libraries |
| **gemmi** | Crystallographic data | Required by meeko |

**Installation command:**
```bash
pip install meeko biopython pandas numpy scipy scikit-learn matplotlib seaborn py3Dmol jupyterlab psutil requests gemmi
```

This will take 3-5 minutes.

---

### 3. Verify Installation

**Run environment check:**
```bash
cat > test_environment.py << 'EOF'
#!/usr/bin/env python3
"""Environment verification test"""

import sys

packages = {
    'RDKit': 'from rdkit import Chem',
    'Meeko': 'import meeko',
    'BioPython': 'from Bio.PDB import PDBParser',
    'Pandas': 'import pandas',
    'NumPy': 'import numpy',
    'SciPy': 'import scipy',
    'scikit-learn': 'from sklearn.cluster import KMeans',
    'Matplotlib': 'import matplotlib.pyplot',
    'Seaborn': 'import seaborn',
    'py3Dmol': 'import py3Dmol',
    'JupyterLab': 'import jupyterlab'
}

print("Testing environment...")
passed = 0
failed = 0

for name, import_cmd in packages.items():
    try:
        exec(import_cmd)
        print(f"✅ {name}")
        passed += 1
    except Exception as e:
        print(f"❌ {name}: {e}")
        failed += 1

print(f"\nResults: {passed}/{passed+failed} passed")

if failed == 0:
    print("✅ Environment ready for research!")
else:
    print(f"⚠️  Fix {failed} package(s)")
EOF

python test_environment.py
rm test_environment.py
```

Expected output: All packages should show ✅

---

## VS Code Configuration

### 1. Install Extensions

**Open VS Code and install these extensions:**

**Essential:**
- **Python** (Microsoft) - Python language support
- **Pylance** (Microsoft) - Fast Python IntelliSense
- **Jupyter** (Microsoft) - Notebook support

**Highly Recommended:**
- **GitLens** - Enhanced Git integration
- **Git Graph** - Visual Git history
- **Rainbow CSV** - CSV file visualization
- **Markdown All in One** - Better markdown editing

**How to install:**
1. Click Extensions icon (left sidebar)
2. Search for extension name
3. Click "Install"

---

### 2. Configure Python Interpreter

1. Open VS Code
2. Press `Cmd+Shift+P`
3. Type: `Python: Select Interpreter`
4. Choose: `rtc-chem` (your conda environment)

This ensures VS Code uses the correct Python with all your packages.

---

### 3. Configure Terminal

VS Code should automatically use your conda environment. Verify:

1. Open integrated terminal: `` Control+` ``
2. You should see `(rtc-chem)` in the prompt

If not:
```bash
conda activate rtc-chem
```

---

## Project Structure

### 1. Navigate to Project Location
```bash
cd ~/Desktop/Botanique/Project
```

---

### 2. Initialize Git Repository
```bash
git init
```

This creates a hidden `.git` folder that tracks all changes.

---

### 3. Create Folder Structure
```bash
mkdir -p data/libraries/zinc15/raw
mkdir -p data/libraries/zinc15/clean
mkdir -p data/libraries/enamine/raw
mkdir -p data/libraries/enamine/clean
mkdir -p data/targets
mkdir -p data/docking/receptors
mkdir -p data/docking/ligands_pdbqt
mkdir -p data/docking/results
mkdir -p data/assays
mkdir -p data/visualization
mkdir -p notebooks
mkdir -p scripts
mkdir -p results/docking
mkdir -p results/assays
mkdir -p results/ordering
mkdir -p results/figures
mkdir -p docs

touch data/libraries/zinc15/raw/.gitkeep
touch data/libraries/zinc15/clean/.gitkeep
touch data/libraries/enamine/raw/.gitkeep
touch data/libraries/enamine/clean/.gitkeep
touch data/targets/.gitkeep
touch data/docking/receptors/.gitkeep
touch data/docking/ligands_pdbqt/.gitkeep
touch data/docking/results/.gitkeep
touch data/assays/.gitkeep
touch data/visualization/.gitkeep
```

**Why this structure?**
```
Project/
├── data/              # All research data (excluded from Git for large files)
│   ├── libraries/     # Chemical compound libraries
│   ├── targets/       # Protein structures (PDB files)
│   ├── docking/       # Docking inputs and results
│   ├── assays/        # Experimental validation data
│   └── visualization/ # Molecular visualizations
├── scripts/           # Python automation scripts
├── notebooks/         # Jupyter analysis notebooks
├── results/           # Final results and figures (for publication)
└── docs/              # Documentation (like this file)
```

---

### 4. Create .gitignore

**What is .gitignore?**  
A file that tells Git which files to ignore (not track). Essential for:
- Large data files
- Temporary files
- System files
```bash
cat > .gitignore << 'EOF'
# Python
__pycache__/
*.pyc
*.pyo
.Python
*.egg-info/

# Jupyter
.ipynb_checkpoints

# Data files (too large for Git)
data/libraries/*/raw/*
data/libraries/*/clean/*
data/docking/results/*
*.pdbqt
*.sdf
*.mol2

# OS files
.DS_Store

# VS Code
.vscode/

# Keep .gitkeep files
!.gitkeep
EOF
```

---

### 5. Create README
```bash
cat > README.md << 'EOF'
# RTC-PPI Inhibitor Discovery Pipeline

**Author:** Olivier Nsekuye  
**Institution:** University of Liège, GIGA-VIN Laboratory  
**Supervisor:** Prof. Jean-Claude Twizere  
**Project:** FRIA 2025 Doctoral Research

## Overview

AI-guided discovery of novel pan-coronavirus inhibitors targeting protein-protein interactions within the SARS-CoV-2 replication-transcription complex (RTC).

## Targets

- NSP7–NSP12
- NSP9–NSP12
- NSP8–NSP7
- NSP12–NSP13
- NSP10–NSP14

## Setup

See `docs/SETUP_GUIDE.md` for complete installation instructions.

## Contact

Olivier Nsekuye - olivier.nsekuye@uliege.be
EOF
```

---

## GitHub Setup

### 1. Create GitHub Account

If you don't have one:
1. Go to https://github.com
2. Sign up with your university email

---

### 2. Create Repository

1. Log into GitHub
2. Click "+" → "New repository"
3. **Name:** `rtc-ppi-inhibitor-discovery_Jan`
4. **Description:** `PhD Research: Pan-coronavirus inhibitor discovery`
5. **Public** ✅
6. **DON'T check any boxes** (you already have files locally)
7. Click "Create repository"

---

### 3. Create Personal Access Token

**Why not use password?**  
GitHub requires tokens for security. Tokens can be revoked if compromised.

**Steps:**
1. Go to: https://github.com/settings/tokens
2. Click "Generate new token (classic)"
3. **Note:** `RTC-PPI Project Mac`
4. **Expiration:** 90 days
5. Check: ✅ `repo` (all sub-boxes)
6. Click "Generate token"
7. **Copy immediately** (you won't see it again!)

---

### 4. Connect Local Repository to GitHub
```bash
# Add remote
git remote add origin https://github.com/nsolly03/rtc-ppi-inhibitor-discovery_Jan.git

# Rename branch to main
git branch -M main

# Push to GitHub
git push -u origin main
```

When prompted:
- **Username:** nsolly03
- **Password:** Paste your token

---

### 5. Verify

Visit: https://github.com/nsolly03/rtc-ppi-inhibitor-discovery_Jan

You should see all your files!

---

## Testing the Environment

### 1. Test Python Environment
```bash
python --version
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "import meeko; print('Meeko OK')"
```

---

### 2. Test First Script

Download PDB structures:
```bash
python scripts/01_download_pdb_structures.py
```

Should download 5 protein structures to `data/targets/`

---

### 3. Test Git Workflow
```bash
bash scripts/git_backup.sh "Test commit"
```

Should push changes to GitHub.

---

## Troubleshooting

### Conda Command Not Found
```bash
# Reinitialize conda
/opt/miniconda3/bin/conda init zsh

# Close and reopen Terminal
```

---

### Package Import Errors
```bash
# Make sure environment is activated
conda activate rtc-chem

# Reinstall package
pip install --upgrade package-name
```

---

### Git Push Fails
```bash
# Pull first
git pull origin main --rebase

# Then push
git push origin main
```

---

### VS Code Can't Find Interpreter

1. Press `Cmd+Shift+P`
2. Type: `Python: Select Interpreter`
3. Click "Refresh"
4. Choose `rtc-chem`

---

## Daily Workflow

### Starting Work
```bash
# 1. Navigate to project
cd ~/Desktop/Botanique/Project

# 2. Activate environment
conda activate rtc-chem

# 3. Open VS Code
code .

# 4. Pull latest changes (if working on multiple computers)
git pull origin main
```

---

### During Work

1. Write code in `scripts/` or notebooks in `notebooks/`
2. Run and test
3. Save frequently (`Cmd+S`)

---

### Ending Work Session
```bash
# Backup to GitHub
bash scripts/git_backup.sh "Description of what you did"
```

---

## Notes for HPC Deployment

**For ultra-large-scale docking:**

Your local Mac is for:
- ✅ Code development
- ✅ Data preparation
- ✅ Result analysis
- ✅ Visualization

CECI HPC cluster is for:
- 🖥️ VirtualFlow screening (billions of molecules)
- 🖥️ AutoDock Vina docking
- 🖥️ Computationally intensive tasks

**To deploy on CECI:**
1. Export conda environment: `conda env export > environment.yml`
2. Transfer to HPC: `scp -r Project username@ceci-cluster:/path/`
3. Recreate environment: `conda env create -f environment.yml`

---

## References

- **RDKit:** https://www.rdkit.org/
- **Meeko:** https://github.com/forlilab/Meeko
- **BioPython:** https://biopython.org/
- **VirtualFlow:** https://virtual-flow.org/
- **CECI:** https://www.ceci-hpc.be/

---

## Changelog

**January 26, 2025**
- Initial environment setup
- Created project structure
- Configured Git and GitHub
- Downloaded initial PDB structures

---

**End of Setup Guide**
EOF
