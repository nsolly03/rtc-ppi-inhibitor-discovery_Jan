cat > docs/TROUBLESHOOTING.md << 'EOF'
# Troubleshooting Log

**Project:** RTC-PPI Inhibitor Discovery  
**Author:** Olivier Nsekuye  
**Started:** January 26, 2025

---

## Issues Encountered and Solutions

### 1. Conda Not Recognized in Terminal

**Date:** January 26, 2025  
**Issue:** After installing Miniconda, `conda` command not found

**Cause:** Terminal needed to be configured to find conda

**Solution:**
```bash
/opt/miniconda3/bin/conda init zsh
# Close and reopen terminal
```

**Prevention:** Always run `conda init` after fresh install

---

### 2. AutoDock Vina Installation Failure

**Date:** January 26, 2025  
**Issue:** `conda install autodock-vina` failed on Apple Silicon (ARM64)

**Error:** `PackagesNotFoundError: autodock-vina not available`

**Cause:** Vina not compiled for Apple Silicon architecture

**Solution:** 
- For local work: Use simple PDBQT conversion (sufficient for testing)
- For production: Use CECI HPC cluster (Vina pre-installed)

**Lesson:** Apple Silicon (M1/M2/M3) has limited support for some chemistry tools

---

### 3. Meeko Receptor Preparation Failed

**Date:** January 26, 2025  
**Issue:** `mk_prepare_receptor.py` failed with "No module named 'prody'"

**Solution:**
```bash
pip install prody
```

**Alternative:** Created `02b_prepare_receptors_simple.py` that works without meeko CLI

---

### 4. Missing NSP10 in 7EDI Structure

**Date:** January 26, 2025  
**Issue:** Visualization showed only NSP14 (purple), missing NSP10

**Cause:** Script only extracted Chain A, excluded Chain B (NSP10)

**Impact:** 
- Could not visualize NSP10-NSP14 interface
- Would have targeted wrong binding site
- Virtual screening would have failed

**Solution:**
```python
# In 02b_prepare_receptors_simple.py, changed:
'7EDI': {'chains': ['A']}  # Wrong - only NSP14

# To:
'7EDI': {'chains': ['A', 'B']}  # Correct - NSP14 + NSP10
```

**Verification:**
```bash
python -c "from Bio.PDB import PDBParser; parser = PDBParser(QUIET=True); structure = parser.get_structure('7EDI', 'data/docking/receptors/7EDI_extracted.pdb'); print([c.id for c in structure.get_chains()])"
```

Expected: `['A', 'B']`

**Lesson:** Always verify biological assembly contains all functional components

**Prevention:**
1. Check PDB database for biological assembly definition
2. Verify chain composition after extraction
3. Visualize structures before proceeding to docking

---

### 5. Jupyter Notebook JSON Parse Error

**Date:** January 26, 2025  
**Issue:** "Unexpected token 'c', 'cat > note'... is not valid JSON"

**Cause:** Attempted to open command text instead of running it in terminal

**Solution:** 
1. Commands starting with `cat >` must be run in **terminal**
2. Only the resulting `.ipynb` file should be opened in VS Code

**Prevention:** 
- Terminal commands go in terminal (bottom panel)
- Notebooks open via left sidebar or `code` command

---

### 6. Git Push Rejected - Remote Contains Work

**Date:** January 26, 2025  
**Issue:** `git push` rejected with "Updates were rejected because the remote contains work"

**Solution:**
```bash
git pull origin main --rebase
git push origin main
```

**Cause:** Changes made on GitHub (via web) not present locally

**Prevention:** Always `git pull` before starting work if using multiple machines

---

## Best Practices Developed

### 1. Structure Verification Workflow
```bash
# Always verify chains after extraction
python << 'EOF'
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure('check', 'path/to/file.pdb')
for chain in structure.get_chains():
    print(f"Chain {chain.id}: {len(list(chain.get_residues()))} residues")
EOF
```

### 2. Daily Work Routine
```bash
# Start of day
cd ~/Desktop/Botanique/Project
conda activate rtc-chem
git pull origin main

# End of day
bash scripts/git_backup.sh "Description of work done"
```

### 3. Before Running Scripts
1. Verify input files exist
2. Check expected output format
3. Run on small test case first
4. Verify output before proceeding

---

## Useful Debugging Commands

### Check Conda Environment
```bash
conda env list
conda list  # Show installed packages
which python  # Verify using correct Python
```

### Check Git Status
```bash
git status
git log --oneline -5
git remote -v
```

### Check File Structure
```bash
tree -L 2  # Overview
ls -lh data/docking/receptors/*.pdbqt  # Check specific files
```

### Python Environment Test
```bash
python -c "import rdkit; print('RDKit OK')"
python -c "import meeko; print('Meeko OK')"
python -c "from Bio.PDB import PDBParser; print('BioPython OK')"
```

---

## When to Update This Log

Add entry when you encounter:
- Errors that took >10 minutes to solve
- Non-obvious solutions
- System-specific issues (Mac vs Linux)
- Issues that might recur
- Important lessons learned

---

**Last Updated:** January 26, 2025
EOF

echo "✅ Troubleshooting log created!"
