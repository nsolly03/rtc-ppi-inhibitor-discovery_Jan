# Quick Reference Guide

## Auto-Backup Script
```bash
# Quick backup with custom message
bash scripts/git_backup.sh "Description of changes"

# Quick backup with automatic timestamp
bash scripts/git_backup.sh
```

Use this at the end of each work session to save your progress to GitHub.

---# Quick Reference Guide

## Conda Commands
```bash
# Activate environment
conda activate rtc-chem

# Deactivate
conda deactivate

# List environments
conda env list

# List installed packages
conda list

# Install package
conda install package-name

# Update package
conda update package-name
```

## Git Commands
```bash
# Check status
git status

# Add files
git add .

# Commit
git commit -m "message"

# Push to GitHub
git push origin main

# Pull from GitHub
git pull origin main

# View history
git log --oneline

# Quick backup
bash scripts/git_backup.sh "message"
```

## VS Code Shortcuts (Mac)

| Action | Shortcut |
|--------|----------|
| Command Palette | `Cmd+Shift+P` |
| Save | `Cmd+S` |
| Open Terminal | ``Control+` `` |
| Find | `Cmd+F` |
| Run Cell (Jupyter) | `Shift+Enter` |
| Comment Line | `Cmd+/` |

## Project Navigation
```bash
# Go to project
cd ~/Desktop/Botanique/Project

# Open in VS Code
code .

# List structure
tree -L 2
# or
ls -R
```

## Python Quick Tests
```bash
# Test imports
python -c "from rdkit import Chem; print('OK')"

# Run script
python scripts/script_name.py

# Start Jupyter
jupyter lab
```
---

## Common Issues & Quick Fixes

### Conda not found
```bash
/opt/miniconda3/bin/conda init zsh
# Restart terminal
```

### Wrong Python environment
```bash
conda activate rtc-chem
which python  # Should show /opt/miniconda3/envs/rtc-chem/bin/python
```

### Git push rejected
```bash
git pull origin main --rebase
git push origin main
```

### Missing package
```bash
pip install package-name
```

### Verify structure chains
```bash
python -c "from Bio.PDB import PDBParser; parser = PDBParser(QUIET=True); structure = parser.get_structure('check', 'file.pdb'); print([c.id for c in structure.get_chains()])"
```

---

## Quality Control Checklist

Before moving to next step:
- [ ] All expected files created
- [ ] File sizes reasonable (not 0 bytes)
- [ ] Chain composition verified
- [ ] Scripts run without errors
- [ ] Changes committed to Git
- [ ] Documentation updated

---
EOF
