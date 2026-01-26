# Quick Reference Guide

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
EOF
