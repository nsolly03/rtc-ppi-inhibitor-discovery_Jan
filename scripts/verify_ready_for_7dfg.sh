#!/bin/bash

echo "========================================================================"
echo "PRE-RUN VERIFICATION: 7DFG Analysis"
echo "========================================================================"
echo ""

# Check file
if [ -f "data/structures/pdb/7DFG.pdb" ]; then
    echo "✓ PDB file exists: 7DFG.pdb"
    ls -lh data/structures/pdb/7DFG.pdb
else
    echo "✗ ERROR: 7DFG.pdb not found!"
    exit 1
fi

echo ""

# Check notebook
if [ -f "notebooks/05_7DFG_interface_discovery.ipynb" ]; then
    echo "✓ Notebook exists: 05_7DFG_interface_discovery.ipynb"
    cell_count=$(grep -o '"cell_type"' notebooks/05_7DFG_interface_discovery.ipynb | wc -l)
    echo "  Cells: $cell_count"
else
    echo "✗ ERROR: Notebook not found!"
    exit 1
fi

echo ""

# Check results directory
if [ -d "data/analysis_results" ]; then
    echo "✓ Results directory exists"
else
    echo "⚠️  Creating results directory..."
    mkdir -p data/analysis_results
fi

echo ""

# Check conda environment
if conda env list | grep -q "rtc-chem"; then
    echo "✓ Conda environment 'rtc-chem' exists"
else
    echo "✗ ERROR: rtc-chem environment not found!"
    exit 1
fi

echo ""

# Check Python packages
echo "Checking Python packages..."
conda activate rtc-chem 2>/dev/null || source activate rtc-chem 2>/dev/null

python3 << 'PYCHECK'
try:
    import py3Dmol
    print("  ✓ py3Dmol")
except:
    print("  ✗ py3Dmol (install: pip install py3Dmol --break-system-packages)")

try:
    from Bio import PDB
    print("  ✓ BioPython")
except:
    print("  ✗ BioPython (install: pip install biopython --break-system-packages)")

try:
    import pandas
    print("  ✓ pandas")
except:
    print("  ✗ pandas")

try:
    import numpy
    print("  ✓ numpy")
except:
    print("  ✗ numpy")
PYCHECK

echo ""
echo "========================================================================"
echo "✅ VERIFICATION COMPLETE"
echo "========================================================================"
echo ""
echo "You're ready to run: notebooks/05_7DFG_interface_discovery.ipynb"
echo ""
echo "To open in VS Code:"
echo "  code notebooks/05_7DFG_interface_discovery.ipynb"
echo ""
echo "To run in Jupyter:"
echo "  jupyter notebook notebooks/05_7DFG_interface_discovery.ipynb"
echo ""
echo "Don't forget to select the 'rtc-chem' kernel!"
