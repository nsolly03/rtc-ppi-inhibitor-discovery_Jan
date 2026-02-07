#!/bin/bash

echo "========================================================================"
echo "DOWNLOADING AND VALIDATING CANDIDATE STRUCTURES"
echo "========================================================================"
echo ""

# Candidate structures to check
NSP14_CANDIDATES="7N0C 7DIY 7MCW"
NSP13_CANDIDATES="6ZSL 7NIO"

TEMP_DIR="data/structures/candidates"
mkdir -p "$TEMP_DIR"

echo "📥 PHASE 1: DOWNLOADING CANDIDATE NSP14 STRUCTURES"
echo "========================================================================"
for pdb in $NSP14_CANDIDATES; do
    echo ""
    echo "Downloading ${pdb}..."
    curl -s "https://files.rcsb.org/download/${pdb}.pdb" -o "$TEMP_DIR/${pdb}.pdb"
    
    if [ -f "$TEMP_DIR/${pdb}.pdb" ]; then
        echo "✓ Downloaded ${pdb}.pdb"
        
        # Validate immediately
        echo ""
        echo "Validating ${pdb}..."
        python3 scripts/validate_structure.py $pdb "$TEMP_DIR/${pdb}.pdb"
        
        echo ""
        echo "------------------------------------------------------------------------"
        echo "Press Enter to continue to next structure..."
        read
    else
        echo "✗ Failed to download ${pdb}.pdb"
    fi
done

echo ""
echo "========================================================================"
echo "📥 PHASE 2: DOWNLOADING CANDIDATE NSP13 STRUCTURES"
echo "========================================================================"
for pdb in $NSP13_CANDIDATES; do
    echo ""
    echo "Downloading ${pdb}..."
    curl -s "https://files.rcsb.org/download/${pdb}.pdb" -o "$TEMP_DIR/${pdb}.pdb"
    
    if [ -f "$TEMP_DIR/${pdb}.pdb" ]; then
        echo "✓ Downloaded ${pdb}.pdb"
        
        # Validate immediately
        echo ""
        echo "Validating ${pdb}..."
        python3 scripts/validate_structure.py $pdb "$TEMP_DIR/${pdb}.pdb"
        
        echo ""
        echo "------------------------------------------------------------------------"
        echo "Press Enter to continue to next structure..."
        read
    else
        echo "✗ Failed to download ${pdb}.pdb"
    fi
done

echo ""
echo "========================================================================"
echo "DOWNLOAD AND VALIDATION COMPLETE"
echo "========================================================================"
echo ""
echo "All candidate structures downloaded to: $TEMP_DIR"
echo ""
echo "Next steps:"
echo "  1. Review validation results above"
echo "  2. Choose best NSP14 structure (7N0C, 7DIY, or 7MCW)"
echo "  3. Choose best NSP13 structure (6ZSL or 7NIO)"
echo "  4. Move chosen structures to data/structures/pdb/"
echo "  5. Delete candidates directory"
