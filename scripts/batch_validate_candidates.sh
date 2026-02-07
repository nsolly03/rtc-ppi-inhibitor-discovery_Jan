#!/bin/bash

echo "========================================================================"
echo "BATCH DOWNLOAD AND VALIDATION - CANDIDATE STRUCTURES"
echo "Date: $(date)"
echo "========================================================================"
echo ""

CANDIDATES_DIR="data/structures/candidates"
mkdir -p "$CANDIDATES_DIR"

NSP14_CANDIDATES="7N0C 7DIY 7MCW"
NSP13_CANDIDATES="6ZSL 7NIO"

echo "DOWNLOADING NSP14 CANDIDATES..."
for pdb in $NSP14_CANDIDATES; do
    curl -s "https://files.rcsb.org/download/${pdb}.pdb" -o "$CANDIDATES_DIR/${pdb}.pdb"
    echo "  ✓ ${pdb}"
done

echo ""
echo "DOWNLOADING NSP13 CANDIDATES..."
for pdb in $NSP13_CANDIDATES; do
    curl -s "https://files.rcsb.org/download/${pdb}.pdb" -o "$CANDIDATES_DIR/${pdb}.pdb"
    echo "  ✓ ${pdb}"
done

echo ""
echo "VALIDATING ALL CANDIDATES..."

ALL_CANDIDATES="$NSP14_CANDIDATES $NSP13_CANDIDATES"

for pdb in $ALL_CANDIDATES; do
    echo ""
    echo "VALIDATING $pdb"
    python3 scripts/validate_structure.py "$pdb" "$CANDIDATES_DIR/${pdb}.pdb"
done

echo ""
echo "DONE."
