#!/bin/bash

echo "========================================================================"
echo "COMPREHENSIVE STRUCTURE VALIDATION REPORT"
echo "Date: $(date)"
echo "========================================================================"
echo ""

for pdb_id in 6W4H 7DFG 6XEZ 7EDI 6W9C; do
    echo ""
    echo "########################################################################"
    echo "# $pdb_id"
    echo "########################################################################"
    echo ""
    
    if [ -f "data/structures/pdb/${pdb_id}.pdb" ]; then
        python3 scripts/validate_structure.py $pdb_id data/structures/pdb/${pdb_id}.pdb
    else
        echo "⚠️  ERROR: File not found: data/structures/pdb/${pdb_id}.pdb"
    fi
    
    echo ""
    echo "------------------------------------------------------------------------"
    echo ""
done

echo "========================================================================"
echo "VALIDATION COMPLETE"
echo "========================================================================"
