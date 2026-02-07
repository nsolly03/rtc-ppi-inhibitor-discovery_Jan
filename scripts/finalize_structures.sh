#!/bin/bash

cd ~/Desktop/Botanique/Project || exit 1

echo "========================================================================"
echo "FINALIZING STRUCTURE SELECTION"
echo "========================================================================"
echo ""

echo "Moving 7N0C (NSP10-NSP14) to main structures..."
cp data/structures/candidates/7N0C.pdb data/structures/pdb/7N0C.pdb
echo "✓ 7N0C.pdb"

echo ""
echo "Moving 6ZSL (NSP13 Helicase) to main structures..."
cp data/structures/candidates/6ZSL.pdb data/structures/pdb/6ZSL.pdb
echo "✓ 6ZSL.pdb"

echo ""
echo "Cleaning up candidates folder..."
rm -rf data/structures/candidates/
echo "✓ Removed candidates directory"

echo ""
echo "========================================================================"
echo "CURRENT PROJECT STRUCTURES"
echo "========================================================================"
echo ""
ls -lh data/structures/pdb/*.pdb
echo ""

echo "Structure count: $(ls data/structures/pdb/*.pdb | wc -l)"
echo ""
echo "Structure selection complete!"

