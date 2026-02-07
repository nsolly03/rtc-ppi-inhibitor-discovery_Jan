#!/usr/bin/env bash
set -e

mkdir -p data/structures/candidates
mkdir -p docs

echo "Downloading NSP14 candidates..."
curl -s https://files.rcsb.org/download/7N0C.pdb -o data/structures/candidates/7N0C.pdb
curl -s https://files.rcsb.org/download/7DIY.pdb -o data/structures/candidates/7DIY.pdb
curl -s https://files.rcsb.org/download/7MCW.pdb -o data/structures/candidates/7MCW.pdb

echo "Downloading NSP13 candidates..."
curl -s https://files.rcsb.org/download/6ZSL.pdb -o data/structures/candidates/6ZSL.pdb
curl -s https://files.rcsb.org/download/7NIO.pdb -o data/structures/candidates/7NIO.pdb

echo "Now validating all..."

OUT="docs/CANDIDATES_VALIDATION_$(date +%Y%m%d).txt"

for pdb in 7N0C 7DIY 7MCW 6ZSL 7NIO; do
    echo "========================================================================" >> "$OUT"
    python3 scripts/validate_structure.py "$pdb" "data/structures/candidates/${pdb}.pdb" >> "$OUT"
    echo "" >> "$OUT"
done

echo "✓ Validation complete!"
echo "Results saved to: $OUT"

