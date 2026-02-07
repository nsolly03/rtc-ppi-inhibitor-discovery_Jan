#!/bin/bash

echo "========================================================================"
echo "FINAL CANDIDATE REVIEW - FULL RCSB VALIDATION"
echo "========================================================================"
echo ""

echo "########################################################################"
echo "# 7N0C - NSP10-NSP14 (Our NSP14 choice)"
echo "########################################################################"
python3 scripts/validate_structure.py 7N0C data/structures/candidates/7N0C.pdb

echo ""
echo ""
echo "########################################################################"
echo "# 6ZSL - NSP13 Helicase (Candidate 1)"
echo "########################################################################"
python3 scripts/validate_structure.py 6ZSL data/structures/candidates/6ZSL.pdb

echo ""
echo ""
echo "########################################################################"
echo "# 7NIO - NSP13 Helicase (Candidate 2)"
echo "########################################################################"
python3 scripts/validate_structure.py 7NIO data/structures/candidates/7NIO.pdb

echo ""
echo "========================================================================"
echo "REVIEW COMPLETE"
echo "========================================================================"
