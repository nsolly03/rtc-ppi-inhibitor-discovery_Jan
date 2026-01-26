#!/bin/bash
# Quick script to add daily log entry
# Usage: bash scripts/add_log_entry.sh

DATE=$(date '+%A, %B %d, %Y')

echo ""
echo "Adding entry for: $DATE"
echo ""

# Open work log
code docs/WORK_LOG.md

echo "✅ Work log opened in VS Code"
echo ""
echo "Add your entry under the date: $DATE"
