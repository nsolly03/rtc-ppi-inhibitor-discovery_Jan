#!/usr/bin/env python3
"""
Add Sections 4-10 to 06_6XEZ_interface_discovery.ipynb
Complete discovery analysis (same as 7DFG but for 6XEZ)
"""

import json

# Read the 7DFG notebook to copy sections 4-10
with open('05_7DFG_interface_discovery.ipynb', 'r') as f:
    notebook_7dfg = json.load(f)

# Read the 6XEZ notebook (sections 1-3 only)
with open('06_6XEZ_interface_discovery.ipynb', 'r') as f:
    notebook_6xez = json.load(f)

print(f"7DFG notebook has {len(notebook_7dfg['cells'])} cells")
print(f"6XEZ notebook has {len(notebook_6xez['cells'])} cells (sections 1-3)")
print()

# Extract sections 4-10 from 7DFG (cells after section 3)
# Assuming section 3 is the 6th cell (markdown+code for sections 1,2,3)
sections_to_copy = notebook_7dfg['cells'][6:]  # Everything after section 3

print(f"Copying {len(sections_to_copy)} cells (sections 4-10) from 7DFG...")

# Replace all "7DFG" with "6XEZ" in the copied sections
for cell in sections_to_copy:
    if cell['cell_type'] == 'code' or cell['cell_type'] == 'markdown':
        if isinstance(cell['source'], list):
            cell['source'] = [line.replace('7DFG', '6XEZ').replace('7dfg', '6xez') for line in cell['source']]
        elif isinstance(cell['source'], str):
            cell['source'] = cell['source'].replace('7DFG', '6XEZ').replace('7dfg', '6xez')

# Add to 6XEZ notebook
notebook_6xez['cells'].extend(sections_to_copy)

# Update the conclusion section to emphasize comparison
for i, cell in enumerate(notebook_6xez['cells']):
    if cell['cell_type'] == 'markdown' and 'Section 10' in cell.get('source', [''])[0]:
        # Found section 10, update it
        notebook_6xez['cells'][i]['source'] = [
            "## Section 10: Conclusions & Comparison with 7DFG\n",
            "\n",
            "### Key Findings\n",
            "\n",
            "This analysis discovered hot spots in **6XEZ** using the same systematic approach as 7DFG.\n",
            "\n",
            "**Comparison Questions:**\n",
            "1. Are the hot spots conserved between 6XEZ and 7DFG?\n",
            "2. Are the hot spot residues identical or shifted?\n",
            "3. Which structure has better-defined interfaces?\n",
            "4. Should we use one structure or both for virtual screening?\n",
            "\n",
            "**Next Steps:**\n",
            "1. **Load 7DFG results** for side-by-side comparison\n",
            "2. **Create comparison notebook** (07_compare_7DFG_6XEZ.ipynb)\n",
            "3. **Decide strategy:** Use best structure OR use both\n",
            "4. **Run fpocket** on selected structure(s)\n",
            "5. **Proceed to 7EDI** (NSP10-NSP14 interface)\n",
            "\n",
            "### Analysis Complete! ✅\n",
            "\n",
            "**Notebook:** `06_6XEZ_interface_discovery.ipynb`  \n",
            "**Results:** `data/analysis_results/6XEZ_*.csv|json`  \n",
            "**Status:** Ready for comparison with 7DFG"
        ]

# Save updated notebook
with open('06_6XEZ_interface_discovery.ipynb', 'w') as f:
    json.dump(notebook_6xez, f, indent=1)

print(f"✓ Added sections 4-10 to 6XEZ notebook")
print(f"  Total cells now: {len(notebook_6xez['cells'])}")
print()
print("=" * 70)
print("✅ 6XEZ NOTEBOOK COMPLETE!")
print("=" * 70)
print()
print("You can now:")
print("  1. Run 06_6XEZ_interface_discovery.ipynb")
print("  2. Compare results with 7DFG")
print("  3. Create comparison notebook (07_compare_7DFG_6XEZ.ipynb)")
