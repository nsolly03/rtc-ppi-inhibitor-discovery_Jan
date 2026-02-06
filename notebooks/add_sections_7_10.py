#!/usr/bin/env python3
"""
Add Sections 7-10 to 05_7DFG_interface_discovery.ipynb
FINAL PART: Grid Boxes, Summary, Export
"""

import json

# Read existing notebook
with open('05_7DFG_interface_discovery.ipynb', 'r') as f:
    notebook = json.load(f)

print(f"Current notebook has {len(notebook['cells'])} cells")

# Sections 7-10 to add
new_cells = [
    # Section 7: Define Docking Grid Boxes
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 7: Define Docking Grid Boxes"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "print(\"=\"*70)\n",
            "print(\"DOCKING GRID BOX PARAMETERS\")\n",
            "print(\"=\"*70)\n",
            "print()\n",
            "\n",
            "grid_boxes = {}\n",
            "\n",
            "# NSP12-NSP7 Grid Box\n",
            "if not df_hot_12_7.empty:\n",
            "    top_spot = df_hot_12_7.iloc[0]\n",
            "    res_12 = model[chain_assignments['nsp12']][top_spot['NSP12_PDB']]\n",
            "    center_12_7 = res_12['CA'].get_coord()\n",
            "    \n",
            "    grid_boxes['NSP12-NSP7'] = {\n",
            "        'center_x': float(center_12_7[0]),\n",
            "        'center_y': float(center_12_7[1]),\n",
            "        'center_z': float(center_12_7[2]),\n",
            "        'size_x': 25.0,\n",
            "        'size_y': 25.0,\n",
            "        'size_z': 25.0,\n",
            "        'hot_spot': f\"{top_spot['NSP12_Res']}{top_spot['NSP12_PDB']}-{top_spot['NSP7_Res']}{top_spot['NSP7_PDB']}\"\n",
            "    }\n",
            "    \n",
            "    print(\"INTERFACE 1: NSP12-NSP7\")\n",
            "    print(f\"  Target: {grid_boxes['NSP12-NSP7']['hot_spot']}\")\n",
            "    print(f\"  Center: ({center_12_7[0]:.3f}, {center_12_7[1]:.3f}, {center_12_7[2]:.3f})\")\n",
            "    print(f\"  Size: 25 × 25 × 25 Ų\")\n",
            "    print()\n",
            "\n",
            "# NSP12-NSP8 Grid Box\n",
            "if not df_hot_12_8.empty:\n",
            "    top_spot = df_hot_12_8.iloc[0]\n",
            "    res_12 = model[chain_assignments['nsp12']][top_spot['NSP12_PDB']]\n",
            "    center_12_8 = res_12['CA'].get_coord()\n",
            "    \n",
            "    grid_boxes['NSP12-NSP8'] = {\n",
            "        'center_x': float(center_12_8[0]),\n",
            "        'center_y': float(center_12_8[1]),\n",
            "        'center_z': float(center_12_8[2]),\n",
            "        'size_x': 25.0,\n",
            "        'size_y': 25.0,\n",
            "        'size_z': 25.0,\n",
            "        'hot_spot': f\"{top_spot['NSP12_Res']}{top_spot['NSP12_PDB']}-{top_spot['NSP8_Res']}{top_spot['NSP8_PDB']}\"\n",
            "    }\n",
            "    \n",
            "    print(\"INTERFACE 2: NSP12-NSP8\")\n",
            "    print(f\"  Target: {grid_boxes['NSP12-NSP8']['hot_spot']}\")\n",
            "    print(f\"  Center: ({center_12_8[0]:.3f}, {center_12_8[1]:.3f}, {center_12_8[2]:.3f})\")\n",
            "    print(f\"  Size: 25 × 25 × 25 Ų\")\n",
            "    print()\n",
            "\n",
            "print(\"=\"*70)\n",
            "print(\"FOR AUTODOCK VINA CONFIG FILES:\")\n",
            "print(\"=\"*70)\n",
            "print()\n",
            "\n",
            "for interface, params in grid_boxes.items():\n",
            "    print(f\"# {interface}\")\n",
            "    print(f\"center_x = {params['center_x']:.3f}\")\n",
            "    print(f\"center_y = {params['center_y']:.3f}\")\n",
            "    print(f\"center_z = {params['center_z']:.3f}\")\n",
            "    print(f\"size_x = {params['size_x']}\")\n",
            "    print(f\"size_y = {params['size_y']}\")\n",
            "    print(f\"size_z = {params['size_z']}\")\n",
            "    print()\n",
            "\n",
            "print(\"=\"*70)"
        ]
    },
    # Section 8: Summary Table
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 8: Summary Table"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Build comprehensive summary\n",
            "summary_data = []\n",
            "\n",
            "summary_data.append(['Property', 'Value'])\n",
            "summary_data.append(['PDB ID', '7DFG'])\n",
            "summary_data.append(['Analysis Date', datetime.now().strftime('%Y-%m-%d')])\n",
            "summary_data.append(['Approach', 'Discovery-based (no prior knowledge)'])\n",
            "summary_data.append(['', ''])\n",
            "\n",
            "# Chain assignments\n",
            "summary_data.append(['CHAIN ASSIGNMENTS', ''])\n",
            "if 'nsp12' in chain_assignments:\n",
            "    summary_data.append(['NSP12', f\"Chain {chain_assignments['nsp12']}'])\n",
            "if 'nsp7' in chain_assignments:\n",
            "    summary_data.append(['NSP7', f\"Chain {chain_assignments['nsp7']}'])\n",
            "if 'nsp8' in chain_assignments:\n",
            "    summary_data.append(['NSP8', f\"Chain {chain_assignments['nsp8']}'])\n",
            "summary_data.append(['', ''])\n",
            "\n",
            "# NSP12-NSP7 Interface\n",
            "summary_data.append(['INTERFACE 1: NSP12-NSP7', ''])\n",
            "if not df_hot_12_7.empty:\n",
            "    top = df_hot_12_7.iloc[0]\n",
            "    summary_data.append(['Hot Spots Found', str(len(df_hot_12_7))])\n",
            "    summary_data.append(['Strongest Hot Spot', f\"{top['NSP12_Res']}{top['NSP12_PDB']}-{top['NSP7_Res']}{top['NSP7_PDB']}'])\n",
            "    summary_data.append(['Distance', f\"{top['Distance']:.2f} Å\"])\n",
            "    summary_data.append(['Type', top['Type']])\n",
            "    if not df_interface_12_7.empty:\n",
            "        summary_data.append(['Interface Residues', f\"{len(df_interface_12_7)} ({nsp12_count} NSP12, {nsp7_count} NSP7)\"])\n",
            "    if 'NSP12-NSP7' in grid_boxes:\n",
            "        gb = grid_boxes['NSP12-NSP7']\n",
            "        summary_data.append(['Grid Box Center', f\"({gb['center_x']:.1f}, {gb['center_y']:.1f}, {gb['center_z']:.1f})\"])\n",
            "        summary_data.append(['Grid Box Size', '25 × 25 × 25 Ų'])\n",
            "else:\n",
            "    summary_data.append(['Status', 'No hot spots found'])\n",
            "summary_data.append(['', ''])\n",
            "\n",
            "# NSP12-NSP8 Interface\n",
            "summary_data.append(['INTERFACE 2: NSP12-NSP8', ''])\n",
            "if not df_hot_12_8.empty:\n",
            "    top = df_hot_12_8.iloc[0]\n",
            "    summary_data.append(['Hot Spots Found', str(len(df_hot_12_8))])\n",
            "    summary_data.append(['Strongest Hot Spot', f\"{top['NSP12_Res']}{top['NSP12_PDB']}-{top['NSP8_Res']}{top['NSP8_PDB']}'])\n",
            "    summary_data.append(['Distance', f\"{top['Distance']:.2f} Å\"])\n",
            "    summary_data.append(['Type', top['Type']])\n",
            "    if not df_interface_12_8.empty:\n",
            "        summary_data.append(['Interface Residues', f\"{len(df_interface_12_8)} ({nsp12_count_8} NSP12, {nsp8_count} NSP8)\"])\n",
            "    if 'NSP12-NSP8' in grid_boxes:\n",
            "        gb = grid_boxes['NSP12-NSP8']\n",
            "        summary_data.append(['Grid Box Center', f\"({gb['center_x']:.1f}, {gb['center_y']:.1f}, {gb['center_z']:.1f})\"])\n",
            "        summary_data.append(['Grid Box Size', '25 × 25 × 25 Ų'])\n",
            "else:\n",
            "    summary_data.append(['Status', 'No hot spots found'])\n",
            "\n",
            "# Display summary\n",
            "df_summary = pd.DataFrame(summary_data[1:], columns=summary_data[0])\n",
            "\n",
            "print(\"=\"*70)\n",
            "print(\"ANALYSIS SUMMARY - 7DFG\")\n",
            "print(\"=\"*70)\n",
            "print()\n",
            "print(df_summary.to_string(index=False))\n",
            "print()\n",
            "print(\"=\"*70)"
        ]
    },
    # Section 9: Export Results
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 9: Export Results"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "print(\"Exporting results...\")\n",
            "print()\n",
            "\n",
            "# 1. Export NSP12-NSP7 hot spots\n",
            "if not df_hot_12_7.empty:\n",
            "    csv_file = RESULTS_DIR / '7DFG_NSP12-NSP7_hotspots.csv'\n",
            "    df_hot_12_7.to_csv(csv_file, index=False)\n",
            "    print(f\"✓ NSP12-NSP7 hot spots: {csv_file}\")\n",
            "\n",
            "# 2. Export NSP12-NSP8 hot spots\n",
            "if not df_hot_12_8.empty:\n",
            "    csv_file = RESULTS_DIR / '7DFG_NSP12-NSP8_hotspots.csv'\n",
            "    df_hot_12_8.to_csv(csv_file, index=False)\n",
            "    print(f\"✓ NSP12-NSP8 hot spots: {csv_file}\")\n",
            "\n",
            "# 3. Export NSP12-NSP7 interface\n",
            "if not df_interface_12_7.empty:\n",
            "    csv_file = RESULTS_DIR / '7DFG_NSP12-NSP7_interface.csv'\n",
            "    df_interface_12_7.to_csv(csv_file, index=False)\n",
            "    print(f\"✓ NSP12-NSP7 interface: {csv_file}\")\n",
            "\n",
            "# 4. Export NSP12-NSP8 interface\n",
            "if not df_interface_12_8.empty:\n",
            "    csv_file = RESULTS_DIR / '7DFG_NSP12-NSP8_interface.csv'\n",
            "    df_interface_12_8.to_csv(csv_file, index=False)\n",
            "    print(f\"✓ NSP12-NSP8 interface: {csv_file}\")\n",
            "\n",
            "# 5. Export analysis JSON\n",
            "analysis_data = {\n",
            "    'pdb_id': '7DFG',\n",
            "    'analysis_date': datetime.now().strftime('%Y-%m-%d'),\n",
            "    'approach': 'discovery-based',\n",
            "    'chains': chain_assignments,\n",
            "    'interfaces': {}\n",
            "}\n",
            "\n",
            "if not df_hot_12_7.empty:\n",
            "    top = df_hot_12_7.iloc[0]\n",
            "    analysis_data['interfaces']['NSP12-NSP7'] = {\n",
            "        'hot_spots_count': len(df_hot_12_7),\n",
            "        'strongest_hot_spot': {\n",
            "            'nsp12_residue': f\"{top['NSP12_Res']}{top['NSP12_PDB']}\",\n",
            "            'nsp7_residue': f\"{top['NSP7_Res']}{top['NSP7_PDB']}\",\n",
            "            'distance': float(top['Distance']),\n",
            "            'type': top['Type']\n",
            "        },\n",
            "        'interface_residues': len(df_interface_12_7) if not df_interface_12_7.empty else 0,\n",
            "        'grid_box': grid_boxes.get('NSP12-NSP7', {})\n",
            "    }\n",
            "\n",
            "if not df_hot_12_8.empty:\n",
            "    top = df_hot_12_8.iloc[0]\n",
            "    analysis_data['interfaces']['NSP12-NSP8'] = {\n",
            "        'hot_spots_count': len(df_hot_12_8),\n",
            "        'strongest_hot_spot': {\n",
            "            'nsp12_residue': f\"{top['NSP12_Res']}{top['NSP12_PDB']}\",\n",
            "            'nsp8_residue': f\"{top['NSP8_Res']}{top['NSP8_PDB']}\",\n",
            "            'distance': float(top['Distance']),\n",
            "            'type': top['Type']\n",
            "        },\n",
            "        'interface_residues': len(df_interface_12_8) if not df_interface_12_8.empty else 0,\n",
            "        'grid_box': grid_boxes.get('NSP12-NSP8', {})\n",
            "    }\n",
            "\n",
            "json_file = RESULTS_DIR / '7DFG_analysis.json'\n",
            "with open(json_file, 'w') as f:\n",
            "    json.dump(analysis_data, f, indent=2)\n",
            "print(f\"✓ Analysis data: {json_file}\")\n",
            "\n",
            "# 6. Export summary\n",
            "summary_file = RESULTS_DIR / '7DFG_summary.csv'\n",
            "df_summary.to_csv(summary_file, index=False)\n",
            "print(f\"✓ Summary table: {summary_file}\")\n",
            "\n",
            "print()\n",
            "print(\"✓ All results exported successfully!\")\n",
            "print()\n",
            "print(f\"Results directory: {RESULTS_DIR.resolve()}\")"
        ]
    },
    # Section 10: Conclusions
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 10: Conclusions\n",
            "\n",
            "### Key Findings\n",
            "\n",
            "This analysis used a **discovery-based approach** to identify druggable hot spots in the 7DFG RdRp holoenzyme without prior knowledge.\n",
            "\n",
            "**Methodology:**\n",
            "- Systematic screening of all charged residue pairs (Lys-Asp, Arg-Asp, Lys-Glu, Arg-Glu)\n",
            "- Distance cutoff: < 5 Å for hot spots, < 10 Å for interface mapping\n",
            "- Grid boxes: 25 × 25 × 25 Ų centered on strongest hot spots\n",
            "\n",
            "**Next Steps:**\n",
            "1. **Week 3-4:** Run fpocket on both interfaces to identify binding pockets\n",
            "2. **Week 4:** Analyze 6XEZ (alternative RdRp structure) and compare with 7DFG\n",
            "3. **Week 4:** Analyze 7EDI (NSP10-NSP14 ExoN complex)\n",
            "4. **Week 5:** Comparative analysis and target prioritization\n",
            "5. **Month 2+:** Virtual screening campaign\n",
            "\n",
            "### Analysis Complete! ✅\n",
            "\n",
            "**Notebook:** `05_7DFG_interface_discovery.ipynb`  \n",
            "**Results:** `data/analysis_results/7DFG_*.csv|json`  \n",
            "**Status:** Ready for fpocket analysis"
        ]
    }
]

print("Adding Sections 7, 8, 9, 10...")

# Add cells to notebook
notebook['cells'].extend(new_cells)

# Save
with open('05_7DFG_interface_discovery.ipynb', 'w') as f:
    json.dump(notebook, f, indent=1)

print(f"✓ Added Sections 7, 8, 9, 10")
print(f"  Total cells now: {len(notebook['cells'])}")
print()
print("=" * 70)
print("✅ NOTEBOOK COMPLETE!")
print("=" * 70)
print()
print("All sections (1-10) have been added to:")
print("  05_7DFG_interface_discovery.ipynb")
print()
print("You can now:")
print("  1. Open the notebook in VS Code or Jupyter")
print("  2. Run all cells to perform the complete 7DFG analysis")
print("  3. Review the discovered hot spots")
print("  4. Proceed with fpocket analysis")
