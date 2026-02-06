#!/usr/bin/env python3
"""
Add Sections 5-10 to 05_7DFG_interface_discovery.ipynb
Part 2: Visualization, Interface Mapping, Grid Boxes, Summary, Export
"""

import json

# Read existing notebook
with open('05_7DFG_interface_discovery.ipynb', 'r') as f:
    notebook = json.load(f)

print(f"Current notebook has {len(notebook['cells'])} cells")

# Sections 5-10 to add
new_cells = [
    # Section 5A Header
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 5A: Visualize NSP12-NSP7 Hot Spots"
        ]
    },
    # Section 5A Code
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "if not df_hot_12_7.empty and len(df_hot_12_7) > 0:\\n",
            "    # Get top hot spot\\n",
            "    top_spot = df_hot_12_7.iloc[0]\\n",
            "    \\n",
            "    print(f\\\"Visualizing NSP12-NSP7 strongest hot spot...\\\")\\n",
            "    print(f\\\"Hot spot: {top_spot['NSP12_Res']}{top_spot['NSP12_PDB']} - {top_spot['NSP7_Res']}{top_spot['NSP7_PDB']}\\\")\\n",
            "    print(f\\\"Distance: {top_spot['Distance']:.2f} Å\\\")\\n",
            "    print()\\n",
            "    \\n",
            "    # Create viewer\\n",
            "    with open(PDB_FILE, 'r') as f:\\n",
            "        pdb_data = f.read()\\n",
            "    \\n",
            "    view2 = py3Dmol.view(width=800, height=600)\\n",
            "    view2.addModel(pdb_data, 'pdb')\\n",
            "    \\n",
            "    # Proteins (semi-transparent)\\n",
            "    view2.setStyle({'chain': chain_assignments['nsp12']}, {'cartoon': {'color': 'cyan', 'opacity': 0.5}})\\n",
            "    view2.setStyle({'chain': chain_assignments['nsp7']}, {'cartoon': {'color': 'magenta', 'opacity': 0.5}})\\n",
            "    if 'nsp8' in chain_assignments:\\n",
            "        view2.setStyle({'chain': chain_assignments['nsp8']}, {'cartoon': {'color': 'yellow', 'opacity': 0.3}})\\n",
            "    \\n",
            "    # Hot spots (large spheres)\\n",
            "    view2.addStyle({'chain': chain_assignments['nsp12'], 'resi': top_spot['NSP12_PDB']},\\n",
            "                   {'sphere': {'color': 'red', 'radius': 2.0}})\\n",
            "    view2.addStyle({'chain': chain_assignments['nsp7'], 'resi': top_spot['NSP7_PDB']},\\n",
            "                   {'sphere': {'color': 'blue', 'radius': 2.0}})\\n",
            "    \\n",
            "    # Labels\\n",
            "    view2.addLabel(f\\\"{top_spot['NSP12_Res']}{top_spot['NSP12_PDB']}\\\\n{top_spot['Distance']:.1f} Å\\\",\\n",
            "                   {'fontColor': 'red', 'fontSize': 14, 'backgroundColor': 'white', 'backgroundOpacity': 0.8},\\n",
            "                   {'chain': chain_assignments['nsp12'], 'resi': top_spot['NSP12_PDB']})\\n",
            "    view2.addLabel(f\\\"{top_spot['NSP7_Res']}{top_spot['NSP7_PDB']}\\\",\\n",
            "                   {'fontColor': 'blue', 'fontSize': 14, 'backgroundColor': 'white', 'backgroundOpacity': 0.8},\\n",
            "                   {'chain': chain_assignments['nsp7'], 'resi': top_spot['NSP7_PDB']})\\n",
            "    \\n",
            "    view2.zoomTo({'chain': chain_assignments['nsp12'], 'resi': top_spot['NSP12_PDB']})\\n",
            "    view2.setBackgroundColor('white')\\n",
            "    \\n",
            "    print(\\\"Visualization:\\\")\\n",
            "    print(f\\\"  🔴 RED sphere = NSP12 {top_spot['NSP12_Res']} {top_spot['NSP12_PDB']}\\\")\\n",
            "    print(f\\\"  🔵 BLUE sphere = NSP7 {top_spot['NSP7_Res']} {top_spot['NSP7_PDB']}\\\")\\n",
            "    print(f\\\"  Distance: {top_spot['Distance']:.2f} Å\\\")\\n",
            "    \\n",
            "    view2.show()\\n",
            "else:\\n",
            "    print(\\\"⚠️  No hot spots to visualize for NSP12-NSP7 interface\\\")"
        ]
    },
    # Section 5B Header
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 5B: Visualize NSP12-NSP8 Hot Spots"
        ]
    },
    # Section 5B Code
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "if not df_hot_12_8.empty and len(df_hot_12_8) > 0:\\n",
            "    # Get top hot spot\\n",
            "    top_spot = df_hot_12_8.iloc[0]\\n",
            "    \\n",
            "    print(f\\\"Visualizing NSP12-NSP8 strongest hot spot...\\\")\\n",
            "    print(f\\\"Hot spot: {top_spot['NSP12_Res']}{top_spot['NSP12_PDB']} - {top_spot['NSP8_Res']}{top_spot['NSP8_PDB']}\\\")\\n",
            "    print(f\\\"Distance: {top_spot['Distance']:.2f} Å\\\")\\n",
            "    print()\\n",
            "    \\n",
            "    # Create viewer\\n",
            "    with open(PDB_FILE, 'r') as f:\\n",
            "        pdb_data = f.read()\\n",
            "    \\n",
            "    view3 = py3Dmol.view(width=800, height=600)\\n",
            "    view3.addModel(pdb_data, 'pdb')\\n",
            "    \\n",
            "    # Proteins (semi-transparent)\\n",
            "    view3.setStyle({'chain': chain_assignments['nsp12']}, {'cartoon': {'color': 'cyan', 'opacity': 0.5}})\\n",
            "    view3.setStyle({'chain': chain_assignments['nsp8']}, {'cartoon': {'color': 'yellow', 'opacity': 0.5}})\\n",
            "    if 'nsp7' in chain_assignments:\\n",
            "        view3.setStyle({'chain': chain_assignments['nsp7']}, {'cartoon': {'color': 'magenta', 'opacity': 0.3}})\\n",
            "    \\n",
            "    # Hot spots (large spheres)\\n",
            "    view3.addStyle({'chain': chain_assignments['nsp12'], 'resi': top_spot['NSP12_PDB']},\\n",
            "                   {'sphere': {'color': 'red', 'radius': 2.0}})\\n",
            "    view3.addStyle({'chain': chain_assignments['nsp8'], 'resi': top_spot['NSP8_PDB']},\\n",
            "                   {'sphere': {'color': 'blue', 'radius': 2.0}})\\n",
            "    \\n",
            "    # Labels\\n",
            "    view3.addLabel(f\\\"{top_spot['NSP12_Res']}{top_spot['NSP12_PDB']}\\\\n{top_spot['Distance']:.1f} Å\\\",\\n",
            "                   {'fontColor': 'red', 'fontSize': 14, 'backgroundColor': 'white', 'backgroundOpacity': 0.8},\\n",
            "                   {'chain': chain_assignments['nsp12'], 'resi': top_spot['NSP12_PDB']})\\n",
            "    view3.addLabel(f\\\"{top_spot['NSP8_Res']}{top_spot['NSP8_PDB']}\\\",\\n",
            "                   {'fontColor': 'blue', 'fontSize': 14, 'backgroundColor': 'white', 'backgroundOpacity': 0.8},\\n",
            "                   {'chain': chain_assignments['nsp8'], 'resi': top_spot['NSP8_PDB']})\\n",
            "    \\n",
            "    view3.zoomTo({'chain': chain_assignments['nsp12'], 'resi': top_spot['NSP12_PDB']})\\n",
            "    view3.setBackgroundColor('white')\\n",
            "    \\n",
            "    print(\\\"Visualization:\\\")\\n",
            "    print(f\\\"  🔴 RED sphere = NSP12 {top_spot['NSP12_Res']} {top_spot['NSP12_PDB']}\\\")\\n",
            "    print(f\\\"  🔵 BLUE sphere = NSP8 {top_spot['NSP8_Res']} {top_spot['NSP8_PDB']}\\\")\\n",
            "    print(f\\\"  Distance: {top_spot['Distance']:.2f} Å\\\")\\n",
            "    \\n",
            "    view3.show()\\n",
            "else:\\n",
            "    print(\\\"⚠️  No hot spots to visualize for NSP12-NSP8 interface\\\")"
        ]
    },
    # Section 6A Header
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 6A: Map NSP12-NSP7 Complete Interface (10 Å)"
        ]
    },
    # Section 6A Code
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "print(\\\"=\\\"*70)\\n",
            "print(\\\"MAPPING NSP12-NSP7 INTERFACE (10 Å cutoff)\\\")\\n",
            "print(\\\"=\\\"*70)\\n",
            "print()\\n",
            "\\n",
            "if not df_hot_12_7.empty and 'nsp12' in chain_assignments and 'nsp7' in chain_assignments:\\n",
            "    # Get hot spot center\\n",
            "    top_spot = df_hot_12_7.iloc[0]\\n",
            "    res_12 = model[chain_assignments['nsp12']][top_spot['NSP12_PDB']]\\n",
            "    center = res_12['CA'].get_coord()\\n",
            "    \\n",
            "    print(f\\\"Center point: NSP12 {top_spot['NSP12_Res']}{top_spot['NSP12_PDB']}\\\")\\n",
            "    print(f\\\"Coordinates: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})\\\")\\n",
            "    print()\\n",
            "    print(\\\"Finding all residues within 10 Å...\\\")\\n",
            "    print()\\n",
            "    \\n",
            "    interface_12_7_full = []\\n",
            "    for chain in [model[chain_assignments['nsp12']], model[chain_assignments['nsp7']]]:\\n",
            "        for res in chain:\\n",
            "            if res.id[0] == ' ' and 'CA' in res:\\n",
            "                ca = res['CA'].get_coord()\\n",
            "                dist = np.linalg.norm(ca - center)\\n",
            "                \\n",
            "                if dist <= 10.0:\\n",
            "                    protein = 'NSP12' if chain.id == chain_assignments['nsp12'] else 'NSP7'\\n",
            "                    interface_12_7_full.append({\\n",
            "                        'Protein': protein,\\n",
            "                        'Chain': chain.id,\\n",
            "                        'Residue': res.get_resname(),\\n",
            "                        'PDB_Num': res.id[1],\\n",
            "                        'Distance': dist\\n",
            "                    })\\n",
            "    \\n",
            "    df_interface_12_7 = pd.DataFrame(interface_12_7_full).sort_values('Distance').reset_index(drop=True)\\n",
            "    \\n",
            "    print(f\\\"✓ Found {len(df_interface_12_7)} interface residues:\\\")\\n",
            "    print()\\n",
            "    print(df_interface_12_7.to_string(index=False))\\n",
            "    print()\\n",
            "    \\n",
            "    nsp12_count = len(df_interface_12_7[df_interface_12_7['Protein'] == 'NSP12'])\\n",
            "    nsp7_count = len(df_interface_12_7[df_interface_12_7['Protein'] == 'NSP7'])\\n",
            "    \\n",
            "    print(\\\"Interface composition:\\\")\\n",
            "    print(f\\\"  NSP12: {nsp12_count} residues\\\")\\n",
            "    print(f\\\"  NSP7: {nsp7_count} residues\\\")\\n",
            "    print(f\\\"  Total: {len(df_interface_12_7)} residues\\\")\\n",
            "else:\\n",
            "    print(\\\"⚠️  Cannot map interface - missing hot spots or chains\\\")\\n",
            "    df_interface_12_7 = pd.DataFrame()\\n",
            "    nsp12_count = nsp7_count = 0\\n",
            "\\n",
            "print()\\n",
            "print(\\\"=\\\"*70)"
        ]
    },
    # Section 6B Header
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## Section 6B: Map NSP12-NSP8 Complete Interface (10 Å)"
        ]
    },
    # Section 6B Code
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "print(\\\"=\\\"*70)\\n",
            "print(\\\"MAPPING NSP12-NSP8 INTERFACE (10 Å cutoff)\\\")\\n",
            "print(\\\"=\\\"*70)\\n",
            "print()\\n",
            "\\n",
            "if not df_hot_12_8.empty and 'nsp12' in chain_assignments and 'nsp8' in chain_assignments:\\n",
            "    # Get hot spot center\\n",
            "    top_spot = df_hot_12_8.iloc[0]\\n",
            "    res_12 = model[chain_assignments['nsp12']][top_spot['NSP12_PDB']]\\n",
            "    center = res_12['CA'].get_coord()\\n",
            "    \\n",
            "    print(f\\\"Center point: NSP12 {top_spot['NSP12_Res']}{top_spot['NSP12_PDB']}\\\")\\n",
            "    print(f\\\"Coordinates: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})\\\")\\n",
            "    print()\\n",
            "    print(\\\"Finding all residues within 10 Å...\\\")\\n",
            "    print()\\n",
            "    \\n",
            "    interface_12_8_full = []\\n",
            "    for chain in [model[chain_assignments['nsp12']], model[chain_assignments['nsp8']]]:\\n",
            "        for res in chain:\\n",
            "            if res.id[0] == ' ' and 'CA' in res:\\n",
            "                ca = res['CA'].get_coord()\\n",
            "                dist = np.linalg.norm(ca - center)\\n",
            "                \\n",
            "                if dist <= 10.0:\\n",
            "                    protein = 'NSP12' if chain.id == chain_assignments['nsp12'] else 'NSP8'\\n",
            "                    interface_12_8_full.append({\\n",
            "                        'Protein': protein,\\n",
            "                        'Chain': chain.id,\\n",
            "                        'Residue': res.get_resname(),\\n",
            "                        'PDB_Num': res.id[1],\\n",
            "                        'Distance': dist\\n",
            "                    })\\n",
            "    \\n",
            "    df_interface_12_8 = pd.DataFrame(interface_12_8_full).sort_values('Distance').reset_index(drop=True)\\n",
            "    \\n",
            "    print(f\\\"✓ Found {len(df_interface_12_8)} interface residues:\\\")\\n",
            "    print()\\n",
            "    print(df_interface_12_8.to_string(index=False))\\n",
            "    print()\\n",
            "    \\n",
            "    nsp12_count_8 = len(df_interface_12_8[df_interface_12_8['Protein'] == 'NSP12'])\\n",
            "    nsp8_count = len(df_interface_12_8[df_interface_12_8['Protein'] == 'NSP8'])\\n",
            "    \\n",
            "    print(\\\"Interface composition:\\\")\\n",
            "    print(f\\\"  NSP12: {nsp12_count_8} residues\\\")\\n",
            "    print(f\\\"  NSP8: {nsp8_count} residues\\\")\\n",
            "    print(f\\\"  Total: {len(df_interface_12_8)} residues\\\")\\n",
            "else:\\n",
            "    print(\\\"⚠️  Cannot map interface - missing hot spots or chains\\\")\\n",
            "    df_interface_12_8 = pd.DataFrame()\\n",
            "    nsp12_count_8 = nsp8_count = 0\\n",
            "\\n",
            "print()\\n",
            "print(\\\"=\\\"*70)"
        ]
    }
]

# Add Sections 7-10 in next continuation...
print("Adding Sections 5A, 5B, 6A, 6B...")

# Add cells to notebook
notebook['cells'].extend(new_cells)

# Save
with open('05_7DFG_interface_discovery.ipynb', 'w') as f:
    json.dump(notebook, f, indent=1)

print(f"✓ Added Sections 5A, 5B, 6A, 6B")
print(f"  Total cells now: {len(notebook['cells'])}")
print()
print("Next: Run add_sections_7_10.py to complete the notebook")
