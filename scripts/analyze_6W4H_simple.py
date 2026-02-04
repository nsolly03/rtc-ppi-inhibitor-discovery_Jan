#!/usr/bin/env python3
"""
Simple analysis of 6W4H NSP10-NSP16 structure
Extract hot spot coordinates based on Trepte et al. 2024
"""

from Bio import PDB
import numpy as np
import json
import os

print("="*70)
print("6W4H NSP10-NSP16 INTERFACE ANALYSIS")
print("="*70)
print()

# Load structure
pdb_file = 'data/structures/pdb/6W4H.pdb'
if not os.path.exists(pdb_file):
    print(f"ERROR: {pdb_file} not found!")
    exit(1)

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('6W4H', pdb_file)
model = structure[0]

# Identify chains (NSP10 is shorter, NSP16 is longer)
chain_info = {}
for chain in model:
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    chain_info[chain.id] = len(residues)
    print(f"Chain {chain.id}: {len(residues)} residues")

print()

nsp10_chain = min(chain_info, key=chain_info.get)
nsp16_chain = max(chain_info, key=chain_info.get)

print(f"NSP10 = Chain {nsp10_chain} ({chain_info[nsp10_chain]} residues)")
print(f"NSP16 = Chain {nsp16_chain} ({chain_info[nsp16_chain]} residues)")
print()

# Get hot spot residues
print("="*70)
print("HOT SPOT RESIDUES")
print("="*70)
print()

try:
    lys93 = model[nsp10_chain][93]
    print(f"NSP10 Lys93 (Chain {nsp10_chain}):")
    print(f"  Residue: {lys93.get_resname()}")
    
    if 'CA' in lys93:
        lys93_ca = lys93['CA'].get_coord()
        print(f"  CA: ({lys93_ca[0]:.3f}, {lys93_ca[1]:.3f}, {lys93_ca[2]:.3f})")
        
        # Center of mass
        atoms = [a.get_coord() for a in lys93.get_atoms()]
        center_lys93 = np.mean(atoms, axis=0)
        print(f"  Center: ({center_lys93[0]:.3f}, {center_lys93[1]:.3f}, {center_lys93[2]:.3f})")
    print()
    
    asp106 = model[nsp16_chain][106]
    print(f"NSP16 Asp106 (Chain {nsp16_chain}):")
    print(f"  Residue: {asp106.get_resname()}")
    
    if 'CA' in asp106:
        asp106_ca = asp106['CA'].get_coord()
        print(f"  CA: ({asp106_ca[0]:.3f}, {asp106_ca[1]:.3f}, {asp106_ca[2]:.3f})")
    print()
    
    # Distance
    distance = np.linalg.norm(lys93_ca - asp106_ca)
    print(f"Distance between hot spots: {distance:.2f} Å")
    print()
    
    if distance < 5.0:
        print("✓ Close proximity - salt bridge likely")
    elif distance < 8.0:
        print("✓ Moderate proximity - interaction possible")
    else:
        print("⚠ Distant - check structure")
    print()
    
except KeyError as e:
    print(f"ERROR: Could not find residue: {e}")
    exit(1)

# Find interface residues (within 10 Å)
print("="*70)
print("INTERFACE RESIDUES (within 10 Å of Lys93)")
print("="*70)
print()

interface = []
for chain in model:
    for res in chain:
        if res.id[0] == ' ' and 'CA' in res:
            ca = res['CA'].get_coord()
            dist = np.linalg.norm(ca - lys93_ca)
            if dist <= 10.0:
                interface.append({
                    'chain': chain.id,
                    'residue': res.get_resname(),
                    'number': res.id[1],
                    'distance': dist
                })

interface.sort(key=lambda x: x['distance'])

print(f"Found {len(interface)} residues:")
print(f"{'Chain':<8}{'Residue':<10}{'Number':<10}{'Distance (Å)':<15}")
print("-"*45)
for r in interface[:20]:  # Show first 20
    print(f"{r['chain']:<8}{r['residue']:<10}{r['number']:<10}{r['distance']:<15.2f}")

if len(interface) > 20:
    print(f"... and {len(interface)-20} more")

print()

nsp10_count = sum(1 for r in interface if r['chain'] == nsp10_chain)
nsp16_count = sum(1 for r in interface if r['chain'] == nsp16_chain)
print(f"NSP10 interface: {nsp10_count} residues")
print(f"NSP16 interface: {nsp16_count} residues")
print()

# Grid box for docking
print("="*70)
print("DOCKING GRID BOX")
print("="*70)
print()
print("Center (Lys93 center of mass):")
print(f"  center_x = {center_lys93[0]:.3f}")
print(f"  center_y = {center_lys93[1]:.3f}")
print(f"  center_z = {center_lys93[2]:.3f}")
print()
print("Size (start with cubic box):")
print("  size_x = 25.0")
print("  size_y = 25.0")
print("  size_z = 25.0")
print()

# Save results
os.makedirs('data/analysis_results', exist_ok=True)

results = {
    'pdb_id': '6W4H',
    'nsp10_chain': nsp10_chain,
    'nsp16_chain': nsp16_chain,
    'hot_spots': {
        'lys93': {
            'chain': nsp10_chain,
            'ca_coord': lys93_ca.tolist(),
            'center': center_lys93.tolist()
        },
        'asp106': {
            'chain': nsp16_chain,
            'ca_coord': asp106_ca.tolist()
        },
        'distance': float(distance)
    },
    'grid_box': {
        'center_x': float(center_lys93[0]),
        'center_y': float(center_lys93[1]),
        'center_z': float(center_lys93[2]),
        'size_x': 25.0,
        'size_y': 25.0,
        'size_z': 25.0
    },
    'interface_residues': len(interface)
}

output_file = 'data/analysis_results/6W4H_analysis.json'
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"✓ Results saved to {output_file}")
print()
print("="*70)
print("ANALYSIS COMPLETE")
print("="*70)
