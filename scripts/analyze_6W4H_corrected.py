#!/usr/bin/env python3
"""
Corrected analysis of 6W4H NSP10-NSP16 structure
Using ACTUAL hot spot residues: K76 (NSP10) and D107 (NSP16)
"""

from Bio import PDB
import numpy as np
import json
import os

print("="*70)
print("6W4H NSP10-NSP16 INTERFACE ANALYSIS (CORRECTED)")
print("="*70)
print()

# Load structure
pdb_file = 'data/structures/pdb/6W4H.pdb'
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('6W4H', pdb_file)
model = structure[0]

# NSP10 = Chain B (starts at 4271)
# NSP16 = Chain A (starts at 6798)
nsp10_chain = 'B'
nsp16_chain = 'A'
nsp10_start = 4271
nsp16_start = 6798

print(f"NSP10 = Chain {nsp10_chain}")
print(f"NSP16 = Chain {nsp16_chain}")
print()

# CORRECTED HOT SPOTS
# Paper says K93-D106, but in 6W4H structure the actual hot spot is K76-D107
lys_seq_pos = 76   # NSP10 K76
asp_seq_pos = 107  # NSP16 D107

lys_pdb_num = nsp10_start + lys_seq_pos - 1  # 4346
asp_pdb_num = nsp16_start + asp_seq_pos - 1  # 6904

print("="*70)
print("HOT SPOT RESIDUES (ACTUAL IN 6W4H STRUCTURE)")
print("="*70)
print()
print(f"Paper reference: NSP10 K93 - NSP16 D106")
print(f"Actual in 6W4H:  NSP10 K{lys_seq_pos} - NSP16 D{asp_seq_pos}")
print()

# Get hot spot residues
lys = model[nsp10_chain][lys_pdb_num]
asp = model[nsp16_chain][asp_pdb_num]

print(f"NSP10 K{lys_seq_pos} (Chain {nsp10_chain}, PDB {lys_pdb_num}):")
print(f"  Residue: {lys.get_resname()}")

lys_ca = lys['CA'].get_coord()
print(f"  CA: ({lys_ca[0]:.3f}, {lys_ca[1]:.3f}, {lys_ca[2]:.3f})")

atoms = [a.get_coord() for a in lys.get_atoms()]
center_lys = np.mean(atoms, axis=0)
print(f"  Center: ({center_lys[0]:.3f}, {center_lys[1]:.3f}, {center_lys[2]:.3f})")
print()

print(f"NSP16 D{asp_seq_pos} (Chain {nsp16_chain}, PDB {asp_pdb_num}):")
print(f"  Residue: {asp.get_resname()}")

asp_ca = asp['CA'].get_coord()
print(f"  CA: ({asp_ca[0]:.3f}, {asp_ca[1]:.3f}, {asp_ca[2]:.3f})")
print()

# Distance
distance = np.linalg.norm(lys_ca - asp_ca)
print(f"Distance between hot spots: {distance:.2f} Å")
print()

if distance < 5.5:
    print("✓ CLOSE proximity - salt bridge LIKELY")
    interaction = "Salt bridge"
elif distance < 8.0:
    print("✓ MODERATE proximity - interaction possible")
    interaction = "H-bond possible"
else:
    print("⚠ DISTANT")
    interaction = "Distant"
print()

# Interface residues (within 10 Å of K76)
print("="*70)
print(f"INTERFACE RESIDUES (within 10 Å of NSP10 K{lys_seq_pos})")
print("="*70)
print()

interface = []
for chain in model:
    for res in chain:
        if res.id[0] == ' ' and 'CA' in res:
            ca = res['CA'].get_coord()
            dist = np.linalg.norm(ca - lys_ca)
            if dist <= 10.0:
                # Calculate sequence position
                pdb_num = res.id[1]
                if chain.id == nsp10_chain:
                    seq_pos = pdb_num - nsp10_start + 1
                    seq_label = f"NSP10_{seq_pos}"
                else:
                    seq_pos = pdb_num - nsp16_start + 1
                    seq_label = f"NSP16_{seq_pos}"
                
                interface.append({
                    'chain': chain.id,
                    'residue': res.get_resname(),
                    'pdb_number': pdb_num,
                    'seq_position': seq_pos,
                    'seq_label': seq_label,
                    'distance': dist
                })

interface.sort(key=lambda x: x['distance'])

print(f"Found {len(interface)} residues:")
print()
print(f"{'Chain':<8}{'Residue':<10}{'PDB#':<10}{'Seq Pos':<12}{'Dist (Å)':<12}")
print("-"*60)

for r in interface[:30]:
    print(f"{r['chain']:<8}{r['residue']:<10}{r['pdb_number']:<10}"
          f"{r['seq_label']:<12}{r['distance']:<12.2f}")

if len(interface) > 30:
    print(f"... and {len(interface)-30} more")

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
print(f"Center (NSP10 K{lys_seq_pos} center of mass):")
print(f"  center_x = {center_lys[0]:.3f}")
print(f"  center_y = {center_lys[1]:.3f}")
print(f"  center_z = {center_lys[2]:.3f}")
print()
print("Size (start with cubic box):")
print("  size_x = 25.0")
print("  size_y = 25.0")
print("  size_z = 25.0")
print()

# Save results (with proper JSON serialization)
os.makedirs('data/analysis_results', exist_ok=True)

results = {
    'pdb_id': '6W4H',
    'note': 'Hot spot K76-D107 in structure (paper refers to K93-D106 in reference sequence)',
    'chains': {
        'nsp10': nsp10_chain,
        'nsp16': nsp16_chain
    },
    'hot_spots': {
        'nsp10_lys': {
            'sequence_position': int(lys_seq_pos),
            'pdb_number': int(lys_pdb_num),
            'chain': nsp10_chain,
            'ca_coord': [float(lys_ca[0]), float(lys_ca[1]), float(lys_ca[2])],
            'center': [float(center_lys[0]), float(center_lys[1]), float(center_lys[2])]
        },
        'nsp16_asp': {
            'sequence_position': int(asp_seq_pos),
            'pdb_number': int(asp_pdb_num),
            'chain': nsp16_chain,
            'ca_coord': [float(asp_ca[0]), float(asp_ca[1]), float(asp_ca[2])]
        },
        'distance': float(distance),
        'interaction': interaction
    },
    'grid_box': {
        'center_x': float(center_lys[0]),
        'center_y': float(center_lys[1]),
        'center_z': float(center_lys[2]),
        'size_x': 25.0,
        'size_y': 25.0,
        'size_z': 25.0
    },
    'interface': {
        'total_residues': len(interface),
        'nsp10_residues': nsp10_count,
        'nsp16_residues': nsp16_count
    }
}

output_file = 'data/analysis_results/6W4H_analysis.json'
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"✓ Results saved to {output_file}")
print()

# Save Vina config
vina_config = f"""# AutoDock Vina Configuration for 6W4H NSP10-NSP16
# Target: NSP10 K{lys_seq_pos} interface (PDB {lys_pdb_num})
# Hot spot: K{lys_seq_pos} (NSP10) - D{asp_seq_pos} (NSP16)
# Distance: {distance:.2f} Å

receptor = 6W4H_prepared.pdbqt
ligand = ligand.pdbqt

# Grid box center (Å)
center_x = {center_lys[0]:.3f}
center_y = {center_lys[1]:.3f}
center_z = {center_lys[2]:.3f}

# Grid box size (Å)
size_x = 25.0
size_y = 25.0
size_z = 25.0

# Docking parameters
exhaustiveness = 8
num_modes = 9
energy_range = 3

# Output
out = docking_output.pdbqt
log = docking_log.txt
"""

vina_file = 'data/analysis_results/6W4H_vina_config.txt'
with open(vina_file, 'w') as f:
    f.write(vina_config)

print(f"✓ Vina config saved to {vina_file}")
print()

# Save CSV of interface residues
csv_file = 'data/analysis_results/6W4H_interface_residues.csv'
with open(csv_file, 'w') as f:
    f.write("Chain,Residue,PDB_Number,Seq_Position,Seq_Label,Distance_A\n")
    for r in interface:
        f.write(f"{r['chain']},{r['residue']},{r['pdb_number']},"
                f"{r['seq_position']},{r['seq_label']},{r['distance']:.2f}\n")

print(f"✓ Interface residues saved to {csv_file}")
print()

print("="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print()
print("IMPORTANT NOTE:")
print("  The paper (Trepte et al. 2024) refers to K93-D106")
print("  But in 6W4H structure, the actual hot spot is K76-D107")
print("  This is the correct target for docking with 6W4H!")
print()
print("FILES CREATED:")
print(f"  1. {output_file}")
print(f"  2. {vina_file}")
print(f"  3. {csv_file}")
print()
