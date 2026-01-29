#!/usr/bin/env python3
"""Find the exact hot spot residues in 6W4H"""

from Bio import PDB
import numpy as np

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('6W4H', 'data/structures/pdb/6W4H.pdb')
model = structure[0]

print("="*70)
print("FINDING HOT SPOT RESIDUES")
print("="*70)
print()

# NSP10 Chain B starts at 4271
# NSP16 Chain A starts at 6798

# Calculate expected PDB residue numbers
nsp10_start = 4271
nsp16_start = 6798

lys93_expected = nsp10_start + 92  # 93-1 = 92 (0-indexed)
asp106_expected = nsp16_start + 105  # 106-1 = 105

print(f"Expected NSP10 Lys93 at PDB residue: {lys93_expected}")
print(f"Expected NSP16 Asp106 at PDB residue: {asp106_expected}")
print()

# Check if they exist
chain_b = model['B']
chain_a = model['A']

print("="*70)
print("CHECKING NSP10 LYS93 (around residue 4363)")
print("="*70)

for resnum in range(4360, 4370):
    try:
        res = chain_b[resnum]
        if res.id[0] == ' ':
            seq_pos = resnum - nsp10_start + 1
            print(f"  PDB {resnum} = NSP10 position {seq_pos}: {res.get_resname()}")
    except KeyError:
        pass

print()

print("="*70)
print("CHECKING NSP16 ASP106 (around residue 6903)")
print("="*70)

for resnum in range(6900, 6910):
    try:
        res = chain_a[resnum]
        if res.id[0] == ' ':
            seq_pos = resnum - nsp16_start + 1
            print(f"  PDB {resnum} = NSP16 position {seq_pos}: {res.get_resname()}")
    except KeyError:
        pass

print()

# Find closest matches
print("="*70)
print("LYSINES IN NSP10 WITH SEQUENCE POSITIONS")
print("="*70)

for res in chain_b:
    if res.id[0] == ' ' and res.get_resname() == 'LYS':
        pdb_num = res.id[1]
        seq_pos = pdb_num - nsp10_start + 1
        print(f"  PDB {pdb_num} = NSP10 K{seq_pos}")

print()

print("="*70)
print("ASPARTATES IN NSP16 WITH SEQUENCE POSITIONS")
print("="*70)

for res in chain_a:
    if res.id[0] == ' ' and res.get_resname() == 'ASP':
        pdb_num = res.id[1]
        seq_pos = pdb_num - nsp16_start + 1
        print(f"  PDB {pdb_num} = NSP16 D{seq_pos}")

print()

# Calculate distances between all LYS-ASP pairs
print("="*70)
print("DISTANCES BETWEEN ALL LYS (NSP10) and ASP (NSP16) PAIRS")
print("="*70)

lysines = [(res.id[1], res) for res in chain_b if res.id[0] == ' ' and res.get_resname() == 'LYS']
asps = [(res.id[1], res) for res in chain_a if res.id[0] == ' ' and res.get_resname() == 'ASP']

closest_pairs = []

for lys_num, lys in lysines:
    lys_seq = lys_num - nsp10_start + 1
    if 'CA' in lys:
        lys_ca = lys['CA'].get_coord()
        
        for asp_num, asp in asps:
            asp_seq = asp_num - nsp16_start + 1
            if 'CA' in asp:
                asp_ca = asp['CA'].get_coord()
                distance = np.linalg.norm(lys_ca - asp_ca)
                
                # Only show close pairs (< 10 Å)
                if distance < 10.0:
                    closest_pairs.append({
                        'lys_pdb': lys_num,
                        'lys_seq': lys_seq,
                        'asp_pdb': asp_num,
                        'asp_seq': asp_seq,
                        'distance': distance
                    })

# Sort by distance
closest_pairs.sort(key=lambda x: x['distance'])

print(f"Found {len(closest_pairs)} LYS-ASP pairs within 10 Å:")
print()
print(f"{'NSP10 Lys':<12}{'NSP16 Asp':<12}{'Distance (Å)':<15}{'Notes'}")
print("-"*70)

for pair in closest_pairs[:10]:  # Show top 10
    lys_label = f"K{pair['lys_seq']}"
    asp_label = f"D{pair['asp_seq']}"
    
    notes = ""
    if pair['lys_seq'] == 93:
        notes += "★ LYS93 "
    if pair['asp_seq'] == 106:
        notes += "★ ASP106 "
    if pair['distance'] < 5.0:
        notes += "SALT BRIDGE"
    
    print(f"{lys_label:<12}{asp_label:<12}{pair['distance']:<15.2f}{notes}")
