#!/usr/bin/env python3
"""
Comprehensive search for ALL Lys-Asp pairs in 6W4H NSP10-NSP16 interface
Check if there are other potential hot spots beyond K76-D107
"""

from Bio import PDB
import numpy as np
import pandas as pd

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('6W4H', 'data/structures/pdb/6W4H.pdb')
model = structure[0]

nsp10_chain = 'B'
nsp16_chain = 'A'
nsp10_start = 4271
nsp16_start = 6798

print("="*80)
print("COMPREHENSIVE LYS-ASP PAIR ANALYSIS")
print("="*80)
print()

# Get all lysines in NSP10
print("Step 1: Finding all LYSINES in NSP10...")
nsp10_lysines = []
for res in model[nsp10_chain]:
    if res.id[0] == ' ' and res.get_resname() == 'LYS':
        if 'CA' in res:
            pdb_num = res.id[1]
            seq_pos = pdb_num - nsp10_start + 1
            nsp10_lysines.append({
                'residue': res,
                'pdb_num': pdb_num,
                'seq_pos': seq_pos,
                'ca_coord': res['CA'].get_coord(),
                'label': f"K{seq_pos}"
            })

print(f"Found {len(nsp10_lysines)} lysines in NSP10:")
for lys in nsp10_lysines:
    print(f"  K{lys['seq_pos']} (PDB {lys['pdb_num']})")
print()

# Get all aspartates in NSP16
print("Step 2: Finding all ASPARTATES in NSP16...")
nsp16_asps = []
for res in model[nsp16_chain]:
    if res.id[0] == ' ' and res.get_resname() == 'ASP':
        if 'CA' in res:
            pdb_num = res.id[1]
            seq_pos = pdb_num - nsp16_start + 1
            nsp16_asps.append({
                'residue': res,
                'pdb_num': pdb_num,
                'seq_pos': seq_pos,
                'ca_coord': res['CA'].get_coord(),
                'label': f"D{seq_pos}"
            })

print(f"Found {len(nsp16_asps)} aspartates in NSP16:")
for asp in nsp16_asps:
    print(f"  D{asp['seq_pos']} (PDB {asp['pdb_num']})")
print()

# Calculate ALL pairwise distances
print("="*80)
print("Step 3: Calculating ALL Lys-Asp distances...")
print("="*80)
print()

all_pairs = []
for lys in nsp10_lysines:
    for asp in nsp16_asps:
        distance = np.linalg.norm(lys['ca_coord'] - asp['ca_coord'])
        all_pairs.append({
            'nsp10_lys': lys['label'],
            'nsp10_pdb': lys['pdb_num'],
            'nsp10_seq': lys['seq_pos'],
            'nsp16_asp': asp['label'],
            'nsp16_pdb': asp['pdb_num'],
            'nsp16_seq': asp['seq_pos'],
            'distance': distance
        })

# Sort by distance
all_pairs.sort(key=lambda x: x['distance'])

# Show ALL pairs
print(f"Total possible Lys-Asp pairs: {len(all_pairs)}")
print()
print("="*80)
print("ALL LYS-ASP PAIRS (sorted by distance)")
print("="*80)
print()
print(f"{'NSP10 Lys':<12}{'PDB':<8}{'NSP16 Asp':<12}{'PDB':<8}{'Distance (Å)':<15}{'Category'}")
print("-"*80)

for pair in all_pairs:
    # Categorize
    if pair['distance'] < 5.0:
        category = "★★★ SALT BRIDGE (STRONG)"
    elif pair['distance'] < 7.0:
        category = "★★ SALT BRIDGE (LIKELY)"
    elif pair['distance'] < 10.0:
        category = "★ H-BOND (POSSIBLE)"
    elif pair['distance'] < 15.0:
        category = "Close (water-mediated?)"
    else:
        category = "Distant"
    
    # Highlight if this is K76-D107
    marker = "→→→" if (pair['nsp10_seq'] == 76 and pair['nsp16_seq'] == 107) else "   "
    
    print(f"{marker} {pair['nsp10_lys']:<12}{pair['nsp10_pdb']:<8}"
          f"{pair['nsp16_asp']:<12}{pair['nsp16_pdb']:<8}"
          f"{pair['distance']:<15.2f}{category}")

print()
print("="*80)
print("INTERFACE PAIRS (Distance < 10 Å)")
print("="*80)
print()

interface_pairs = [p for p in all_pairs if p['distance'] < 10.0]
print(f"Found {len(interface_pairs)} Lys-Asp pairs within 10 Å")
print()

if len(interface_pairs) > 0:
    print(f"{'Rank':<6}{'Pair':<20}{'Distance (Å)':<15}{'Assessment'}")
    print("-"*70)
    
    for i, pair in enumerate(interface_pairs, 1):
        pair_label = f"{pair['nsp10_lys']}-{pair['nsp16_asp']}"
        
        if pair['distance'] < 5.0:
            assessment = "VERY STRONG - Primary hot spot"
        elif pair['distance'] < 7.0:
            assessment = "STRONG - Secondary hot spot"
        else:
            assessment = "Weak - Tertiary contact"
        
        marker = "★" if (pair['nsp10_seq'] == 76 and pair['nsp16_seq'] == 107) else " "
        
        print(f"{marker} {i:<6}{pair_label:<20}{pair['distance']:<15.2f}{assessment}")

print()
print("="*80)
print("ANALYSIS SUMMARY")
print("="*80)
print()

# Count by category
salt_bridge_strong = len([p for p in all_pairs if p['distance'] < 5.0])
salt_bridge_likely = len([p for p in all_pairs if 5.0 <= p['distance'] < 7.0])
hbond_possible = len([p for p in all_pairs if 7.0 <= p['distance'] < 10.0])

print(f"Salt bridges (< 5.0 Å):     {salt_bridge_strong} pairs")
print(f"Salt bridges (5-7 Å):       {salt_bridge_likely} pairs")
print(f"H-bonds possible (7-10 Å):  {hbond_possible} pairs")
print(f"Total interface pairs:      {len(interface_pairs)} pairs")
print()

if salt_bridge_strong > 0:
    print("STRONGEST INTERACTIONS:")
    for pair in all_pairs[:salt_bridge_strong]:
        print(f"  {pair['nsp10_lys']}-{pair['nsp16_asp']}: {pair['distance']:.2f} Å")
    print()

# Check if K76-D107 is the best
k76_d107 = [p for p in all_pairs if p['nsp10_seq'] == 76 and p['nsp16_seq'] == 107][0]
rank = all_pairs.index(k76_d107) + 1

print(f"K76-D107 (paper's reference as K93-D106):")
print(f"  Distance: {k76_d107['distance']:.2f} Å")
print(f"  Rank: #{rank} out of {len(all_pairs)} total pairs")
print(f"  Rank in interface: #{interface_pairs.index(k76_d107) + 1} out of {len(interface_pairs)}")

if rank == 1:
    print(f"  ✓ This IS the strongest Lys-Asp interaction!")
else:
    strongest = all_pairs[0]
    print(f"  ⚠ NOTE: Strongest pair is {strongest['nsp10_lys']}-{strongest['nsp16_asp']} ({strongest['distance']:.2f} Å)")

print()
print("="*80)
print("RECOMMENDATION")
print("="*80)
print()

if salt_bridge_strong > 1:
    print(f"Multiple strong salt bridges detected ({salt_bridge_strong} pairs < 5 Å)")
    print("Consider targeting ALL of these in docking:")
    for pair in [p for p in all_pairs if p['distance'] < 5.0]:
        print(f"  • {pair['nsp10_lys']}-{pair['nsp16_asp']}: {pair['distance']:.2f} Å")
    print()
    print("Strategy: Use larger grid box to encompass all hot spots")
elif salt_bridge_strong == 1:
    print("Single dominant salt bridge detected")
    print(f"Primary target: {all_pairs[0]['nsp10_lys']}-{all_pairs[0]['nsp16_asp']} ({all_pairs[0]['distance']:.2f} Å)")
else:
    print("No strong salt bridges (< 5 Å) detected")
    if salt_bridge_likely > 0:
        print(f"Weaker salt bridges present ({salt_bridge_likely} pairs at 5-7 Å)")

print()

# Save detailed results
import csv
csv_file = 'data/analysis_results/6W4H_all_lys_asp_pairs.csv'
with open(csv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['nsp10_lys', 'nsp10_pdb', 'nsp10_seq', 
                                            'nsp16_asp', 'nsp16_pdb', 'nsp16_seq', 'distance'])
    writer.writeheader()
    writer.writerows(all_pairs)

print(f"✓ Detailed results saved to: {csv_file}")
print()
