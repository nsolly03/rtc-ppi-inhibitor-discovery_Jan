#!/usr/bin/env python3
"""Check actual residue numbering in 6W4H"""

from Bio import PDB

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('6W4H', 'data/structures/pdb/6W4H.pdb')
model = structure[0]

# NSP10 = Chain B (116 residues)
# NSP16 = Chain A (299 residues)

print("="*70)
print("CHAIN B (NSP10) RESIDUES")
print("="*70)

chain_b = model['B']
residues_b = [r for r in chain_b.get_residues() if r.id[0] == ' ']

print(f"Total residues: {len(residues_b)}")
print()
print("First 10 residues:")
for i, res in enumerate(residues_b[:10]):
    print(f"  {i+1}. Residue number: {res.id[1]}, Name: {res.get_resname()}")

print()
print("Last 10 residues:")
for i, res in enumerate(residues_b[-10:]):
    print(f"  {len(residues_b)-10+i+1}. Residue number: {res.id[1]}, Name: {res.get_resname()}")

print()
print("="*70)
print("CHAIN A (NSP16) RESIDUES")
print("="*70)

chain_a = model['A']
residues_a = [r for r in chain_a.get_residues() if r.id[0] == ' ']

print(f"Total residues: {len(residues_a)}")
print()
print("First 10 residues:")
for i, res in enumerate(residues_a[:10]):
    print(f"  {i+1}. Residue number: {res.id[1]}, Name: {res.get_resname()}")

print()
print("Last 10 residues:")
for i, res in enumerate(residues_a[-10:]):
    print(f"  {len(residues_a)-10+i+1}. Residue number: {res.id[1]}, Name: {res.get_resname()}")

print()

# Search for Lysine in Chain B (NSP10)
print("="*70)
print("SEARCHING FOR LYSINE IN CHAIN B (NSP10)")
print("="*70)

lysines = []
for res in residues_b:
    if res.get_resname() == 'LYS':
        lysines.append((res.id[1], res))

print(f"Found {len(lysines)} lysine residues in NSP10:")
for resnum, res in lysines:
    print(f"  Residue number: {resnum}, Name: LYS")

print()

# Search for Aspartate in Chain A (NSP16)
print("="*70)
print("SEARCHING FOR ASPARTATE IN CHAIN A (NSP16)")
print("="*70)

asps = []
for res in residues_a:
    if res.get_resname() == 'ASP':
        asps.append((res.id[1], res))

print(f"Found {len(asps)} aspartate residues in NSP16:")
for resnum, res in asps:
    print(f"  Residue number: {resnum}, Name: ASP")
