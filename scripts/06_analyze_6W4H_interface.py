#!/usr/bin/env python3
"""
Analyze 6W4H NSP10-NSP16 structure and identify hot spot residues
Based on Trepte et al. 2024 findings:
- NSP10 Lys93 (K93)
- NSP16 Asp106 (D106)
"""

import os
from Bio import PDB
import numpy as np

def analyze_6w4h():
    """Analyze 6W4H structure and extract hot spot information."""
    
    print("=" * 80)
    print("6W4H NSP10-NSP16 STRUCTURE ANALYSIS")
    print("=" * 80)
    print()
    
    # Parse PDB structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('6W4H', 'data/structures/pdb/6W4H.pdb')
    
    # Get model (usually only one)
    model = structure[0]
    
    # Print basic information
    print("STRUCTURE INFORMATION:")
    print("-" * 80)
    chains = list(model.get_chains())
    print(f"Number of chains: {len(chains)}")
    print()
    
    chain_info = {}
    for chain in chains:
        chain_id = chain.id
        residues = list(chain.get_residues())
        num_residues = len(residues)
        
        # Get residue range
        res_numbers = [res.id[1] for res in residues if res.id[0] == ' ']
        if res_numbers:
            res_range = f"{min(res_numbers)}-{max(res_numbers)}"
        else:
            res_range = "N/A"
        
        chain_info[chain_id] = {
            'num_residues': num_residues,
            'range': res_range,
            'residues': residues
        }
        
        print(f"Chain {chain_id}:")
        print(f"  Number of residues: {num_residues}")
        print(f"  Residue range: {res_range}")
        print()
    
    # Identify which chain is NSP10 and which is NSP16
    # NSP10 is shorter (~139 residues), NSP16 is longer (~298 residues)
    nsp10_chain = None
    nsp16_chain = None
    
    for chain_id, info in chain_info.items():
        if info['num_residues'] < 200:
            nsp10_chain = chain_id
            print(f"Chain {chain_id} identified as NSP10 (shorter chain)")
        else:
            nsp16_chain = chain_id
            print(f"Chain {chain_id} identified as NSP16 (longer chain)")
    
    print()
    print("=" * 80)
    print("HOT SPOT RESIDUES (from Trepte et al. 2024)")
    print("=" * 80)
    print()
    
    # Find NSP10 Lys93
    if nsp10_chain:
        chain = model[nsp10_chain]
        try:
            lys93 = chain[93]
            print(f"NSP10 Lys93 (Chain {nsp10_chain}, residue 93):")
            print(f"  Residue name: {lys93.get_resname()}")
            
            # Get CA (alpha carbon) coordinates
            if 'CA' in lys93:
                ca_coord = lys93['CA'].get_coord()
                print(f"  CA coordinates: X={ca_coord[0]:.3f}, Y={ca_coord[1]:.3f}, Z={ca_coord[2]:.3f}")
            
            # Get all atom coordinates
            atoms = list(lys93.get_atoms())
            print(f"  Number of atoms: {len(atoms)}")
            
            # Calculate center of mass
            coords = np.array([atom.get_coord() for atom in atoms])
            center = coords.mean(axis=0)
            print(f"  Center of mass: X={center[0]:.3f}, Y={center[1]:.3f}, Z={center[2]:.3f}")
            print()
            
        except KeyError:
            print(f"WARNING: Lys93 not found in chain {nsp10_chain}")
            print()
    
    # Find NSP16 Asp106
    if nsp16_chain:
        chain = model[nsp16_chain]
        try:
            asp106 = chain[106]
            print(f"NSP16 Asp106 (Chain {nsp16_chain}, residue 106):")
            print(f"  Residue name: {asp106.get_resname()}")
            
            # Get CA coordinates
            if 'CA' in asp106:
                ca_coord = asp106['CA'].get_coord()
                print(f"  CA coordinates: X={ca_coord[0]:.3f}, Y={ca_coord[1]:.3f}, Z={ca_coord[2]:.3f}")
            
            # Get all atom coordinates
            atoms = list(asp106.get_atoms())
            print(f"  Number of atoms: {len(atoms)}")
            
            # Calculate center of mass
            coords = np.array([atom.get_coord() for atom in atoms])
            center = coords.mean(axis=0)
            print(f"  Center of mass: X={center[0]:.3f}, Y={center[1]:.3f}, Z={center[2]:.3f}")
            print()
            
        except KeyError:
            print(f"WARNING: Asp106 not found in chain {nsp16_chain}")
            print()
    
    # Calculate distance between hot spots
    if nsp10_chain and nsp16_chain:
        try:
            lys93 = model[nsp10_chain][93]
            asp106 = model[nsp16_chain][106]
            
            if 'CA' in lys93 and 'CA' in asp106:
                lys93_ca = lys93['CA'].get_coord()
                asp106_ca = asp106['CA'].get_coord()
                distance = np.linalg.norm(lys93_ca - asp106_ca)
                
                print("=" * 80)
                print("HOT SPOT INTERACTION DISTANCE")
                print("=" * 80)
                print(f"Distance between Lys93 (NSP10) and Asp106 (NSP16): {distance:.2f} Å")
                print()
                
                if distance < 5.0:
                    print("✓ Hot spots are in CLOSE proximity (salt bridge likely)")
                elif distance < 8.0:
                    print("✓ Hot spots are in MODERATE proximity (hydrogen bond possible)")
                else:
                    print("⚠ Hot spots are DISTANT (check chain assignment)")
                print()
        except KeyError:
            print("Could not calculate distance - residues not found")
            print()
    
    # Find residues within 10 Å of Lys93 (for grid box definition)
    print("=" * 80)
    print("INTERFACE RESIDUES NEAR LYS93 (within 10 Å)")
    print("=" * 80)
    print()
    
    if nsp10_chain:
        try:
            lys93 = model[nsp10_chain][93]
            lys93_ca = lys93['CA'].get_coord()
            
            nearby_residues = []
            
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard residue
                        if 'CA' in residue:
                            ca = residue['CA'].get_coord()
                            distance = np.linalg.norm(ca - lys93_ca)
                            
                            if distance <= 10.0:
                                nearby_residues.append({
                                    'chain': chain.id,
                                    'resname': residue.get_resname(),
                                    'resnum': residue.id[1],
                                    'distance': distance
                                })
            
            # Sort by distance
            nearby_residues.sort(key=lambda x: x['distance'])
            
            print(f"Found {len(nearby_residues)} residues within 10 Å of Lys93:")
            print()
            print(f"{'Chain':<8}{'Residue':<12}{'Number':<10}{'Distance (Å)':<15}")
            print("-" * 50)
            
            for res in nearby_residues:
                print(f"{res['chain']:<8}{res['resname']:<12}{res['resnum']:<10}{res['distance']:<15.2f}")
            
            print()
            
            # Separate by chain
            nsp10_nearby = [r for r in nearby_residues if r['chain'] == nsp10_chain]
            nsp16_nearby = [r for r in nearby_residues if r['chain'] == nsp16_chain]
            
            print(f"NSP10 (Chain {nsp10_chain}): {len(nsp10_nearby)} residues")
            print(f"NSP16 (Chain {nsp16_chain}): {len(nsp16_nearby)} residues")
            print()
            
        except KeyError:
            print("Could not find nearby residues - Lys93 not found")
            print()
    
    # Calculate center for docking grid box
    print("=" * 80)
    print("SUGGESTED DOCKING GRID BOX")
    print("=" * 80)
    print()
    
    if nsp10_chain:
        try:
            lys93 = model[nsp10_chain][93]
            atoms = list(lys93.get_atoms())
            coords = np.array([atom.get_coord() for atom in atoms])
            center = coords.mean(axis=0)
            
            print("Grid box center (based on Lys93 center of mass):")
            print(f"  center_x = {center[0]:.3f}")
            print(f"  center_y = {center[1]:.3f}")
            print(f"  center_z = {center[2]:.3f}")
            print()
            print("Suggested grid box size (from Trepte et al. 2024):")
            print("  size_x = 25")
            print("  size_y = 25")
            print("  size_z = 25")
            print()
            print("Note: Trepte et al. used asymmetric box: 75.6 × 16.8 × 17.6 Å")
            print("Start with 25 Å cubic box, adjust if needed")
            print()
            
        except KeyError:
            print("Could not calculate grid box center - Lys93 not found")
            print()
    
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print("Next steps:")
    print("1. Visualize structure in PyMOL to confirm hot spots")
    print("2. Use grid box coordinates for docking (Week 5+)")
    print("3. Proceed to pocket identification with fpocket (Week 3-4)")
    print()

if __name__ == "__main__":
    analyze_6w4h()
