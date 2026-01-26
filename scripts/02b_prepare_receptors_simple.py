#!/usr/bin/env python3
"""
Simple receptor preparation - converts PDB to PDBQT format
Author: Olivier Nsekuye
"""

import os
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select


class ChainSelector(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def extract_chains(pdb_file, chain_ids, output_file):
    print(f"\nExtracting chains {chain_ids}...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, ChainSelector(chain_ids))
    print(f"   Saved to {output_file}")


def prepare_receptor_simple(pdb_file, output_pdbqt):
    """Convert PDB to basic PDBQT format"""
    print(f"\nPreparing receptor...")
    
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()
    
    with open(output_pdbqt, 'w') as out:
        out.write("REMARK  Receptor prepared for docking\n")
        out.write(f"REMARK  Source: {os.path.basename(pdb_file)}\n")
        
        for line in pdb_lines:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                
                # Assign basic AutoDock atom type
                if atom_name.startswith('C'):
                    ad_type = 'C'
                elif atom_name.startswith('N'):
                    ad_type = 'N'
                elif atom_name.startswith('O'):
                    ad_type = 'O'
                elif atom_name.startswith('S'):
                    ad_type = 'S'
                else:
                    ad_type = 'C'
                
                out.write(line.rstrip() + f"  0.00  0.00    {ad_type}\n")
    
    print(f"   Saved to {output_pdbqt}")


def get_receptor_info(pdbqt_file):
    if not os.path.exists(pdbqt_file):
        print(f"   File not found: {pdbqt_file}")
        return
    
    with open(pdbqt_file, 'r') as f:
        atom_lines = [l for l in f if l.startswith('ATOM')]
    
    if atom_lines:
        x_coords = [float(line[30:38]) for line in atom_lines]
        y_coords = [float(line[38:46]) for line in atom_lines]
        z_coords = [float(line[46:54]) for line in atom_lines]
        
        center_x = sum(x_coords) / len(x_coords)
        center_y = sum(y_coords) / len(y_coords)
        center_z = sum(z_coords) / len(z_coords)
        
        print(f"\nReceptor Information:")
        print(f"   Atoms: {len(atom_lines)}")
        print(f"   Center: ({center_x:.1f}, {center_y:.1f}, {center_z:.1f})")


def main():
    print("=" * 70)
    print("RECEPTOR PREPARATION")
    print("=" * 70)
    
    targets = {
        '7DFG': {
            'description': 'NSP12-NSP7-NSP8 (RdRp complex)',
            'chains': ['A', 'B', 'C']
        },
        '6XEZ': {
            'description': 'NSP12-NSP7-NSP8 (alternative conformation)',
            'chains': ['A', 'B', 'C']
        },
        '6W4H': {
            'description': 'NSP10-NSP16 (2-O-MTase)',
            'chains': ['A', 'B']
        },
        '7EDI': {
            'description': 'NSP10-NSP14 (ExoN complex)',
            'chains': ['A', 'B']
        },
        '6W9C': {
            'description': 'NSP13 (Helicase)',
            'chains': ['A', 'B', 'C']
        }
    }
    
    input_dir = 'data/targets'
    output_dir = 'data/docking/receptors'
    
    success = 0
    
    for pdb_id, info in targets.items():
        print(f"\n{'='*70}")
        print(f"Processing: {pdb_id} - {info['description']}")
        print('='*70)
        
        try:
            clean_pdb = f"{input_dir}/{pdb_id}_clean.pdb"
            
            if not os.path.exists(clean_pdb):
                print(f"   File not found: {clean_pdb}")
                continue
            
            extracted_pdb = f"{output_dir}/{pdb_id}_extracted.pdb"
            extract_chains(clean_pdb, info['chains'], extracted_pdb)
            
            output_pdbqt = f"{output_dir}/{pdb_id}_receptor.pdbqt"
            prepare_receptor_simple(extracted_pdb, output_pdbqt)
            
            get_receptor_info(output_pdbqt)
            
            print(f"\nCompleted {pdb_id}!")
            success += 1
            
        except Exception as e:
            print(f"\nError: {e}")
    
    print("\n" + "=" * 70)
    print(f"COMPLETE! Prepared {success}/5 receptors")
    print("=" * 70)


if __name__ == "__main__":
    main()
