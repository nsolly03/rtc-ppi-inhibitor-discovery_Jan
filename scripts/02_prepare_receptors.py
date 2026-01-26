#!/usr/bin/env python3
"""
Prepare receptor structures for docking
- Extract relevant chains
- Convert to PDBQT format
- Add hydrogens and charges

Author: Olivier Nsekuye
Lab: GIGA-VIN, University of Liège
Date: January 2025
"""

import os
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
import subprocess


class ChainSelector(Select):
    """Select specific chains from PDB structure"""
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def extract_chains(pdb_file, chain_ids, output_file):
    """
    Extract specific chains from PDB structure.
    
    Args:
        pdb_file (str): Input PDB file
        chain_ids (list): List of chain IDs to extract (e.g., ['A', 'B'])
        output_file (str): Output PDB file
    """
    print(f"\n🔍 Extracting chains {chain_ids} from {os.path.basename(pdb_file)}...")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, ChainSelector(chain_ids))
    
    # Count atoms in extracted structure
    with open(output_file, 'r') as f:
        atoms = [line for line in f if line.startswith('ATOM')]
    
    print(f"   ✅ Extracted {len(atoms)} atoms")
    print(f"   ✅ Saved to {output_file}")


def prepare_receptor_meeko(pdb_file, output_pdbqt):
    """
    Prepare receptor using Meeko (adds hydrogens, charges).
    
    Args:
        pdb_file (str): Input PDB file
        output_pdbqt (str): Output PDBQT file
        
    Returns:
        bool: True if successful
    """
    print(f"\n⚙️  Preparing receptor with Meeko...")
    print(f"   Adding hydrogens and charges...")
    
    # Use mk_prepare_receptor.py from Meeko
    cmd = f"mk_prepare_receptor.py -i {pdb_file} -o {output_pdbqt}"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"   ✅ Receptor prepared successfully")
        return True
    else:
        print(f"   ⚠️  Meeko preparation encountered issues:")
        if result.stderr:
            print(f"   {result.stderr[:200]}")
        return False


def get_receptor_info(pdbqt_file):
    """
    Show information about prepared receptor.
    
    Args:
        pdbqt_file (str): PDBQT file path
    """
    if not os.path.exists(pdbqt_file):
        print(f"\n   ⚠️  File not found: {pdbqt_file}")
        return
    
    size_kb = os.path.getsize(pdbqt_file) / 1024
    
    # Count atoms and get center
    with open(pdbqt_file, 'r') as f:
        lines = f.readlines()
    
    atom_lines = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]
    
    # Calculate center of mass (approximate)
    if atom_lines:
        x_coords = [float(line[30:38]) for line in atom_lines]
        y_coords = [float(line[38:46]) for line in atom_lines]
        z_coords = [float(line[46:54]) for line in atom_lines]
        
        center_x = sum(x_coords) / len(x_coords)
        center_y = sum(y_coords) / len(y_coords)
        center_z = sum(z_coords) / len(z_coords)
    
    print(f"\n📊 Receptor Information:")
    print(f"   File: {os.path.basename(pdbqt_file)}")
    print(f"   Size: {size_kb:.1f} KB")
    print(f"   Atoms: {len(atom_lines)}")
    print(f"   Approximate center: ({center_x:.1f}, {center_y:.1f}, {center_z:.1f})")


def main():
    """Main execution function"""
    
    print("=" * 70)
    print("🧬 RECEPTOR PREPARATION FOR DOCKING")
    print("   Author: Olivier Nsekuye")
    print("   Lab: GIGA-VIN, University of Liège")
    print("=" * 70)
    
    # Define which chains to use for each target
    # Based on your FRIA proposal - targeting PPI interfaces
    targets = {
        '7DFG': {
            'description': 'NSP12-NSP7-NSP8 (RdRp complex)',
            'chains': ['A', 'B', 'C'],  # NSP12 (A), NSP7 (C), NSP8 (B)
            'note': 'Main target - contains NSP12-NSP7, NSP12-NSP8, NSP7-NSP8 interfaces'
        },
        '6W4H': {
            'description': 'NSP10-NSP16 (2-O-MTase complex)',
            'chains': ['A', 'B'],  # NSP16 (A), NSP10 (B)
            'note': 'Target NSP10-NSP16 interface'
        },
        '7EDI': {
            'description': 'NSP10-NSP14 (ExoN complex)',
            'chains': ['A'],  # NSP14 with NSP10 (selecting first chain)
            'note': 'Target NSP10-NSP14 interface'
        }
    }
    
    input_dir = 'data/targets'
    output_dir = 'data/docking/receptors'
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    successful = 0
    failed = 0
    
    for pdb_id, info in targets.items():
        print(f"\n{'='*70}")
        print(f"Processing: {pdb_id} - {info['description']}")
        print(f"Note: {info['note']}")
        print('='*70)
        
        try:
            # Input file (cleaned PDB)
            clean_pdb = f"{input_dir}/{pdb_id}_clean.pdb"
            
            if not os.path.exists(clean_pdb):
                print(f"   ⚠️  File not found: {clean_pdb}")
                print(f"   Run 01_download_pdb_structures.py first")
                failed += 1
                continue
            
            # Extract relevant chains
            extracted_pdb = f"{output_dir}/{pdb_id}_extracted.pdb"
            extract_chains(clean_pdb, info['chains'], extracted_pdb)
            
            # Prepare receptor (convert to PDBQT)
            output_pdbqt = f"{output_dir}/{pdb_id}_receptor.pdbqt"
            success = prepare_receptor_meeko(extracted_pdb, output_pdbqt)
            
            if success:
                # Show info about prepared receptor
                get_receptor_info(output_pdbqt)
                print(f"\n✅ {pdb_id} receptor preparation complete!")
                successful += 1
            else:
                print(f"\n⚠️  {pdb_id} preparation had issues")
                failed += 1
            
        except Exception as e:
            print(f"\n❌ Error processing {pdb_id}: {e}")
            failed += 1
            continue
    
    print("\n" + "=" * 70)
    print(f"🎉 RECEPTOR PREPARATION COMPLETE!")
    print("=" * 70)
    print(f"\n📊 Summary:")
    print(f"   ✅ Successful: {successful}/{successful+failed}")
    print(f"   ❌ Failed: {failed}/{successful+failed}")
    print(f"\n📁 Receptors saved in: {output_dir}/")
    print(f"\n📋 Next steps:")
    print(f"   1. Identify binding pockets (manual or fpocket)")
    print(f"   2. Define docking grid boxes for each interface")
    print(f"   3. Prepare ligand libraries for screening")
    print(f"\n💡 Tips:")
    print(f"   - Use PyMOL or Chimera to visualize receptors")
    print(f"   - Check the approximate centers for grid box placement")
    print(f"   - For VirtualFlow: use these PDBQT files as receptors")
    

if __name__ == "__main__":
    main()
