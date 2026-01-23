#!/usr/bin/env python3
"""
Download and prepare PDB structures for RTC-PPI docking
Author: Olivier Nsekuye
Lab: GIGA-VIN, University of Liège
Date: January 2025
"""

import os
import requests
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select


class ProteinCleaner(Select):
    """Remove water molecules and heteroatoms from PDB structure"""
    def accept_residue(self, residue):
        # Only keep standard amino acid residues
        return residue.id[0] == ' '


def download_pdb(pdb_id, output_dir):
    """
    Download PDB structure from RCSB PDB database.
    
    Args:
        pdb_id (str): 4-letter PDB code
        output_dir (str): Directory to save PDB file
        
    Returns:
        str: Path to downloaded PDB file
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = f"{output_dir}/{pdb_id}.pdb"
    
    print(f"\n📥 Downloading {pdb_id} from RCSB PDB...")
    print(f"   URL: {url}")
    
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(output_path, 'w') as f:
            f.write(response.text)
        print(f"   ✅ Saved to {output_path}")
        return output_path
    else:
        raise Exception(f"   ❌ Failed to download {pdb_id} (Status: {response.status_code})")


def clean_structure(pdb_file, output_file):
    """
    Remove water molecules and heteroatoms from PDB structure.
    
    Args:
        pdb_file (str): Input PDB file
        output_file (str): Output cleaned PDB file
    """
    print(f"\n🧹 Cleaning structure...")
    print(f"   Removing: water molecules, ions, ligands")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, ProteinCleaner())
    
    print(f"   ✅ Cleaned structure saved to {output_file}")


def get_structure_info(pdb_file):
    """
    Extract basic information from PDB structure.
    
    Args:
        pdb_file (str): PDB file path
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    print(f"\n📊 Structure Information:")
    
    # Count chains and residues
    chains = list(structure.get_chains())
    print(f"   Chains: {len(chains)}")
    
    for chain in chains:
        residues = list(chain.get_residues())
        print(f"   Chain {chain.id}: {len(residues)} residues")


def main():
    """Main execution function"""
    
    print("=" * 70)
    print("🧬 RTC-PPI STRUCTURE PREPARATION")
    print("   Author: Olivier Nsekuye")
    print("   Lab: GIGA-VIN, University of Liège")
    print("=" * 70)
    
    # Define your target structures
    # These are the NSP complexes mentioned in your FRIA proposal
    targets = {
        '7DFG': 'NSP12-NSP7-NSP8 (RdRp complex - main target)',
        '6XEZ': 'NSP12-NSP7-NSP8 (alternative conformation)',
        '7EDI': 'NSP10-NSP14 (ExoN proofreading complex)',
        '6W4H': 'NSP10-NSP16 (2-O-MTase complex)',
        '6W9C': 'NSP13 (Helicase)'
    }
    
    output_dir = 'data/targets'
    
    for pdb_id, description in targets.items():
        print(f"\n{'='*70}")
        print(f"Processing: {pdb_id} - {description}")
        print('='*70)
        
        try:
            # Download structure
            raw_pdb = download_pdb(pdb_id, output_dir)
            
            # Clean structure
            clean_pdb = f"{output_dir}/{pdb_id}_clean.pdb"
            clean_structure(raw_pdb, clean_pdb)
            
            # Show info
            get_structure_info(clean_pdb)
            
            print(f"\n✅ {pdb_id} processing complete!")
            
        except Exception as e:
            print(f"\n❌ Error processing {pdb_id}: {e}")
            continue
    
    print("\n" + "=" * 70)
    print("🎉 ALL STRUCTURES DOWNLOADED AND CLEANED!")
    print("=" * 70)
    print(f"\n📁 Files saved in: {output_dir}/")
    print("\n📋 Next steps:")
    print("   1. Identify binding pockets (fpocket or manual)")
    print("   2. Prepare receptors for docking (add hydrogens, charges)")
    print("   3. Define docking grid boxes")
    print("   4. Prepare ligand libraries for screening")
    

if __name__ == "__main__":
    main()