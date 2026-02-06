#!/usr/bin/env python3
"""
Download all 5 coronavirus RTC structures from RCSB PDB
"""

import requests
from pathlib import Path
from Bio import PDB
import sys

# Structure IDs to download
STRUCTURES = {
    '6W4H': 'NSP10-NSP16 (2\'-O-MTase)',
    '7DFG': 'NSP12-NSP7-NSP8 (RdRp)',
    '6XEZ': 'NSP12-NSP7-NSP8 (RdRp alt)',
    '7EDI': 'NSP10-NSP14 (ExoN)',
    '6W9C': 'NSP13 (Helicase)'
}

# Output directory
PDB_DIR = Path('data/structures/pdb')
PDB_DIR.mkdir(exist_ok=True, parents=True)

def download_pdb(pdb_id):
    """Download PDB file from RCSB"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_file = PDB_DIR / f"{pdb_id}.pdb"
    
    print(f"Downloading {pdb_id}...", end=' ')
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        with open(output_file, 'w') as f:
            f.write(response.text)
        
        # Verify it's a valid PDB
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(output_file))
        num_chains = len(list(structure[0].get_chains()))
        
        print(f"✓ Downloaded ({num_chains} chains)")
        return True
        
    except Exception as e:
        print(f"✗ Failed: {e}")
        return False

def main():
    print("="*70)
    print("DOWNLOADING CORONAVIRUS RTC STRUCTURES")
    print("="*70)
    print()
    
    success_count = 0
    
    for pdb_id, description in STRUCTURES.items():
        print(f"{pdb_id}: {description}")
        
        output_file = PDB_DIR / f"{pdb_id}.pdb"
        
        if output_file.exists():
            print(f"  Already exists, skipping...")
            success_count += 1
        else:
            if download_pdb(pdb_id):
                success_count += 1
        print()
    
    print("="*70)
    print(f"DOWNLOAD COMPLETE: {success_count}/{len(STRUCTURES)} structures")
    print("="*70)
    
    if success_count == len(STRUCTURES):
        print("\n✓ All structures ready for analysis!")
        return 0
    else:
        print(f"\n⚠ {len(STRUCTURES) - success_count} structures failed")
        return 1

if __name__ == '__main__':
    sys.exit(main())
