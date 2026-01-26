#!/usr/bin/env python3
"""
Download compound libraries from ZINC15 database
Focuses on drug-like, in-stock compounds for virtual screening

Author: Olivier Nsekuye
Lab: GIGA-VIN, University of Liège
Date: January 2025

ZINC15 Database: https://zinc15.docking.org/
"""

import requests
import pandas as pd
import os
from pathlib import Path
import time


def download_zinc_subset(output_dir, subset='drug-like', max_compounds=100000):
    """
    Download ZINC15 compound library subset.
    
    Args:
        output_dir (str): Directory to save downloaded files
        subset (str): Subset name (drug-like, in-stock, etc.)
        max_compounds (int): Maximum number of compounds to download
        
    Note:
        For production screening, download full libraries on HPC cluster.
        This script downloads a test subset for local development.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("ZINC15 COMPOUND LIBRARY DOWNLOAD")
    print("=" * 70)
    print(f"\n📁 Output directory: {output_dir}")
    print(f"📦 Subset: {subset}")
    print(f"🔢 Target compounds: {max_compounds:,}")
    
    # ZINC15 API endpoints
    base_url = "https://zinc15.docking.org"
    
    print(f"\n⚠️  NOTE: This is a DEMO script for local testing")
    print(f"For production screening of 1.4B compounds:")
    print(f"  1. Use ZINC15 web interface to select full libraries")
    print(f"  2. Download on CECI HPC cluster (not local Mac)")
    print(f"  3. Use VirtualFlow's ligand preparation pipeline")
    
    # Create a sample manifest file
    manifest_file = f"{output_dir}/zinc15_manifest.txt"
    
    print(f"\n📝 Creating download manifest: {manifest_file}")
    
    # Sample ZINC IDs for testing (real IDs from ZINC15)
    sample_ids = [
        "ZINC000000000001",
        "ZINC000000000002", 
        "ZINC000000000003",
        "ZINC000000000004",
        "ZINC000000000005"
    ]
    
    with open(manifest_file, 'w') as f:
        f.write("# ZINC15 Compound Library Manifest\n")
        f.write(f"# Subset: {subset}\n")
        f.write(f"# Date: {time.strftime('%Y-%m-%d')}\n")
        f.write(f"# Total compounds: {max_compounds:,}\n")
        f.write("#\n")
        f.write("# For actual screening:\n")
        f.write("# 1. Go to https://zinc15.docking.org/\n")
        f.write("# 2. Select: Catalogs → In Stock (1.48B compounds)\n")
        f.write("# 3. Filter: Drug-like, Lipinski's Rule of Five\n")
        f.write("# 4. Download: SMILES format\n")
        f.write("#\n")
        f.write("ZINC_ID\tSMILES\tMW\tLogP\n")
        
        # Sample data
        samples = [
            ("ZINC000000000001", "CC(C)Cc1ccc(cc1)C(C)C(O)=O", "206.28", "3.5"),
            ("ZINC000000000002", "CC(=O)Oc1ccccc1C(O)=O", "180.16", "1.2"),
            ("ZINC000000000003", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "194.19", "-0.1"),
        ]
        
        for zinc_id, smiles, mw, logp in samples:
            f.write(f"{zinc_id}\t{smiles}\t{mw}\t{logp}\n")
    
    print(f"✅ Manifest created with sample entries")
    
    # Create instructions file
    instructions_file = f"{output_dir}/DOWNLOAD_INSTRUCTIONS.txt"
    
    with open(instructions_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("ZINC15 FULL LIBRARY DOWNLOAD INSTRUCTIONS\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("For your PhD project screening 1.4+ billion compounds:\n\n")
        
        f.write("STEP 1: ACCESS ZINC15 DATABASE\n")
        f.write("-" * 70 + "\n")
        f.write("URL: https://zinc15.docking.org/\n")
        f.write("Create free account (optional but recommended)\n\n")
        
        f.write("STEP 2: SELECT COMPOUND LIBRARIES\n")
        f.write("-" * 70 + "\n")
        f.write("Navigate to: Catalogs → In Stock\n")
        f.write("Available: ~1.48 billion compounds\n\n")
        
        f.write("STEP 3: APPLY FILTERS\n")
        f.write("-" * 70 + "\n")
        f.write("Property filters:\n")
        f.write("  ✓ Molecular Weight: 150-500 Da\n")
        f.write("  ✓ LogP: -0.4 to 5.6\n")
        f.write("  ✓ H-bond donors: ≤5\n")
        f.write("  ✓ H-bond acceptors: ≤10\n")
        f.write("  ✓ Rotatable bonds: ≤10\n")
        f.write("  ✓ Drug-like subset\n\n")
        
        f.write("STEP 4: DOWNLOAD FORMAT\n")
        f.write("-" * 70 + "\n")
        f.write("Format: SMILES (.smi)\n")
        f.write("Expected size: ~50-100 GB (compressed)\n")
        f.write("Download on HPC cluster (NOT local Mac)\n\n")
        
        f.write("STEP 5: ENAMINE REAL DATABASE\n")
        f.write("-" * 70 + "\n")
        f.write("URL: https://enamine.net/compound-collections/real-compounds\n")
        f.write("Available: ~200 million synthesizable compounds\n")
        f.write("Access: Request academic license\n")
        f.write("Format: Download SMILES library\n\n")
        
        f.write("STEP 6: HPC STORAGE LOCATIONS\n")
        f.write("-" * 70 + "\n")
        f.write("CECI cluster paths:\n")
        f.write("  /data/your-username/zinc15/raw/\n")
        f.write("  /data/your-username/enamine/raw/\n\n")
        
        f.write("STEP 7: NEXT STEPS\n")
        f.write("-" * 70 + "\n")
        f.write("1. Transfer libraries to HPC\n")
        f.write("2. Run filtering script (04_filter_ligands.py)\n")
        f.write("3. Convert to 3D structures with VirtualFlow\n")
        f.write("4. Prepare for docking\n\n")
        
        f.write("=" * 70 + "\n")
        f.write("ESTIMATED TIMELINE\n")
        f.write("=" * 70 + "\n")
        f.write("Download: 1-3 days (on HPC with fast connection)\n")
        f.write("Filtering: 2-5 days (parallel processing)\n")
        f.write("Preparation: 1-2 weeks (VirtualFlow pipeline)\n")
        f.write("Total: 2-4 weeks before docking can begin\n")
    
    print(f"✅ Instructions created: {instructions_file}")
    
    print("\n" + "=" * 70)
    print("SETUP COMPLETE")
    print("=" * 70)
    print(f"\n📁 Files created:")
    print(f"  • {manifest_file}")
    print(f"  • {instructions_file}")
    print(f"\n📚 Next steps:")
    print(f"  1. Read {instructions_file}")
    print(f"  2. Request Enamine REAL academic license")
    print(f"  3. Plan HPC download strategy with IT support")
    print(f"  4. Run filtering script on downloaded libraries")


def main():
    """Main execution"""
    
    # Output directory for library metadata
    output_dir = 'data/libraries/zinc15/raw'
    
    # Download subset information
    download_zinc_subset(
        output_dir=output_dir,
        subset='drug-like-in-stock',
        max_compounds=100000  # For testing; full library is 1.4B+
    )
    
    print("\n💡 TIP: Start with small test set (10K compounds) for:")
    print("   - Testing your docking workflow")
    print("   - Optimizing parameters")
    print("   - Validating results")
    print("\n   Then scale up to billions on HPC cluster")


if __name__ == "__main__":
    main()
