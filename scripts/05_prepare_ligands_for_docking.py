#!/usr/bin/env python3
"""
Prepare ligands for molecular docking
Converts SMILES to 3D structures and PDBQT format

Author: Olivier Nsekuye
Lab: GIGA-VIN, University of Liège
Date: January 2025

Pipeline:
1. SMILES → 3D structure (RDKit)
2. Energy minimization (MMFF94)
3. Add hydrogens (pH 7.4)
4. Convert to PDBQT (Meeko)
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import os
from pathlib import Path


def smiles_to_3d(smiles, optimize=True):
    """
    Convert SMILES to 3D structure.
    
    Args:
        smiles (str): SMILES string
        optimize (bool): Run energy minimization
        
    Returns:
        RDKit Mol object or None
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    
    if result != 0:
        # Embedding failed, try with different parameters
        result = AllChem.EmbedMolecule(
            mol, 
            randomSeed=42,
            useRandomCoords=True
        )
        
        if result != 0:
            return None
    
    # Energy minimization
    if optimize:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except:
            pass  # Continue even if optimization fails
    
    return mol


def prepare_ligand(mol, output_pdbqt):
    """
    Prepare ligand for docking (convert to PDBQT).
    
    Args:
        mol: RDKit Mol object
        output_pdbqt (str): Output PDBQT file path
        
    Returns:
        bool: True if successful
    """
    try:
        # Prepare molecule with Meeko
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        
        # Write PDBQT using updated API
        pdbqt_string = preparator.write_pdbqt_string()
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
        
        return True
        
    except Exception as e:
        print(f"   ⚠️  Meeko preparation failed: {e}")
        return False        


def process_library(input_file, output_dir, max_compounds=10):
    """
    Process compound library for docking.
    
    Args:
        input_file (str): Filtered compounds TSV
        output_dir (str): Output directory for PDBQT files
        max_compounds (int): Maximum compounds to process (for testing)
    """
    print("\n" + "=" * 70)
    print("LIGAND PREPARATION FOR DOCKING")
    print("=" * 70)
    print(f"\n📥 Input: {input_file}")
    print(f"📁 Output: {output_dir}")
    print(f"🔢 Processing: {max_compounds} compounds (test mode)")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Read filtered library
    if not Path(input_file).exists():
        print(f"\n❌ Input file not found: {input_file}")
        print(f"   Run 04_filter_ligands.py first")
        return
    
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return
    
    print(f"\n📊 Total compounds in library: {len(df):,}")
    
    # Process compounds
    print(f"\n⚙️  Processing compounds...")
    
    success_count = 0
    failed_3d = 0
    failed_pdbqt = 0
    
    # Process subset
    n_process = min(max_compounds, len(df))
    
    for idx, row in df.head(n_process).iterrows():
        zinc_id = row.get('ZINC_ID', f'compound_{idx}')
        smiles = row.get('SMILES', '')
        
        print(f"   [{idx+1}/{n_process}] {zinc_id}...", end='')
        
        if not smiles:
            print(" ❌ No SMILES")
            continue
        
        # Convert to 3D
        mol = smiles_to_3d(smiles)
        
        if mol is None:
            print(" ❌ 3D generation failed")
            failed_3d += 1
            continue
        
        # Prepare for docking
        output_pdbqt = f"{output_dir}/{zinc_id}.pdbqt"
        
        if prepare_ligand(mol, output_pdbqt):
            print(" ✅")
            success_count += 1
        else:
            print(" ❌ PDBQT conversion failed")
            failed_pdbqt += 1
    
    # Statistics
    print("\n" + "=" * 70)
    print("PREPARATION STATISTICS")
    print("=" * 70)
    print(f"\n✅ Successfully prepared: {success_count}/{n_process}")
    print(f"❌ 3D generation failed: {failed_3d}")
    print(f"❌ PDBQT conversion failed: {failed_pdbqt}")
    print(f"📈 Success rate: {success_count/n_process*100:.1f}%")
    
    if success_count > 0:
        print(f"\n📁 PDBQT files saved in: {output_dir}/")
        print(f"   Files ready for docking!")


def main():
    """Main execution"""
    
    print("=" * 70)
    print("🧬 RTC-PPI INHIBITOR DISCOVERY")
    print("   Ligand Preparation Pipeline")
    print("=" * 70)
    
    # Create demo compounds if filtered file doesn't exist
    input_file = 'data/libraries/zinc15/clean/filtered_compounds.tsv'
    
    if not Path(input_file).exists():
        print(f"\n⚠️  Filtered library not found")
        print(f"   Creating demo compounds for testing...")
        
        # Create demo library with PPI-suitable compounds
        demo_data = {
            'ZINC_ID': ['ZINC_DEMO_001', 'ZINC_DEMO_002', 'ZINC_DEMO_003'],
            'SMILES': [
                # Larger drug-like molecules suitable for PPIs
                'CC(C)Cc1ccc(cc1)C(C)C(=O)NCCC(=O)O',  # ~277 Da
                'COc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N',  # ~306 Da
                'Cc1ccc(cc1)S(=O)(=O)Nc2ccc(cc2C)C(=O)O'  # ~305 Da
            ],
            'MW': [277.36, 306.32, 305.35],
            'LogP': [2.5, 2.1, 2.8]
        }
        
        demo_df = pd.DataFrame(demo_data)
        Path(input_file).parent.mkdir(parents=True, exist_ok=True)
        demo_df.to_csv(input_file, sep='\t', index=False)
        
        print(f"   ✅ Demo library created: {input_file}")
    
    # Process library
    process_library(
        input_file=input_file,
        output_dir='data/docking/ligands_pdbqt',
        max_compounds=10  # Test with 10 compounds
    )
    
    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print("\n1. For production (HPC cluster):")
    print("   - Process full filtered library (millions of compounds)")
    print("   - Use parallel processing (100+ cores)")
    print("   - Expected time: 1-2 weeks")
    
    print("\n2. Quality control:")
    print("   - Visualize prepared structures")
    print("   - Check PDBQT format validity")
    print("   - Verify 3D geometry")
    
    print("\n3. Ready for docking:")
    print("   - Receptors: data/docking/receptors/*.pdbqt")
    print("   - Ligands: data/docking/ligands_pdbqt/*.pdbqt")
    print("   - Deploy VirtualFlow on CECI HPC")
    
    print("\n💡 TIP: Test docking workflow locally with 10-100 compounds")
    print("   before scaling to billions on HPC")


if __name__ == "__main__":
    main()
