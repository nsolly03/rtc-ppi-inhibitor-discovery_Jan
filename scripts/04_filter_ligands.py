#!/usr/bin/env python3
"""
Filter compound libraries based on drug-like properties
Applies Lipinski's Rule of Five and removes PAINS

Author: Olivier Nsekuye
Lab: GIGA-VIN, University of Liège
Date: January 2025

Filtering Criteria:
- Lipinski's Rule of Five (drug-likeness)
- PAINS removal (Pan-Assay Interference Compounds)
- Custom filters for PPI inhibitors
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, FilterCatalog
from pathlib import Path
import sys


def calculate_descriptors(smiles):
    """
    Calculate molecular descriptors.
    
    Args:
        smiles (str): SMILES string
        
    Returns:
        dict: Molecular descriptors or None if invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    return {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotBonds': Descriptors.NumRotatableBonds(mol),
        'Rings': Descriptors.RingCount(mol),
        'AromaticRings': Descriptors.NumAromaticRings(mol)
    }


def lipinski_rule_of_five(descriptors):
    """
    Check Lipinski's Rule of Five.
    
    Criteria:
    - MW ≤ 500 Da
    - LogP ≤ 5
    - H-bond donors ≤ 5
    - H-bond acceptors ≤ 10
    
    Args:
        descriptors (dict): Molecular descriptors
        
    Returns:
        bool: True if passes all criteria
    """
    if descriptors is None:
        return False
    
    return (
        descriptors['MW'] <= 500 and
        descriptors['LogP'] <= 5 and
        descriptors['HBD'] <= 5 and
        descriptors['HBA'] <= 10
    )


def ppi_optimized_filter(descriptors):
    """
    Optimized filters for PPI inhibitors.
    
    PPIs typically require:
    - Higher MW (350-700 Da) - larger binding surface
    - More rings (≥3) - rigidity
    - LogP 2-5 - moderate lipophilicity
    
    Args:
        descriptors (dict): Molecular descriptors
        
    Returns:
        bool: True if suitable for PPI targeting
    """
    if descriptors is None:
        return False
    
    return (
        350 <= descriptors['MW'] <= 700 and
        2 <= descriptors['LogP'] <= 5 and
        descriptors['Rings'] >= 3 and
        descriptors['RotBonds'] <= 10 and
        descriptors['TPSA'] <= 140
    )


def has_pains(smiles):
    """
    Check for PAINS (Pan-Assay Interference Compounds).
    
    PAINS are molecules that frequently show activity in assays
    but are actually assay artifacts (false positives).
    
    Args:
        smiles (str): SMILES string
        
    Returns:
        bool: True if contains PAINS substructure
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return True
    
    # Load PAINS filter catalog
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    
    # Check if molecule matches any PAINS patterns
    return catalog.HasMatch(mol)


def filter_library(input_file, output_file, use_ppi_filter=True):
    """
    Filter compound library.
    
    Args:
        input_file (str): Input file (TSV/CSV with SMILES)
        output_file (str): Output file for filtered compounds
        use_ppi_filter (bool): Use PPI-optimized filters vs standard Lipinski
    """
    print("\n" + "=" * 70)
    print("COMPOUND LIBRARY FILTERING")
    print("=" * 70)
    print(f"\n📥 Input: {input_file}")
    print(f"📤 Output: {output_file}")
    print(f"🔬 Filter type: {'PPI-optimized' if use_ppi_filter else 'Lipinski Ro5'}")
    
    # Read library
    print(f"\n📖 Reading library...")
    
    if not Path(input_file).exists():
        print(f"❌ Error: Input file not found: {input_file}")
        return
    
    try:
        df = pd.read_csv(input_file, sep='\t', comment='#')
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return
    
    print(f"   Total compounds: {len(df):,}")
    
    # Filter compounds
    print(f"\n🔍 Filtering compounds...")
    
    results = []
    passed_lipinski = 0
    passed_ppi = 0
    passed_pains = 0
    invalid_smiles = 0
    
    for idx, row in df.iterrows():
        smiles = row.get('SMILES', '')
        
        if not smiles:
            invalid_smiles += 1
            continue
        
        # Calculate descriptors
        descriptors = calculate_descriptors(smiles)
        
        if descriptors is None:
            invalid_smiles += 1
            continue
        
        # Check PAINS
        if has_pains(smiles):
            continue
        
        passed_pains += 1
        
        # Apply drug-like filter
        if use_ppi_filter:
            if not ppi_optimized_filter(descriptors):
                continue
            passed_ppi += 1
        else:
            if not lipinski_rule_of_five(descriptors):
                continue
            passed_lipinski += 1
        
        # Add to results
        result = row.to_dict()
        result.update(descriptors)
        results.append(result)
        
        # Progress indicator
        if (idx + 1) % 100 == 0:
            print(f"   Processed: {idx + 1:,} compounds", end='\r')
    
    print(f"   Processed: {len(df):,} compounds ✓")
    
    # Create output dataframe
    if results:
        output_df = pd.DataFrame(results)
        
        # Save filtered library
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n✅ Filtered library saved: {output_file}")
    else:
        print(f"\n⚠️  No compounds passed filters!")
        return
    
    # Statistics
    print("\n" + "=" * 70)
    print("FILTERING STATISTICS")
    print("=" * 70)
    print(f"\n📊 Input compounds: {len(df):,}")
    print(f"❌ Invalid SMILES: {invalid_smiles:,}")
    print(f"✅ Passed PAINS filter: {passed_pains:,}")
    
    if use_ppi_filter:
        print(f"✅ Passed PPI filter: {passed_ppi:,}")
        print(f"📈 Success rate: {passed_ppi/len(df)*100:.1f}%")
    else:
        print(f"✅ Passed Lipinski Ro5: {passed_lipinski:,}")
        print(f"📈 Success rate: {passed_lipinski/len(df)*100:.1f}%")
    
    print(f"\n📁 Output compounds: {len(results):,}")
    
    # Property distribution
    print("\n" + "=" * 70)
    print("PROPERTY DISTRIBUTION (Filtered Set)")
    print("=" * 70)
    
    print(f"\nMolecular Weight:")
    print(f"   Mean: {output_df['MW'].mean():.1f} Da")
    print(f"   Range: {output_df['MW'].min():.1f} - {output_df['MW'].max():.1f} Da")
    
    print(f"\nLogP:")
    print(f"   Mean: {output_df['LogP'].mean():.1f}")
    print(f"   Range: {output_df['LogP'].min():.1f} - {output_df['LogP'].max():.1f}")
    
    print(f"\nH-bond Donors:")
    print(f"   Mean: {output_df['HBD'].mean():.1f}")
    
    print(f"\nH-bond Acceptors:")
    print(f"   Mean: {output_df['HBA'].mean():.1f}")
    
    print(f"\nRotatable Bonds:")
    print(f"   Mean: {output_df['RotBonds'].mean():.1f}")


def main():
    """Main execution"""
    
    print("=" * 70)
    print("🧬 RTC-PPI INHIBITOR DISCOVERY")
    print("   Ligand Library Filtering")
    print("=" * 70)
    
    # Filter ZINC15 manifest (demo)
    filter_library(
        input_file='data/libraries/zinc15/raw/zinc15_manifest.txt',
        output_file='data/libraries/zinc15/clean/filtered_compounds.tsv',
        use_ppi_filter=True  # Use PPI-optimized filters
    )
    
    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print("\n1. For full library filtering:")
    print("   - Download complete ZINC15/Enamine libraries on HPC")
    print("   - Run this script with parallel processing")
    print("   - Expected output: ~50-100 million compounds")
    
    print("\n2. Convert filtered compounds:")
    print("   - Run: python scripts/05_prepare_ligands_for_docking.py")
    print("   - Converts SMILES → 3D structures → PDBQT format")
    
    print("\n3. Quality control:")
    print("   - Visualize sample compounds")
    print("   - Check property distributions")
    print("   - Validate conversions")


if __name__ == "__main__":
    main()
