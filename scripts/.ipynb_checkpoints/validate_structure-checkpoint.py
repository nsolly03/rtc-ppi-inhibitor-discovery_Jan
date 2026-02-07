#!/usr/bin/env python3
"""
Validate PDB structure against official PDB metadata
"""

import requests
from Bio import PDB
import sys
from pathlib import Path

def fetch_pdb_metadata(pdb_id):
    """Fetch official metadata from RCSB PDB"""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"⚠️  Could not fetch metadata for {pdb_id}")
            return None
    except Exception as e:
        print(f"⚠️  Error fetching metadata: {e}")
        return None

def fetch_polymer_entities(pdb_id):
    """Fetch polymer entity information"""
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    entities = []
    entity_num = 1
    
    while True:
        try:
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_num}"
            response = requests.get(entity_url)
            if response.status_code == 200:
                entities.append(response.json())
                entity_num += 1
            else:
                break
        except:
            break
    
    return entities

def analyze_local_structure(pdb_file):
    """Analyze local PDB file"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('struct', str(pdb_file))
    model = structure[0]
    
    chains_info = []
    for chain in model:
        residues = [r for r in chain.get_residues() if r.id[0] == ' ']
        if residues:
            chains_info.append({
                'chain_id': chain.id,
                'num_residues': len(residues),
                'first_res': residues[0].id[1],
                'last_res': residues[-1].id[1]
            })
    
    return chains_info

def validate_structure(pdb_id, pdb_file):
    """Complete validation"""
    print("="*70)
    print(f"VALIDATING STRUCTURE: {pdb_id}")
    print("="*70)
    print()
    
    # 1. Analyze local file
    print("📂 LOCAL FILE ANALYSIS:")
    local_chains = analyze_local_structure(pdb_file)
    for chain in local_chains:
        print(f"  Chain {chain['chain_id']}: {chain['num_residues']} residues (PDB {chain['first_res']}-{chain['last_res']})")
    print()
    
    # 2. Fetch official metadata
    print("🌐 FETCHING OFFICIAL PDB METADATA...")
    metadata = fetch_pdb_metadata(pdb_id)
    
    if metadata:
        print(f"✓ Title: {metadata.get('struct', {}).get('title', 'N/A')}")
        print(f"✓ Release Date: {metadata.get('rcsb_accession_info', {}).get('initial_release_date', 'N/A')}")
        print()
        
        # Get citation info
        if 'citation' in metadata:
            citations = metadata['citation']
            if citations:
                first_citation = citations[0]
                print("📄 PRIMARY CITATION:")
                print(f"  Title: {first_citation.get('title', 'N/A')}")
                print(f"  Authors: {first_citation.get('rcsb_authors', ['N/A'])[0]} et al.")
                print(f"  Journal: {first_citation.get('journal_abbrev', 'N/A')} ({first_citation.get('year', 'N/A')})")
                print()
    
    # 3. Fetch polymer entities
    print("🧬 OFFICIAL POLYMER ENTITIES:")
    entities = fetch_polymer_entities(pdb_id)
    
    if entities:
        for entity in entities:
            entity_desc = entity.get('rcsb_polymer_entity', {}).get('pdbx_description', 'N/A')
            entity_chains = entity.get('entity_poly', {}).get('pdbx_strand_id', 'N/A')
            entity_length = entity.get('entity_poly', {}).get('rcsb_sample_sequence_length', 0)
            
            print(f"  Entity: {entity_desc}")
            print(f"    Chains: {entity_chains}")
            print(f"    Length: {entity_length} residues")
            print()
    else:
        print("  ⚠️  Could not fetch entity information")
        print()
    
    # 4. Cross-validation
    print("="*70)
    print("✅ VALIDATION SUMMARY:")
    print("="*70)
    print(f"✓ Local file has {len(local_chains)} chains")
    if metadata:
        print(f"✓ Metadata retrieved successfully")
    else:
        print(f"⚠️  Metadata retrieval failed - manual verification needed")
    
    if entities:
        print(f"✓ Found {len(entities)} polymer entities in PDB database")
        print()
        print("📋 RECOMMENDED ACTIONS:")
        print("  1. Compare local chain count with official entity count")
        print("  2. Match chain IDs with entity descriptions")
        print("  3. Verify residue counts match expected protein lengths")
        print("  4. Read the primary citation for experimental details")
    
    print()
    print("="*70)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python validate_structure.py <PDB_ID> <path_to_pdb_file>")
        sys.exit(1)
    
    pdb_id = sys.argv[1].upper()
    pdb_file = Path(sys.argv[2])
    
    if not pdb_file.exists():
        print(f"Error: File not found: {pdb_file}")
        sys.exit(1)
    
    validate_structure(pdb_id, pdb_file)
