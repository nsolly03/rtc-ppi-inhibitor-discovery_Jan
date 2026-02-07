#!/usr/bin/env python3
"""
validate_structure.py

Validates a local PDB structure against expected NSP content
by analyzing chain lengths and basic composition.

Usage:
    python3 scripts/validate_structure.py <PDB_ID> <PDB_FILE>
"""

import sys
from Bio.PDB import PDBParser


# Expected NSP sizes (approximate, SARS-CoV-2)
EXPECTED_SIZES = {
    "NSP12": (880, 940),   # RdRp
    "NSP8":  (90, 130),    # Cofactor
    "NSP7":  (60, 80),     # Cofactor
    "NSP13": (580, 620),   # Helicase
    "NSP10": (120, 150),   # Cofactor
    "NSP14": (500, 540),   # Exonuclease
}


def analyze_local_structure(pdb_file):
    """
    Parse PDB file and return chain length information.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("STRUCTURE", pdb_file)

    # ✅ Robust model selection (NO assumption about model ID)
    model = next(structure.get_models())

    chain_info = {}

    for chain in model:
        residues = [
            res for res in chain
            if res.id[0] == " "  # exclude hetero/water
        ]

        if not residues:
            continue

        pdb_numbers = [res.id[1] for res in residues]

        chain_info[chain.id] = {
            "length": len(residues),
            "pdb_start": min(pdb_numbers),
            "pdb_end": max(pdb_numbers),
        }

    return chain_info


def assign_nsp_by_size(chain_info):
    """
    Assign NSP identities based on chain length heuristics.
    """
    assignments = {}
    used_chains = set()

    for nsp, (min_len, max_len) in EXPECTED_SIZES.items():
        for chain_id, info in chain_info.items():
            if chain_id in used_chains:
                continue

            if min_len <= info["length"] <= max_len:
                assignments[nsp] = chain_id
                used_chains.add(chain_id)
                break

    return assignments


def validate_structure(pdb_id, pdb_file):
    print("=" * 70)
    print(f"STRUCTURE INFORMATION - {pdb_id}")
    print("=" * 70)
    print()

    chain_info = analyze_local_structure(pdb_file)

    for chain_id, info in chain_info.items():
        print(
            f"Chain {chain_id}: {info['length']} residues "
            f"(PDB {info['pdb_start']}-{info['pdb_end']})"
        )

    print()
    print("Analyzing chain sizes...")
    print()

    assignments = assign_nsp_by_size(chain_info)

    for chain_id, info in chain_info.items():
        assigned = None
        for nsp, cid in assignments.items():
            if cid == chain_id:
                assigned = nsp
                break

        if assigned:
            print(f"Chain {chain_id} ({info['length']} res) → {assigned}")
        else:
            if info["length"] < 50:
                print(f"Chain {chain_id} ({info['length']} res) → RNA/Substrate/Ligand")
            else:
                print(f"Chain {chain_id} ({info['length']} res) → Unknown")

    print()
    print("=" * 70)
    print("FINAL ASSIGNMENTS:")
    print("=" * 70)

    for key in ["NSP12", "NSP7", "NSP8", "NSP13", "NSP10", "NSP14"]:
        if key in assignments:
            print(f"  {key} = Chain {assignments[key]}")
        else:
            print(f"  {key} = NOT FOUND ❌")

    print()
    print("Target Interfaces for Discovery:")

    if "NSP12" in assignments and "NSP7" in assignments:
        print("  1. NSP12-NSP7 interface")
    if "NSP12" in assignments and "NSP8" in assignments:
        print("  2. NSP12-NSP8 interface")
    if "NSP13" in assignments:
        print("  3. NSP13 helicase interfaces")

    print("=" * 70)
    print()


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 validate_structure.py <PDB_ID> <PDB_FILE>")
        sys.exit(1)

    pdb_id = sys.argv[1]
    pdb_file = sys.argv[2]

    validate_structure(pdb_id, pdb_file)


if __name__ == "__main__":
    main()

