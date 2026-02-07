# Structure Validation Protocol

Use this checklist BEFORE analyzing any structure to avoid mistakes.

---

## Pre-Analysis Validation Checklist

### Step 1: Official PDB Metadata
- [ ] Run validation script: `python3 scripts/validate_structure.py <PDB_ID> <file>`
- [ ] Read official title from RCSB PDB
- [ ] Check publication citation
- [ ] Review experimental method and resolution

### Step 2: Literature Review
- [ ] Find primary citation on PubMed
- [ ] Read abstract and methods
- [ ] Identify exact protein composition
- [ ] Note any special features (ligands, RNA, modifications)

### Step 3: Manual Inspection
- [ ] Open structure in PyMOL or Chimera
- [ ] Visually inspect all chains
- [ ] Check for missing residues
- [ ] Verify chain assignments

### Step 4: Cross-Reference
- [ ] Compare local chain count with PDB entity count
- [ ] Match chain IDs with protein names
- [ ] Verify residue counts match expected lengths
- [ ] Check for symmetric copies (homo-oligomers)

---

## Chain Assignment Validation

### Method 1: Residue Count (Automated)
```python
# NSP12 (RdRp): ~900 residues
# NSP13 (Helicase): ~600 residues
# NSP10: ~140 residues
# NSP16: ~300 residues
# NSP14: ~520 residues
# NSP8: ~180-200 residues
# NSP7: ~70 residues
```

### Method 2: Sequence Alignment (Gold Standard)
```bash
# Extract FASTA from PDB
grep "^SEQRES" structure.pdb

# BLAST against UniProt SARS-CoV-2 proteome
# Confirm protein identity
```

### Method 3: Literature (Most Reliable)
- Always check the paper FIRST
- Trust published assignments over automated detection
- Note any discrepancies

---

## Common Mistakes to Avoid

❌ **Assuming similar structures are identical**
- 7DFG ≠ 6XEZ (6XEZ has nsp13!)
- Always validate each structure independently

❌ **Relying only on residue counts**
- Multiple proteins can have similar sizes
- Use literature + sequence alignment

❌ **Ignoring small chains**
- RNA, ligands, ions are important
- May indicate functional state

❌ **Not reading the paper**
- Paper explains what's in the structure
- Describes experimental conditions
- Notes any special features

---

## Validation Record

### 6W4H ✅ VALIDATED
- **Source:** Trepte et al. (2024) Mol Syst Biol
- **Composition:** NSP10 (Chain B) + NSP16 (Chain A)
- **Status:** Confirmed via literature
- **Hot Spot:** K76-D107 (validated)

### 7DFG ⏳ NEEDS VALIDATION
- **Source:** [Citation needed]
- **Composition:** [To be confirmed]
- **Status:** Pending validation
- **Action:** Run validation script + literature review

### 6XEZ ⚠️ COMPLEX STRUCTURE
- **Source:** Chen et al. (2020) Cell
- **Composition:** NSP12 + NSP7 + NSP8 + 2×NSP13 + RNA
- **Status:** More complex than expected
- **Note:** Contains helicase (nsp13) - different from simple RdRp
- **Decision:** Skip for now OR analyze nsp13 interfaces separately

### 7EDI ⏳ NEEDS VALIDATION
- **Source:** [Citation needed]
- **Composition:** [To be confirmed]
- **Status:** Pending

### 6W9C ⏳ NEEDS VALIDATION
- **Source:** [Citation needed]
- **Composition:** [To be confirmed]
- **Status:** Pending

---

## When to Reject a Structure

🚫 Reject if:
- Composition doesn't match project goals
- Too complex for current analysis
- Poor resolution (>3.5 Å for drug design)
- Missing critical regions
- Contains mutations/variants

✅ Accept if:
- Composition matches target interface
- Good resolution (<3.0 Å ideal)
- Complete structure (no large missing loops)
- Wild-type sequence
- Well-documented in literature

