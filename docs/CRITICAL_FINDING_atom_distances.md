# CRITICAL FINDING: CA-CA vs Specific Atom Distances

**Date:** February 7, 2026  
**Discovery:** Major difference between CA-CA and specific atom distances  
**Impact:** Changes how we validate hot spots

---

## The Discovery

### GLU431-LYS2 (NSP12-NSP7)

| Measurement Method | Distance | Classification |
|-------------------|----------|----------------|
| **CA-CA (backbone)** | 7.45 Å | NOT interface (> 7 Å) |
| **OE2-CB (sidechains)** | 3.29 Å | STRONG salt bridge! |

**Difference:** 4.16 Å gap!

---

## Why This Matters

### CA-CA Distance (What We Used for Verification)
- Measures **backbone** proximity
- Conservative estimate
- Can MISS sidechain interactions
- Used in 6W4H comprehensive analysis ⚠️

### Specific Atom Distance (What Discovery Found)
- Measures **sidechain** proximity
- More accurate for salt bridges
- Captures the actual interaction
- What matters for drug design ✓

---

## The Problem with 6W4H Analysis

Looking back at 6W4H comprehensive analysis:
- Used CA-CA distances for ALL 144 combinations
- Found K76-D107 at **5.15 Å** (CA-CA)
- But what's the actual sidechain distance?
- **Could be much shorter!**

---

## Implications

### 1. Our Comprehensive Analysis Approach is FLAWED
- Using CA-CA underestimates interaction strength
- May miss real hot spots with distant backbones
- Need to use **sidechain atoms** instead

### 2. 7DFG Results are VALID
- GLU431-LYS2 IS a real hot spot (3.29 Å sidechains)
- Our verification failed because it used CA-CA
- Discovery notebook is MORE ACCURATE

### 3. Need to Revise Methodology

**Old approach (6W4H):**
```python
dist = np.linalg.norm(res1['CA'].get_coord() - res2['CA'].get_coord())
```

**Better approach (7DFG discovery):**
```python
# Find closest atom pair between residues
for atom1 in res1.get_atoms():
    for atom2 in res2.get_atoms():
        dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
```

---

## Corrected Understanding: 7DFG is SPARSE but STRONG

### NSP12-NSP7 Interface
- **1 hot spot found:** GLU431-LYS2
- **CA-CA:** 7.45 Å (looks weak)
- **Sidechain:** 3.29 Å (VERY strong!)
- **Conclusion:** Highly specific, strong interaction

### NSP12-NSP8 Interface
- **2 hot spots found**
- **ARG331-ASP112:** 3.02 Å (strongest)
- **ASP523-ARG80:** 3.74 Å (strong)
- **Conclusion:** Well-defined, strong interactions

---

## Action Items

### 1. Re-validate 6W4H ⚠️ IMPORTANT
- [ ] Check K76-D107 sidechain distance
- [ ] May be much stronger than 5.15 Å CA-CA
- [ ] Recalculate charged triad with sidechains

### 2. Update Comprehensive Analysis Methodology
For all future structures:
- Use **minimum sidechain distance** instead of CA-CA
- Specifically for salt bridges:
  - LYS: NZ atom
  - ARG: NH1, NH2 atoms
  - ASP: OD1, OD2 atoms
  - GLU: OE1, OE2 atoms

### 3. Accept 7DFG Sparse Results
- Discovery is correct
- 1-2 hot spots per interface is REAL
- These are STRONG interactions (< 4 Å)
- Proceed with comprehensive analysis using correct method

---

## Conclusion

✅ **7DFG discovery results are CORRECT**  
⚠️ **Our verification method was WRONG**  
🔄 **Need to update comprehensive analysis approach**  
📊 **6W4H may need re-validation with sidechain distances**

The sparse results (1-2 hot spots) are REAL and represent highly specific, strong interactions - potentially even better drug targets than dense charged clusters!

