# 7DFG Discovery Results - COMPLETE

**Date:** February 7, 2026  
**Status:** ✅ Discovery completed  
**Notebook:** `05_7DFG_interface_discovery.ipynb`

---

## Summary Statistics

| Interface | Atomic Contacts | Residue Pairs | Hot Spots Found | Strongest |
|-----------|----------------|---------------|-----------------|-----------|
| NSP12-NSP7 | 383 | 55 | **1** ⚠️ | GLU431-LYS2 (3.29 Å) |
| NSP12-NSP8 | 1,176 | 169 | **2** | ARG331-ASP112 (3.02 Å) |

---

## Interface 1: NSP12-NSP7 ⚠️ SPARSE

### Hot Spots Discovered: 1 (ONLY ONE!)

**Primary (and ONLY) Hot Spot:**
- **Residues:** GLU431 (NSP12) ←→ LYS2 (NSP7)
- **Distance:** 3.29 Å
- **Atoms:** OE2 (GLU) - CB (LYS)
- **Type:** Salt bridge (strong, < 4 Å)
- **Rank:** #1 of 1 (by default!)

### Analysis
- ✅ 55 unique residue pairs within 5 Å
- ⚠️ Only **1** charged interaction found
- ⚠️ **NO secondary hot spots!**
- 📊 This is unusual - suggests:
  1. Very specific interaction
  2. Or missed interactions (need comprehensive check)
  3. Interface may not be primarily electrostatic

### Critical Questions
1. **Why only 1?** 
   - Is this really the only charged pair < 5 Å?
   - Could there be K-E or R-D pairs we missed?
2. **What about 5-7 Å range?**
   - Weak salt bridges?
   - H-bonds?
3. **Is GLU431-LYS2 validated?**
   - Need ALL pairwise distances to confirm

---

## Interface 2: NSP12-NSP8 ✓ BETTER DEFINED

### Hot Spots Discovered: 2

**Primary Hot Spot:**
- **Residues:** ARG331 (NSP12) ←→ ASP112 (NSP8)
- **Distance:** 3.02 Å ⭐ (VERY STRONG)
- **Atoms:** NH1 (ARG) - O (ASP)
- **Type:** Salt bridge (strong, < 4 Å)
- **Rank:** #1 of 2

**Secondary Hot Spot:**
- **Residues:** ASP523 (NSP12) ←→ ARG80 (NSP8)
- **Distance:** 3.74 Å
- **Atoms:** OD2 (ASP) - NH2 (ARG)
- **Type:** Salt bridge (strong, < 4 Å)
- **Rank:** #2 of 2

### Analysis
- ✅ 169 unique residue pairs within 5 Å
- ✅ 2 charged interactions found
- ✅ Both are strong salt bridges (< 4 Å)
- ⚠️ Still only 2 - relatively sparse
- 📊 Gap between #1 and #2: 0.72 Å

### Critical Questions
1. **Why only 2?**
   - Expected more for 169 residue pairs
   - Need to check ALL charged pairs
2. **Is there a charged triad?**
   - Could ARG331 interact with other ASP?
   - Could ASP523 interact with other ARG?
3. **Which NSP8 chain?**
   - Chain B or G?
   - Do they have different hot spots?

---

## Comparison with 6W4H

| Metric | 6W4H (NSP10-NSP16) | 7DFG NSP12-NSP7 | 7DFG NSP12-NSP8 |
|--------|-------------------|-----------------|-----------------|
| **Hot Spots Found** | 3 | **1** ⚠️ | 2 |
| **Primary Distance** | 5.15 Å | **3.29 Å** ✓ | **3.02 Å** ⭐ |
| **Primary Type** | Salt bridge | Salt bridge | Salt bridge |
| **Charged Triad?** | Yes (K76-K78-D107) | ❓ Unknown | ❓ Unknown |
| **Interface Pairs** | 144 total | 55 pairs | 169 pairs |

**Key Observations:**
1. ✅ 7DFG hot spots are **STRONGER** (shorter distances)
2. ⚠️ 7DFG has **FEWER** hot spots (more selective?)
3. ❓ Need comprehensive analysis to validate

---

## Red Flags & Concerns

### 🚨 NSP12-NSP7 Interface
- **ONLY 1 hot spot found** - This is suspicious!
- Possible issues:
  1. Discovery cutoff too strict (< 5 Å)?
  2. Missed K-E and R-D combinations?
  3. Interface not primarily electrostatic?

### ⚠️ NSP12-NSP8 Interface  
- Only 2 hot spots for 169 pairs seems low
- Need to verify these are truly the only ones

### 🔍 Chain Confusion
- NSP8 has TWO chains (B and G)
- Discovery notebook may have only analyzed one
- Need to check BOTH chains B AND G separately

---

## Immediate Actions Required

### 1. Verify Discovery Results ⚠️ URGENT
Before comprehensive analysis, check:
- [ ] Did Section 4A screen ALL charge combinations?
  - LYS-ASP ✓
  - LYS-GLU ❓
  - ARG-ASP ❓
  - ARG-GLU ❓
- [ ] Did Section 4B check BOTH NSP8 chains (B and G)?

### 2. Expand Search Criteria
- [ ] Check 5-7 Å range (weak salt bridges)
- [ ] Check 7-10 Å range (H-bonds)
- [ ] Look for ALL charged combinations

### 3. Run Comprehensive Validation
Create `09_7DFG_comprehensive_validation.ipynb` to:
- [ ] Find ALL LYS + ARG in NSP12
- [ ] Find ALL ASP + GLU in NSP7
- [ ] Find ALL ASP + GLU in NSP8 (both B and G)
- [ ] Calculate complete distance matrix
- [ ] Rank ALL interactions
- [ ] Validate discovered hot spots

---

## Hypothesis: Sparse vs Dense Interfaces

**6W4H (Dense):**
- 3 hot spots in small interface
- Multiple lysines clustering around D107
- Charged triad formation

**7DFG (Sparse):**
- 1-2 hot spots per interface
- Very specific interactions
- No obvious clustering (yet)

**Implications:**
- Sparse interfaces may be harder to drug (fewer targets)
- Or could be easier (highly specific single target)
- Need comprehensive analysis to decide!

---

## Next Steps

### Priority 1: Verify Discovery (Today)
Check notebook code to ensure ALL combinations screened

### Priority 2: Comprehensive Analysis (Tomorrow)
- Create exhaustive pairwise distance matrix
- Include ALL charge types (K, R, D, E)
- Check 0-10 Å range (not just < 5 Å)

### Priority 3: Strategic Decision
Based on comprehensive results:
- If truly sparse → Focus on those 1-2 hot spots
- If hidden hot spots found → Target clusters
- Compare with 6W4H to inform strategy

