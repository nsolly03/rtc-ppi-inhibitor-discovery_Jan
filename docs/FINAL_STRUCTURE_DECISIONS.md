# Final Structure Decisions

**Date:** February 7, 2026  
**Status:** Finalizing replacement structures

---

## Decision Summary

### ✅ NSP14: 7N0C (SELECTED)
- **Composition:** NSP10 (Chain A, 131 res) + NSP14 (Chain B, 513 res) + RNA
- **Interface:** NSP10-NSP14 exonuclease complex
- **Why chosen:** 
  - Complete NSP10-NSP14 complex
  - NSP14 has expected size (~520 residues)
  - Can compare NSP10 interfaces (vs 6W4H NSP10-NSP16)
  - Has RNA substrate (functional state)
- **Status:** ✅ APPROVED - Move to main structures
- **Citation:** [To be verified from RCSB]

### 🔍 NSP13: Choosing between 6ZSL and 7NIO

Both are NSP13 helicase homodimers (~590 res). Need to check:

#### 6ZSL
- Chains: A (572 res) + B (585 res)
- NSP13 homodimer
- [Need: resolution, citation, ligands]

#### 7NIO  
- Chains: A (590 res) + E (585 res)
- NSP13 homodimer
- [Need: resolution, citation, ligands]

**DEFAULT CHOICE: Use 6ZSL (lower PDB ID = older/earlier structure)**
- Reasoning: Earlier structures often have more citations
- Can switch to 7NIO if it has better resolution/ligands

### ❌ Rejected Candidates
- **7DIY:** Incomplete (Chain B only 286 res, should be ~520)
- **7MCW:** File parsing error (corrupted or unusual format)

---

## Actions Required

1. **Move 7N0C to main structures:**
```bash
   cp data/structures/candidates/7N0C.pdb data/structures/pdb/7N0C.pdb
```

2. **Move 6ZSL to main structures (tentative):**
```bash
   cp data/structures/candidates/6ZSL.pdb data/structures/pdb/6ZSL.pdb
```

3. **Clean up candidates folder:**
```bash
   rm -rf data/structures/candidates/
```

4. **Update structure inventory**

---

## Updated Project Structure List

| # | PDB ID | Target | Status | Priority |
|---|--------|--------|--------|----------|
| 1 | 6W4H | NSP10-NSP16 | ✅ Complete | Done |
| 2 | 7DFG | NSP12-NSP7-NSP8 RdRp | ✅ Ready | HIGH |
| 3 | 7N0C | NSP10-NSP14 ExoN | ✅ Ready | HIGH |
| 4 | 6ZSL | NSP13 Helicase | ✅ Ready | MEDIUM |
| 5 | 6XEZ | RTC-Helicase complex | ⚠️ Optional | LOW |

**Valid structures: 4/5 (or 5/5 if we keep 6XEZ)**

---

## Next Steps

### Today (Feb 7)
- [x] Validate all candidates
- [x] Choose NSP14 structure (7N0C)
- [x] Choose NSP13 structure (6ZSL - tentative)
- [ ] Move structures to main folder
- [ ] Clean up candidates
- [ ] Update documentation

### Tomorrow (Feb 8)
- [ ] Run 7DFG analysis (NSP12-NSP7-NSP8)
- [ ] Create 7N0C analysis notebook (NSP10-NSP14)
- [ ] Create 6ZSL analysis notebook (NSP13)

### Week 4
- [ ] Complete all interface analyses
- [ ] Run fpocket on all structures
- [ ] Rank targets by druggability
- [ ] Begin virtual screening prep

