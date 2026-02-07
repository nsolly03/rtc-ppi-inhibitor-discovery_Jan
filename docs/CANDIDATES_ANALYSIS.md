# Candidate Structures Analysis

**Date:** February 7, 2026  
**Status:** Reviewing candidates for NSP14 and NSP13

---

## NSP14 Candidates (3 structures)

### 7N0C ✅ BEST CANDIDATE
- **Chains:** A (NSP10, 131 res) + B (NSP14, 513 res) + RNA
- **Composition:** NSP10-NSP14 complex ✓
- **Status:** ✅ PERFECT - Has both NSP10 and NSP14
- **Interface:** NSP10-NSP14 (can compare with 6W4H NSP10-NSP16)
- **Recommendation:** **USE THIS ONE**

### 7DIY ⚠️ INCOMPLETE
- **Chains:** A (NSP10, 131 res) + B (Unknown, 286 res)
- **Issue:** Chain B only 286 residues (NSP14 should be ~520 res)
- **Status:** ⚠️ Likely partial structure or fragment
- **Recommendation:** Skip - incomplete

### 7MCW ❌ ERROR
- **Status:** ❌ File parsing error (StopIteration)
- **Issue:** Corrupted file or unusual format
- **Recommendation:** Skip - technical issue

**DECISION: Use 7N0C for NSP10-NSP14**

---

## NSP13 Candidates (2 structures)

### 6ZSL - NSP13 dimer/multimer
- **Chains:** A (572 res) + B (585 res)
- **Composition:** NSP13 helicase (likely homodimer)
- **Status:** ✅ Valid NSP13 structure
- **Note:** Two chains of similar size = homodimer
- **Interface:** NSP13-NSP13 (homo-oligomer)

### 7NIO - NSP13 dimer/multimer
- **Chains:** A (590 res) + E (585 res)
- **Composition:** NSP13 helicase (likely homodimer)
- **Status:** ✅ Valid NSP13 structure
- **Note:** Two chains of similar size = homodimer
- **Interface:** NSP13-NSP13 (homo-oligomer)

**DECISION: Need to check RCSB metadata to choose between 6ZSL and 7NIO**

---

## Next Steps

1. ✅ **7N0C** - Move to main structures (NSP10-NSP14)
2. 🔍 **6ZSL vs 7NIO** - Need to validate with RCSB to choose best
3. ❌ Delete 7DIY (incomplete) and 7MCW (corrupted)

