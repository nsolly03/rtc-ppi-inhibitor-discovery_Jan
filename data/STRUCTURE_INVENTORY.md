# Structure Inventory - FINAL

**Last Updated:** February 7, 2026, 12:00 CAT  
**Status:** ✅ All 5 structures validated and finalized

---

## Complete Structure List (5/5)

### 1. 6W4H - NSP10-NSP16 Methyltransferase ✅ COMPLETE
- **Location:** `data/structures/pdb/6W4H.pdb`
- **Composition:** NSP10 (Chain B, 116 res) + NSP16 (Chain A, 299 res)
- **Citation:** Rosas-Lemus et al. (2020) Sci. Signal.
- **Title:** 1.80 Angstrom Resolution Crystal Structure
- **Resolution:** 1.80 Å
- **Release:** 2020-03-18
- **Status:** ✅ ANALYZED - Hot spot K76-D107 validated
- **Notebook:** `02_analyze_6W4H_interface.ipynb`
- **Results:** `data/analysis_results/6W4H_*`
- **Interface:** 22 residues (16 NSP10, 6 NSP16)
- **Grid Box:** Center (75.883, 11.641, 10.087), Size 25×25×25 Ų

---

### 2. 7DFG - NSP12-NSP7-NSP8 RdRp Holoenzyme 🔄 READY
- **Location:** `data/structures/pdb/7DFG.pdb`
- **Composition:** NSP12 (Chain A, 908 res) + NSP7 (Chain C, 69 res) + NSP8 (Chains B,G, 117+106 res) + RNA (Chains P,T)
- **Citation:** Yin et al. (2021) To be published
- **Title:** Structure bound to favipiravir
- **Release:** 2021-11-17
- **Status:** 🔄 NOTEBOOK READY - Awaiting analysis
- **Notebook:** `05_7DFG_interface_discovery.ipynb` (complete, 10 sections)
- **Target Interfaces:** NSP12-NSP7 and NSP12-NSP8
- **Note:** Has 2 NSP8 copies (Chains B and G)
- **Priority:** HIGH - Run this next!

---

### 3. 7N0C - NSP10-NSP14 Exonuclease Complex 🆕 NEW
- **Location:** `data/structures/pdb/7N0C.pdb`
- **Composition:** NSP10 (Chain A, 131 res) + NSP14 (Chain B, 513 res) + RNA (Chains T,D)
- **Citation:** [To be determined from RCSB]
- **Status:** 🆕 VALIDATED - Needs notebook creation
- **Target Interface:** NSP10-NSP14
- **Comparison:** Can compare NSP10 binding with 6W4H (NSP10-NSP16)
- **Priority:** HIGH
- **Next Step:** Create `07_7N0C_interface_discovery.ipynb`

---

### 4. 6ZSL - NSP13 Helicase 🆕 NEW
- **Location:** `data/structures/pdb/6ZSL.pdb`
- **Composition:** NSP13 homodimer (Chains A,B, ~575-585 res each)
- **Citation:** [To be determined from RCSB]
- **Status:** 🆕 VALIDATED - Needs notebook creation
- **Target Interface:** NSP13-NSP13 (homodimer interface)
- **Note:** Homo-oligomeric interface (different from heterodimers)
- **Priority:** MEDIUM
- **Next Step:** Create `08_6ZSL_interface_discovery.ipynb`

---

### 5. 6XEZ - RTC-Helicase Complex ⚠️ COMPLEX
- **Location:** `data/structures/pdb/6XEZ.pdb`
- **Composition:** NSP12 (Chain A, 926 res) + NSP7 (Chain C, 73 res) + NSP8 (Chains B,D, 186 res) + NSP13 (Chains E,F, 596 res) + RNA (Chains P,T)
- **Citation:** Chen et al. (2020) Cell
- **Title:** Helicase-Polymerase Coupling in SARS-CoV-2 RTC
- **Release:** 2020-07-29
- **Status:** ⚠️ COMPLEX - Multiple interfaces available
- **Possible Targets:**
  - NSP12-NSP7 (compare with 7DFG)
  - NSP12-NSP8 (compare with 7DFG)
  - NSP13-NSP8 (helicase-RTC coupling) ← UNIQUE
  - NSP13-NSP12 (helicase-polymerase) ← UNIQUE
- **Priority:** LOW (use for comparison or advanced analysis)
- **Decision:** Analyze later after completing simpler structures

---

## Analysis Pipeline Status

| Structure | Validated | Notebook | Analysis | fpocket | Results |
|-----------|-----------|----------|----------|---------|---------|
| 6W4H | ✅ | ✅ | ✅ | ⏳ | ✅ |
| 7DFG | ✅ | ✅ | ⏳ | ⏳ | ⏳ |
| 7N0C | ✅ | ⏳ | ⏳ | ⏳ | ⏳ |
| 6ZSL | ✅ | ⏳ | ⏳ | ⏳ | ⏳ |
| 6XEZ | ✅ | ⏳ | ⏳ | ⏳ | ⏳ |

**Progress:** 1/5 complete (20%)

---

## Week 3-4 Timeline

### This Week (Week 3 - Feb 5-11)
- [x] Validate all structures
- [x] Replace incorrect structures (7EDI, 6W9C)
- [x] Add 7N0C and 6ZSL
- [ ] Run 7DFG analysis
- [ ] Create 7N0C notebook
- [ ] Create 6ZSL notebook

### Next Week (Week 4 - Feb 12-18)
- [ ] Complete all interface analyses
- [ ] Run fpocket on all structures
- [ ] Comparative analysis (NSP10 in 6W4H vs 7N0C)
- [ ] Rank targets by druggability
- [ ] Begin virtual screening prep

---

## Removed Structures (Learning Record)

### ❌ 7EDI - SPIKE PROTEIN (Not NSP14)
- **Was Expected:** NSP10-NSP14 ExoN
- **Actually Was:** Spike glycoprotein (B.1.1.7 variant)
- **Citation:** Yang et al. (2021) Nat. Struct. Mol. Biol.
- **Lesson:** Always validate before downloading

### ❌ 6W9C - PAPAIN-LIKE PROTEASE (Not NSP13)
- **Was Expected:** NSP13 Helicase
- **Actually Was:** NSP3 PLpro (papain-like protease)
- **Citation:** Osipiuk et al. (2020)
- **Lesson:** PDB IDs can be misleading

---

## Key Comparisons Enabled

1. **NSP10 Interfaces:** 6W4H (NSP10-NSP16) vs 7N0C (NSP10-NSP14)
   - Same protein (NSP10), different partners
   - Can identify conserved NSP10 binding modes

2. **RdRp Complexes:** 7DFG vs 6XEZ
   - Both have NSP12-NSP7-NSP8
   - 6XEZ additionally has NSP13 helicase

3. **Homodimer vs Heterodimer:**
   - 6ZSL: NSP13-NSP13 (homodimer)
   - All others: Heterodimeric interfaces

---

## Summary Statistics

- **Total Structures:** 5
- **Validated:** 5/5 (100%)
- **Analyzed:** 1/5 (20%)
- **Protein-Protein Interfaces:** 6-8 (depending on 6XEZ usage)
- **Target Proteins:** NSP10, NSP12, NSP13, NSP14, NSP16, NSP7, NSP8
- **Total Chains:** 22 chains across 5 structures

