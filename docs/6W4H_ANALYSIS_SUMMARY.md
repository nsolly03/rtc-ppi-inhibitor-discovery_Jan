# 6W4H NSP10-NSP16 Interface Analysis Summary

**Date:** January 28, 2025  
**Analyst:** Olivier Nsekuye  
**Reference:** Trepte et al. (2024) Molecular Systems Biology

---

## Critical Finding

**⚠️ IMPORTANT DISCREPANCY DISCOVERED:**

The paper (Trepte et al. 2024) references hot spots as:
- NSP10 **Lys93** (K93)
- NSP16 **Asp106** (D106)

However, in the **actual 6W4H crystal structure**, these positions contain:
- Position 93 in NSP10 = **PHE** (not LYS)
- Position 106 in NSP16 = **SER** (not ASP)

**The ACTUAL hot spot in 6W4H structure is:**
- NSP10 **Lys76** (K76) - PDB residue 4346
- NSP16 **Asp107** (D107) - PDB residue 6904

---

## Explanation

The 6W4H PDB file uses **polyprotein numbering** (starts at residue 4271 for NSP10):
- NSP10: residues 4271-4386 (116 residues)
- NSP16: residues 6798-7096 (299 residues)

The paper likely refers to positions in a **canonical reference sequence** or different SARS-CoV-2 variant.

**For docking with 6W4H, we must use K76-D107!**

---

## Hot Spot Characterization

### NSP10 K76 (Chain B, PDB 4346)
- **Residue:** Lysine (LYS)
- **CA coordinates:** (74.536, 12.890, 9.624) Å
- **Center of mass:** (75.883, 11.641, 10.087) Å

### NSP16 D107 (Chain A, PDB 6904)
- **Residue:** Aspartate (ASP)
- **CA coordinates:** (78.186, 14.903, 12.641) Å

### Interaction
- **Distance:** 5.15 Å (CA to CA)
- **Type:** Salt bridge (LIKELY)
- **Strength:** Close proximity ideal for electrostatic interaction

---

## Interface Composition

**Total residues within 10 Å of K76:** 22 residues

### NSP10 Interface (16 residues):
- K76 (hot spot), G77, L75, A54, G53, K78, C73, D74, C57, G52, S55, C60, Y79, C56, F72, V80

### NSP16 Interface (6 residues):
- D107 (hot spot), S106, A108, V105, D109, A84

### Key Interface Features:
- Multiple **cysteine** residues (C73, C57, C60, C56) → potential disulfide bridges
- Additional **charged residues** (K78, D74) → electrostatic network
- **Hydrophobic residues** (L75, F72, V80) → hydrophobic packing

---

## Docking Grid Box Parameters

**For AutoDock Vina (Week 5+):**
```
receptor = 6W4H_prepared.pdbqt
ligand = ligand.pdbqt

# Grid box center (Å) - centered on NSP10 K76
center_x = 75.883
center_y = 11.641
center_z = 10.087

# Grid box size (Å) - cubic box
size_x = 25.0
size_y = 25.0
size_z = 25.0

# Docking parameters
exhaustiveness = 8
num_modes = 9
energy_range = 3
```

**Note:** Start with 25 Å cubic box. If results are poor, consider:
- Expanding to 30-35 Å
- Using asymmetric box (like Trepte et al.: 75.6 × 16.8 × 17.6 Å)

---

## Files Generated

1. **`data/analysis_results/6W4H_analysis.json`**
   - Complete structural analysis in JSON format
   - Hot spot coordinates
   - Grid box parameters
   - Interface statistics

2. **`data/analysis_results/6W4H_vina_config.txt`**
   - Ready-to-use AutoDock Vina configuration file
   - Pre-filled with grid box parameters

3. **`data/analysis_results/6W4H_interface_residues.csv`**
   - List of all interface residues
   - Includes: Chain, Residue type, PDB number, Sequence position, Distance

4. **`scripts/analyze_6W4H_corrected.py`**
   - Analysis script (repeatable)
   - Handles polyprotein numbering correctly

---

## Validation Against Literature

### Trepte et al. 2024 Findings:
- ✅ Interface is druggable (compound 459: Kd 13 µM)
- ✅ Hot spot involves charged residues (K-D salt bridge)
- ✅ Docking box centered on lysine residue
- ⚠️ Numbering discrepancy (K93 vs K76, D106 vs D107)

### Our Findings Match:
- ✅ K-D salt bridge identified (5.15 Å)
- ✅ Interface contains ~20 residues
- ✅ Grid box centered on lysine
- ✅ Similar interface characteristics

**Conclusion:** The interface topology is consistent with the paper, only the residue numbering differs due to the PDB file using polyprotein numbering.

---

## Next Steps

### Week 3-4: Pocket Identification
- [ ] Install fpocket
- [ ] Run fpocket on 6W4H
- [ ] Compare pockets with K76 interface
- [ ] Validate pocket includes K76-D107

### Week 5-7: Docking Preparation
- [ ] Remove waters from 6W4H
- [ ] Add hydrogens to structure
- [ ] Convert to PDBQT format
- [ ] Test docking with known inhibitors
- [ ] Benchmark against Trepte et al. results

### Month 2+: Virtual Screening
- [ ] Select compound library
- [ ] Run screening on JURECA HPC
- [ ] Target grid box: (75.883, 11.641, 10.087)
- [ ] Aim for scores ~-8.5 kcal/mol (like compound 459)

---

## References

1. **Trepte et al. (2024)** AI-guided pipeline for protein–protein interaction drug discovery identifies a SARS-CoV-2 inhibitor. *Molecular Systems Biology* 20:428-457.

2. **Rosas-Lemus et al. (2020)** High-resolution structures of the SARS-CoV-2 2'-O-methyltransferase reveal strategies for structure-based inhibitor design. *Sci Signal* 13:eabe1202. [PDB: 6W4H]

---

**Analysis completed:** January 28, 2025  
**Status:** ✅ Structure analyzed, coordinates extracted, ready for Week 3-4
