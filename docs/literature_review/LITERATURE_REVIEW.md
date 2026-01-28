# Literature Review: Coronavirus RTC Protein-Protein Interactions

**Author:** Olivier Nsekuye  
**Project:** RTC-PPI Inhibitor Discovery  
**Started:** January 27, 2025

---

## Review Objectives

1. Understand RTC protein structures and functions
2. Identify protein-protein interface residues
3. Find known binding sites and inhibitors (if any)
4. Determine conservation across coronaviruses
5. Identify hot spots and druggable pockets

---

## Papers by Target

### NSP12-NSP7-NSP8 (RdRp Complex)

#### Paper 1: [Title]
- **Authors:**
- **Journal/Year:**
- **PDB Structures:**
- **Key Findings:**
- **Interface Residues:**
- **Hot Spots:**
- **Notes:**

---

### NSP10-NSP14 (ExoN Complex)

#### Paper 1: [Title]
- **Authors:**
- **Journal/Year:**
- **PDB Structures:**
- **Key Findings:**
- **Interface Residues:**
- **Hot Spots:**
- **Notes:**

---

### NSP10-NSP16 (MTase Complex)

#### Paper 1: AI-guided pipeline for protein–protein interaction drug discovery identifies a SARS-CoV-2 inhibitor

- **Authors:** Trepte, Secker, Olivet, Blavier, et al.
- **Journal/Year:** Molecular Systems Biology, 2024
- **PDB Structures:** 6W4H (main), 6WVN (alternative)
- **Key Findings:**
  - Developed complete AI-guided pipeline: binary PPI mapping + machine learning + AlphaFold-Multimer + virtual screening
  - Mapped 29 high-confidence SARS-CoV-2 PPIs
  - Targeted NSP10-NSP16 methyltransferase complex interface
  - Virtual screening of ~350M compounds (Enamine REAL)
  - Identified compound 459: Kd 12.97 µM (NSP10 binding), IC50 9.2 µM (PPI disruption), IC50 39.5 µM (antiviral)
  - ~50% inhibition of methyltransferase activity at 100 µM
  - No cytotoxicity observed
  
- **Interface Residues:**
  - **NSP10:** **Lys93 (critical hot spot)**, Val, Met, Phe, Ser, Cys, Arg, His, Tyr (flexible interface residues)
  - **NSP16:** **Asp106 (critical hot spot)** - forms salt bridge with NSP10 Lys93
  
- **Hot Spots:**
  - **NSP10 Lys93** - lowest ΔG, critical for binding (K93E mutation abolishes interaction)
  - **NSP16 Asp106** - lowest ΔG on NSP16 side (D106K mutation abolishes interaction)
  - Salt bridge Lys93-Asp106 is PRIMARY interaction
  
- **Notes:**
  - **HIGHLY RELEVANT** - directly targets one of my primary complexes
  - Provides validated hot spots to focus my docking (Lys93 region)
  - Compound 459 serves as benchmark (micromolar affinity)
  - Proves NSP10-NSP16 interface is druggable with small molecules
  - Complete workflow from structure prediction to experimental validation
  - Interface highly conserved across coronaviruses → pan-coronavirus potential
  - Virtual screening used VirtualFlow platform
  - Docking box centered on NSP10 Lys93: 75.647 × 16.822 × 17.631 Å
  - 12 flexible residues used in re-docking
  - Top hits: ~-8.5 kcal/mol docking scores
  - See detailed notes: `notes/Paper_01_Trepte_2024_NSP10-NSP16_Inhibitor.md`
---

## Summary Tables

### Interface Residues by Target


| Target | Protein 1 Residues | Protein 2 Residues | Interface Size | Conservation | Reference |
|--------|-------------------|-------------------|----------------|--------------|-----------|
| NSP10-NSP16 | **Lys93 (hot spot)**, Val, Met, Phe, Ser, Cys, Arg, His, Tyr | **Asp106 (hot spot)** | Measured by PDBePISA | High (pan-coronavirus) | Trepte et al. 2024 |



| Target | Protein 1 Residues | Protein 2 Residues | Interface Size | Conservation |
|--------|-------------------|-------------------|----------------|--------------|
| NSP12-NSP7 | | | | |
| NSP12-NSP8 | | | | |
| NSP10-NSP14 | | | | |
| NSP10-NSP16 | | | | |

### Hot Spots Identified

| Target | Residue | Type | Contribution | Mutational Data | Reference |
|--------|---------|------|--------------|-----------------|-----------|
| NSP10-NSP16 | NSP10 Lys93 | Salt bridge | Lowest ΔG (critical) | K93E → binding abolished | Trepte et al. 2024 |
| NSP10-NSP16 | NSP16 Asp106 | Salt bridge | Lowest ΔG (critical) | D106K → binding abolished | Trepte et al. 2024 |


| Target | Residue | Type | Contribution | Mutational Data |
|--------|---------|------|--------------|-----------------|
| | | | | |

### Known Inhibitors

| Target | Compound | IC50 / Kd | Binding Mode | Reference |
|--------|----------|-----------|--------------|-----------|
| NSP10-NSP16 | Compound 459 | Kd 12.97 µM (NSP10), IC50 9.2 µM (PPI), IC50 39.5 µM (antiviral) | Binds NSP10 near Lys93, blocks NSP16 interaction | Trepte et al. 2024 |
| NSP10-NSP16 | NSP10-derived peptides | Variable (µM-mM) | Interface mimicry | Wang et al. 2015 |
| NSP10-NSP16 | Short peptides | Variable (µM-mM) | Interface disruption | Ke et al. 2012 |

| Target | Compound | IC50 | Binding Mode | Reference |
|--------|----------|------|--------------|-----------|
| | | | | |

---

## Key Insights

### Druggability Assessment

**NSP10-NSP16 Complex:**
- **VALIDATED AS DRUGGABLE** (Trepte et al. 2024)
- Compound 459 demonstrates small molecules CAN bind PPI interface
- Interface characterized as "groove-shaped" - challenging but accessible
- Critical hot spot: NSP10 Lys93 ↔ NSP16 Asp106 salt bridge
- Micromolar inhibitor achievable through virtual screening
- Target binding site: NSP10 interface around Lys93
- Docking box: 75.647 × 16.822 × 17.631 Å (large, asymmetric)
- Successful hits: docking scores ~-8.5 kcal/mol

### Conservation Analysis

**NSP10-NSP16 Interface:**
- **Highly conserved across coronaviruses** (Trepte et al. 2024, Lugari et al. 2010)
- Lys93 position conserved across coronavirus groups
- Pan-coronavirus targeting potential
- Low mutation frequency expected (interface under structural constraint)

### Structural Features

**NSP10-NSP16 Complex (PDB: 6W4H):**
- AlphaFold-Multimer predictions **match experimental structure** with high accuracy
- Low PAE (Predicted Alignment Error) at interface
- Interface dominated by Lys93-Asp106 salt bridge
- Supporting hydrophobic interactions
- Flexible interface residues (12 residues) used in virtual screening
- MTase activity dependent on NSP10-NSP16 interaction

### Challenges Identified

**For NSP10-NSP16 Targeting:**
- Interface is groove-shaped (not deep pocket) - more challenging for small molecule binding
- Micromolar affinity achieved (compound 459: 13 µM Kd) - needs optimization to nanomolar
- Moderate enzymatic inhibition (50% at 100 µM) - not complete block
- Gap between cellular PPI IC50 (9.2 µM) and antiviral IC50 (39.5 µM)
- PPI interfaces generally harder to drug than enzyme active sites
- Large interface area requires larger, more complex molecules
- Charged interactions (Lys93-Asp106) may require polar scaffolds

**Opportunities:**
- Validated druggable target with experimental proof
- Multiple validation assays established (MST, BRET, MTase, antiviral)
- Ultra-large virtual screening feasible (~350M compounds)
- Combination therapy potential (compound 459 + AZ1 additive effects)
- Pan-coronavirus potential due to conservation
### Druggability Assessment
- 

### Conservation Analysis
- 

### Structural Features
- 

### Challenges Identified
- 

---

**Last Updated:** January 27, 2025
