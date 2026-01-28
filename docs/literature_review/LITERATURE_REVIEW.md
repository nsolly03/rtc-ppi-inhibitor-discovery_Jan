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
### NSP9 (RTC RNA-binding component)

#### Paper 2: Nanobodies against SARS-CoV-2 non-structural protein Nsp9 inhibit viral replication by targeting innate immunity

- **Authors:** Venit, Blavier, Maseko, Shu, Esposito, **Twizere**, Percipalle, et al.
- **Journal/Year:** bioRxiv preprint, 2023 (NOT PEER-REVIEWED)
- **PDB Structures:** None new (refers to previous NSP9 structures)
- **Key Findings:**
  - Developed anti-NSP9 nanobody (2NSP23) delivered as LNP-encapsulated mRNA
  - **EC50 = 1.8 µM** (comparable to Remdesivir), CC50 >30 µM
  - >90% inhibition of Wuhan, Delta, Mu, Omicron; ~60% inhibition of Beta variant
  - Mechanism: Stabilizes non-functional tetrameric NSP9 (functional form is monomer/dimer)
  - Prevents RTC assembly → blocks viral RNA replication
  - RNA-seq: Rescues host cell transcriptome disrupted by infection
  - **Activates innate immune response genes BEFORE infection** (prophylactic potential)
  
- **NSP9 Function:**
  - RNA-binding protein essential for RTC assembly
  - Functions as monomer/dimer in replication complex
  - Dimerization interface overlaps with NSP12 binding site
  - Tetramerization = non-functional state
  
- **Interface Residues:**
  - **NOT detailed in this paper** - refers to Esposito et al. 2021 for epitope mapping
  - Nanobody binds at/near dimerization interface
  - Binding prevents functional oligomeric state
  
- **Hot Spots:**
  - **NOT explicitly defined** - need to read Esposito et al. 2021
  - Dimerization interface residues critical for function
  
- **Notes:**
  - **PREPRINT** - not yet peer-reviewed
  - **Prof. Twizere is co-author** - direct collaboration opportunity!
  - Validates RTC targeting approach
  - NSP9 highly conserved across coronaviruses (low mutation rate)
  - Different approach (nanobody) but proves concept
  - Dual mechanism: RTC disruption + immune activation
  - Could be alternative/complementary target to NSP10-NSP16
  - Immune genes activated: RSAD2, OAS1/2/3, IFIT1/2/3, ISG15, MX1
  - Mitochondrial genes rescued: COX5B, NDUFS8, MRPL9, MRPL24
  - See detailed notes: `notes/Paper_02_Venit_2023_NSP9_Nanobody_Inhibitor.md`
  - **MUST READ:** Esposito et al. 2021 for structural details

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


| Target | Compound | IC50 / Kd | Binding Mode | Reference |
|--------|----------|-----------|--------------|-----------|
| NSP10-NSP16 | Compound 459 | Kd 12.97 µM (NSP10), IC50 9.2 µM (PPI), IC50 39.5 µM (antiviral) | Binds NSP10 near Lys93, blocks NSP16 interaction | Trepte et al. 2024 |
| NSP10-NSP16 | NSP10-derived peptides | Variable (µM-mM) | Interface mimicry | Wang et al. 2015 |
| NSP10-NSP16 | Short peptides | Variable (µM-mM) | Interface disruption | Ke et al. 2012 |
| **NSP9** | **2NSP23 Nanobody** | **EC50 1.8 µM (antiviral)** | **Stabilizes non-functional NSP9 tetramer, prevents RTC assembly** | **Venit et al. 2023 (preprint)** |
| **NSP9** | **2NSP90 Nanobody** | **Not detailed** | **Similar to 2NSP23** | **Venit et al. 2023 (preprint)** |



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

**NSP9 (RTC RNA-binding protein):**
- **VALIDATED AS DRUGGABLE** (Venit et al. 2023 preprint)
- 2NSP23 nanobody demonstrates NSP9 can be targeted effectively
- Mechanism: Stabilization of non-functional tetrameric state
- Oligomerization interface is druggable
- EC50 1.8 µM comparable to approved antivirals (Remdesivir)
- Target site: **Dimerization interface** (overlaps with NSP12 binding)
- Success with nanobodies suggests small molecules could work too
- Dual mechanism possible: Direct RTC disruption + immune activation

### Conservation Analysis

**NSP10-NSP16 Interface:**
- **Highly conserved across coronaviruses** (Trepte et al. 2024, Lugari et al. 2010)
- Lys93 position conserved across coronavirus groups
- Pan-coronavirus targeting potential
- Low mutation frequency expected (interface under structural constraint)

**NSP9:**
- **Highly conserved across coronaviruses** (Venit et al. 2023, Abbasian et al. 2023)
- **"Extremely low degree of mutagenicity"**
- 2NSP23 works on multiple variants (>90% inhibition)
- Beta variant slightly less susceptible (~60%)
- Pan-coronavirus targeting potential
- Less prone to resistance than Spike protein


### Structural Features

**NSP10-NSP16 Complex (PDB: 6W4H):**
- AlphaFold-Multimer predictions **match experimental structure** with high accuracy
- Low PAE (Predicted Alignment Error) at interface
- Interface dominated by Lys93-Asp106 salt bridge
- Supporting hydrophobic interactions
- Flexible interface residues (12 residues) used in virtual screening
- MTase activity dependent on NSP10-NSP16 interaction


**NSP9:**
- Small RNA-binding protein (~12 kDa)
- Oligomerization-dependent function:
  - **Monomer/dimer = FUNCTIONAL** for RTC assembly
  - **Tetramer = NON-FUNCTIONAL** (stabilized by 2NSP23)
- Dimerization interface = NSP12 binding site
- Located in RTC catalytic center
- Specific epitope mapping by NMR (Esposito et al. 2021)
- Strong tendency to oligomerize in vitro

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



**For NSP9 Targeting:**
- Nanobody approach different from small molecules
- Need structural details of epitope (read Esposito et al. 2021)
- Oligomerization interface may be flat/large
- Beta variant shows partial resistance
- Clinical translation of mRNA-LNP delivery uncertain
- Preprint status - not yet peer-reviewed

**Opportunities:**
- Validated druggable target with experimental proof
- Dual mechanism: RTC + immune activation
- Pan-coronavirus potential (conservation)
- Could combine with NSP10-NSP16 targeting
- Prof. Twizere collaboration possible
- Prophylactic potential (immune priming)


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
