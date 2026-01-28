# Paper Analysis: AI-guided pipeline for protein–protein interaction drug discovery identifies a SARS-CoV-2 inhibitor

**Date Read:** January 28, 2025  
**Reviewer:** Olivier Nsekuye  
**Paper Number:** 01

---

## Citation

**Authors:** Philipp Trepte, Christopher Secker, Julien Olivet, Jeremy Blavier, et al.  
**Title:** AI-guided pipeline for protein–protein interaction drug discovery identifies a SARS-CoV-2 inhibitor  
**Journal:** Molecular Systems Biology  
**Year:** 2024  
**Volume:** 20, Issue 4, Pages 428-457  
**DOI:** 10.1038/s44320-024-00019-8  
**PMID:** Not yet assigned  

---

## Summary (3-5 sentences)

This paper describes a comprehensive AI-guided pipeline combining experimental binary PPI mapping (LuTHy assay), machine learning-based interaction scoring (maSVM algorithm), AlphaFold-Multimer structure prediction, and ultra-large virtual screening (VirtualFlow) to discover PPI inhibitors. The authors mapped 29 high-confidence SARS-CoV-2 protein-protein interactions and targeted the NSP10-NSP16 methyltransferase complex interface. Through virtual screening of ~350 million compounds from the Enamine REAL library, they identified compound 459 that binds NSP10 (Kd ~13 µM), disrupts the NSP10-NSP16 interaction (IC50 9.2 µM), inhibits methyltransferase activity (~50%), and reduces SARS-CoV-2 replication (IC50 39.5 µM) without cytotoxicity. The study validates that PPI interfaces can be successfully targeted with small molecules and provides a complete blueprint for structure-based PPI drug discovery.

---

## Relevance to Project

**Relevance Score:** 5/5 ⭐⭐⭐⭐⭐

**Why relevant:**
- **DIRECTLY targets NSP10-NSP16** - one of my primary target complexes (PDB: 6W4H)
- **Validates the druggability** of the NSP10-NSP16 interface
- **Identifies critical hot spot residues:** Lys93 (NSP10) and Asp106 (NSP16)
- **Provides successful inhibitor** (compound 459) as benchmark for my work
- **Uses identical computational approach:** AlphaFold-Multimer + Virtual Screening
- **Demonstrates complete workflow** from structure prediction to experimental validation
- **Proves feasibility** of targeting coronavirus RTC PPI interfaces with small molecules

**How I will use this:**
1. **Target the same interface** (NSP10 Lys93 region) in my virtual screening
2. **Benchmark my predictions** against their validated hot spots
3. **Design grid box** based on their docking box coordinates
4. **Compare my AlphaFold-Multimer predictions** to their structures
5. **Use their validation strategy** (MTase assay → binding → cellular → antiviral)
6. **Cite as proof-of-concept** that RTC PPI interfaces are druggable
7. **Reference their workflow** in my Methods section

---

## Key Findings

### Main Results

1. **Developed maSVM algorithm** for scoring quantitative PPI data with improved sensitivity/specificity
2. **Mapped SARS-CoV-2 interactome:** Identified 29 high-confidence interactions (>95% probability) including known interactions (NSP8-NSP12, NSP10-NSP16, NSP10-NSP14) and 15 novel interactions
3. **AlphaFold-Multimer predictions:** Successfully predicted structures for 23 SARS-CoV-2 PPIs, validated 15 with >95% probability
4. **NSP10-NSP16 hot spots identified:**
   - **NSP10 Lys93** (lysine 93) - lowest ΔG, critical for binding
   - **NSP16 Asp106** (aspartate 106) - forms salt bridge with Lys93
   - Mutagenesis validation: K93E and D106K mutations strongly reduce binding
5. **Virtual screening results:**
   - Screened ~350 million compounds (Enamine REAL)
   - Targeted NSP10 interface at Lys93
   - Identified 15 synthesizable compounds
   - Top hit: **Compound 459**
6. **Compound 459 characterization:**
   - Binds NSP10: Kd = 12.97 µM (MST assay)
   - Disrupts NSP10-NSP16: IC50 = 9.2 µM (LuTHy-BRET)
   - Inhibits MTase activity: ~50% at 100 µM
   - Inhibits SARS-CoV-2 replication: IC50 = 39.5 µM
   - No cytotoxicity observed
   - **Additive effects** with AZ1 (USP25 inhibitor)

### Structural Information

**PDB IDs mentioned:**
- **6W4H** - SARS-CoV-2 NSP10-NSP16 complex (MAIN structure used)
- **6WVN** - Alternative SARS-CoV-2 NSP10-NSP16 structure
- **7DFG, 6XEZ** - NSP12-NSP7-NSP8 complex
- **7EDI, 7DIY** - NSP10-NSP14 complex
- **6YYT, 7EIZ** - NSP8-NSP12 complex
- **7PKU** - NSP3-N complex
- **6VYO** - Nucleocapsid N protein
- **6XDC** - ORF3a structure
- **6VYB** - Spike protein

**Resolution:** Not explicitly stated for 6W4H in this paper (refer to original structure paper)

**Chains/proteins:**
- NSP10: Amino acids 23-145 (construct used for protein production)
- NSP16: Full-length co-expressed with NSP10
- Complex: NSP10-NSP16 heterodimer

**Ligands co-crystallized:**
- SAM (S-adenosylmethionine) - cofactor for MTase activity
- Cap-0 RNA substrate used in MTase assays

### Interface Analysis

**Interface residues identified:**

**NSP10 side:**
- **Lys93 (K93)** - CRITICAL HOT SPOT
  - Lowest solvation-free energy (ΔG)
  - Forms contacts with NSP16 Asp106
  - Targeted in virtual screening
- Flexible residues for docking: Val, Met, Phe, Ser, Cys, Cys, Arg, His, Tyr, Lys, Lys, His (12 residues total)

**NSP16 side:**
- **Asp106 (D106)** - CRITICAL HOT SPOT
  - Lowest ΔG on NSP16 side
  - Forms salt bridge with NSP10 Lys93
  - Mutation D106K abolishes interaction

**Buried surface area:** Measured by PDBePISA (specific values in supplementary data)

**Key interactions:**
- **Salt bridge:** NSP10 Lys93 ↔ NSP16 Asp106 (PRIMARY interaction)
- Multiple supporting hydrophobic contacts
- Hydrogen bonding network

**Hot Spots:**

**Critical residues for binding (validated by mutagenesis):**
1. **NSP10 Lys93** → K93E mutation = strong reduction in binding
2. **NSP16 Asp106** → D106K mutation = strong reduction in binding

**Mutational data:**
- **NSP10 K93E:** BRET50 dramatically reduced (Fig. 5E)
- **NSP16 D106K:** BRET50 = 0.02 (essentially no binding)
- Both mutations maintain protein stability (no expression/folding issues)

**Conservation:**
- Interface region **highly conserved** among coronaviruses (cited from Lugari et al. 2010)
- Lys93 position conserved across coronavirus groups
- Makes interface attractive for **pan-coronavirus targeting**

### Druggability Assessment

**Pockets mentioned:**
- **NSP10 interface groove** around Lys93
- Characterized as "groove-shaped" rather than deep pocket
- Challenging but druggable (proven by compound 459)

**Druggability score (if given):**
- Not explicitly provided as numerical score
- **Implicit validation:** Successful identification of micromolar inhibitor

**Known inhibitors:**
1. **Compound 459 (this study):**
   - Kd = 12.97 µM (binding to NSP10)
   - IC50 = 9.2 µM (PPI disruption)
   - IC50 = 39.5 µM (antiviral)
   - Mechanism: Binds NSP10 near Lys93, blocks NSP16 interaction
   - Chemical structure: Heterocyclic scaffold (see Fig. 6D)
   - Synthesis: Enamine Ltd.

2. **Previous peptide inhibitors** (cited):
   - Wang et al. (2015) - NSP10-derived peptides
   - Ke et al. (2012) - Short peptides targeting interface
   - Generally lower potency and poor drug-like properties

3. **Other small molecules** (from screening):
   - 14 additional compounds tested (lower activity than 459)
   - 3 compounds showed significant MTase inhibition
   - Compound 459 most potent

4. **Combination therapy:**
   - **AZ1** (USP25 inhibitor) + Compound 459
   - Additive effect on SARS-CoV-2 replication
   - Suggests potential for multi-target approaches

---

## Figures/Tables of Interest

### Figure 5: Predicting NSP10-NSP16 complex and targeting interface by virtual screening

**5A - PAE heatmap:**
- Very low predicted alignment error
- High confidence in interface prediction
- **Key insight:** AlphaFold-Multimer reliably predicts this complex

**5B - Structural overlay:**
- 5 AlphaFold models vs experimental structure (6W4H)
- Excellent overlap
- **Key insight:** AFM predictions match experimental structure

**5C - ΔG per residue:**
- **CRITICAL GRAPH**
- Shows solvation-free energy for each amino acid
- **Lys93 and Asp106 have lowest ΔG**
- **Action:** Extract these exact residues in my structures

**5D - Close-up of Lys93-Asp106 interaction:**
- Shows hydrogen bonds and salt bridge
- **Action:** Visualize this in PyMOL for 6W4H

**5E - Mutagenesis validation:**
- BRET50 curves for WT vs mutants
- K93E: BRET50,max = 0.06 (vs 0.16 for WT)
- D106K: BRET50,max = 0.02 (essentially no binding)
- **Key insight:** Validates computational hot spot predictions

**5F - Docking box visualization:**
- **CRITICAL IMAGE**
- Shows exact docking box placement on NSP10
- Centered on Lys93
- **Action:** Use similar positioning for my grid box

**5G - Virtual screening workflow:**
- 350M ligands → 10M re-docking → clustering → 20 compounds → 15 synthesized
- **Action:** Follow similar workflow for my screening

**5H - Docking scores:**
- Top 100 hits: mean ~-8.5 kcal/mol
- **Action:** Compare my docking scores to these values

### Figure 6: Compound 459 validation

**6A - MTase assay schematic:**
- Cap-0 to Cap-1 conversion
- SAM to SAH
- Luminescence readout
- **Action:** Understand this for future validation planning

**6B - MTase inhibition:**
- Heatmap of 15 compounds tested
- 3 hits: 459, 468, 469
- Compound 459: ~50% inhibition at 100 µM
- **Benchmark:** This is level of inhibition to aim for

**6C - Compound 459 docked:**
- Shows binding pose on NSP10
- Near Lys93 as designed
- **Action:** Compare my top hits to this binding mode

**6D - Chemical structure:**
- Heterocyclic scaffold with phenyl groups
- **Action:** Note structural features for future hit optimization

**6E - MST assay principle:**
- Microscale thermophoresis schematic
- **Action:** Consider this assay for future validation

**6F - MST traces:**
- Fluorescence curves at different compound concentrations
- **Action:** Understand readout for binding assays

**6G - Binding curve:**
- **Kd = 12.97 µM**
- **Benchmark:** My hits should target similar or better affinity

**6H - Cellular PPI disruption:**
- **IC50 = 9.2 µM** (LuTHy-BRET)
- Concentration-dependent reduction
- **Benchmark:** Cellular activity in low micromolar range

**6I - Antiviral assay:**
- **IC50 = 39.5 µM** (SARS-CoV-2 replication)
- Nanoluciferase reporter
- **Benchmark:** Proof that PPI disruption → antiviral effect

**6J - Combination with AZ1:**
- Additive effects observed
- Statistical significance shown
- **Key insight:** Multi-target approaches may be beneficial

### Supplementary Figures Relevant:

**Appendix Fig S7:**
- BRET50 values correlate with ΔG
- Validates computational predictions
- **Action:** Check correlation in my analysis

**Fig EV3:**
- Additional SARS-CoV-2 interactions validated
- Shows broader applicability
- **Context:** NSP10-NSP16 is one of many RTC PPIs

---

## Methods Relevant to My Project

### Techniques used:

1. **AlphaFold-Multimer prediction:**
   - ColabFold interface (version 1.2.0 or 1.3.0)
   - Parameters:
use_amber: 'no'
 template_mode: 'none'
 msa_mode: 'MMseq2 (UniRef+Environmental)'
 pair_mode: 'unpaired+paired'
 model_type: 'auto' (AlphaFold2-multimer-v2)
 num_models: 5
 num_recycles: 3
 rank_by: 'auto'
 stop_at_score: 100

- **Action:** Use these exact parameters for my predictions

2. **PDBePISA analysis:**
   - Extracts interface area (iA) and ΔG values
   - Python wrapper: PisaPy available
   - **Action:** Install and use PDBePISA for my structures

3. **Virtual screening (VirtualFlow):**
   - **Primary docking:** Quick Vina 2 (exhaustiveness = 1)
   - **Re-docking:** AutoDock Vina + Smina Vinardo (exhaustiveness = 1)
   - **Docking box:** 75.647 × 16.822 × 17.631 Å
   - **Flexible residues:** 12 amino acids at interface
   - Library: Enamine REAL (~350M compounds)
   - **Action:** Consider VirtualFlow if accessible, otherwise AutoDock Vina locally

4. **maSVM machine learning algorithm:**
   - Multi-adaptive support vector machine
   - Trained on reference sets (positive and random pairs)
   - Predicts interaction probability
   - **Note:** More relevant for experimental PPI mapping, less for my project

5. **Mutagenesis validation:**
   - Site-directed mutagenesis (K93E, D106K)
   - LuTHy-BRET donor saturation assays
   - **Action:** Understand method for interpreting hot spot validation

6. **Microscale thermophoresis (MST):**
   - Fluorescently labeled NSP10
   - Measures binding affinity (Kd)
   - **Action:** Potential future validation method

7. **Methyltransferase assay (MTase-Glo):**
   - Detects SAM → SAH conversion
   - Luminescence-based
   - Primary screen for inhibitor activity
   - **Action:** Critical assay to plan for experimental validation

8. **Cellular antiviral assay:**
   - icSARS-CoV-2-nanoluciferase reporter
   - HEK293-ACE2 cells
   - Measures viral replication inhibition
   - **Action:** Ultimate validation step

### Could I use these methods?

**YES - Immediately applicable:**
- ✅ AlphaFold-Multimer prediction (ColabFold)
- ✅ PDBePISA interface analysis
- ✅ AutoDock Vina docking (similar to VirtualFlow)
- ✅ Flexible docking approach
- ✅ Chemical clustering for hit selection

**MAYBE - Depends on resources:**
- ⚠️ VirtualFlow platform (need HPC access, already have JURECA)
- ⚠️ Enamine REAL library (~350M compounds, may be too large initially)

**NO - Requires experimental facilities:**
- ❌ LuTHy assay (Prof. Twizere's lab could do this)
- ❌ MST assay (need equipment)
- ❌ MTase assay (need reagents)
- ❌ Antiviral assay (need BSL-3 facility)

**Recommendation:**
- Use computational methods NOW (AlphaFold, PDBePISA, Docking)
- Plan experimental validation with Prof. Twizere's lab for LATER
- Focus on identifying top computational hits first

---

## Questions Raised

1. **How to define optimal docking box size for groove-shaped interfaces?**
   - Their box: 75.647 × 16.822 × 17.631 Å (large and asymmetric)
   - Is 25 Å cubic box sufficient or should I use larger/asymmetric?

2. **What is the threshold for "good" docking scores for PPI interfaces?**
   - Their top 100: ~-8.5 kcal/mol
   - How does this compare to typical enzyme active site docking?

3. **Should I also target NSP16 side (Asp106) in parallel screening?**
   - They focused on NSP10 Lys93 side
   - Might NSP16 side offer alternative binding pockets?

4. **How to handle flexible docking computationally?**
   - They used 12 flexible residues
   - Significant computational cost increase
   - Worth it for initial screening or only for re-docking top hits?

5. **What chemical properties should I filter for PPI inhibitors?**
   - Compound 459 properties not fully detailed
   - Need PPI-specific filters (likely larger MW, more hydrophobic than standard drugs)

6. **How was the Enamine REAL library accessed?**
   - Academic license? (they mention free academic access)
   - Can I get access for my project?

7. **What is the relationship between in vitro Kd and cellular IC50?**
   - Kd = 13 µM vs IC50 = 9.2 µM (close but IC50 lower)
   - Cell permeability factors?

8. **Why does antiviral IC50 (39.5 µM) differ from PPI IC50 (9.2 µM)?**
   - Cellular context differences?
   - Additional cellular factors affecting viral replication?

---

## Action Items

### Immediate (Week 2):
- [x] Read and annotate this paper completely
- [ ] **Download 6W4H structure** from PDB
- [ ] **Visualize Lys93-Asp106 interaction** in PyMOL
- [ ] **Extract exact coordinates** of Lys93 and surrounding residues
- [ ] **Identify all residues within 10 Å** of Lys93 for grid box definition
- [ ] **Compare 6W4H to my AlphaFold-Multimer predictions** (when generated)
- [ ] **Document hot spot residues** in summary table

### Week 3-4:
- [ ] **Install PDBePISA** or use web interface
- [ ] **Run PDBePISA analysis** on 6W4H to extract iA and ΔG values
- [ ] **Validate** that Lys93 and Asp106 have lowest ΔG in my analysis
- [ ] **Define docking box** centered on Lys93 (similar to their approach)
- [ ] **Test docking box size** (25 Å vs larger asymmetric box)
- [ ] **Identify 12 flexible residues** at interface (if using flexible docking)

### Week 5-7:
- [ ] **Apply for Enamine REAL academic license**
- [ ] **Set up docking workflow** (AutoDock Vina or VirtualFlow)
- [ ] **Test docking protocol** with known inhibitors (compound 459 if structure available)
- [ ] **Benchmark my docking scores** against their -8.5 kcal/mol threshold

### Month 2+:
- [ ] **Run virtual screening** on selected library
- [ ] **Chemical clustering** of top hits
- [ ] **Select representative compounds** (15-30 like they did)
- [ ] **Coordinate with Prof. Twizere** for experimental validation
- [ ] **Plan MTase assay** as primary experimental screen

### For Thesis:
- [ ] **Cite this paper** in Introduction (proof PPI interfaces druggable)
- [ ] **Reference in Methods** (AlphaFold parameters, docking approach)
- [ ] **Compare results** (my hot spots vs Lys93/Asp106)
- [ ] **Benchmark hits** (my docking scores vs their -8.5 kcal/mol; my predicted affinities vs 13 µM)
- [ ] **Use validation strategy** as template for future work

---

## Related Papers to Read

### PRIORITY 1 - Structural papers for my targets:

1. **Rosas-Lemus et al. (2020)** - SARS-CoV-2 NSP10-NSP16 2'-O-methyltransferase structure
   - PDB: 6W4H, 6WVN
   - **Sci Signal 13:eabe1202**
   - **Action:** READ NEXT - This is the main structural paper for 6W4H

2. **Lin et al. (2021)** - SARS-CoV-2 NSP10-NSP14 exonuclease structure
   - PDB: 7EDI
   - **Nucleic Acids Res 49:gkab320**
   - **Action:** High priority for NSP10-NSP14 interface

3. **Gao et al. (2020)** - NSP12-NSP7-NSP8 RdRp complex structure
   - PDB: 7DFG (likely, need to verify)
   - **Action:** High priority for NSP12 interfaces

### PRIORITY 2 - Functional and mechanistic:

4. **Decroly et al. (2011)** - NSP10/NSP16 MTase functional analysis
   - **PLoS Pathog 7:e1002059**
   - **Action:** Understand enzyme mechanism for validation assays

5. **Chen et al. (2011)** - SARS-CoV-1 NSP10/NSP16 mechanism
   - **PLoS Pathog 7:e1002294**
   - **Action:** Compare SARS-CoV-1 vs CoV-2 interface conservation

6. **Daffis et al. (2010)** - Importance of 2'-O methylation for viral replication
   - **Nature 468:452-456**
   - **Action:** Understand why targeting NSP10-NSP16 blocks replication

### PRIORITY 3 - Previous inhibitor studies:

7. **Wang et al. (2015)** - NSP10-derived peptide inhibitors
   - **J Virol 89:8416-8427**
   - **Action:** Learn from previous inhibitor approaches

8. **Ke et al. (2012)** - Short peptides targeting NSP10/NSP16
   - **Virus Res 167:322-328**
   - **Action:** Compare to compound 459 approach

9. **Lugari et al. (2010)** - NSP10-NSP16 interface conservation
   - **J Biol Chem 285:33230-33241**
   - **Action:** Understand pan-coronavirus conservation

### PRIORITY 4 - Virtual screening methodology:

10. **Gorgulla et al. (2020)** - VirtualFlow platform
    - **Nature 580:663-668**
    - **Action:** Understand VirtualFlow if planning to use it

11. **Gorgulla et al. (2021)** - SARS-CoV-2 virtual screening (includes NSP10-NSP16)
    - **iScience 24:102021**
    - **Action:** Compare their NSP10-NSP16 hits to compound 459

### PRIORITY 5 - AlphaFold and structure prediction:

12. **Evans et al. (2022)** - AlphaFold-Multimer
    - **Preprint bioRxiv**
    - **Action:** Technical details on AFM algorithm

13. **Jumper et al. (2021)** - AlphaFold2
    - **Nature 596:583-589**
    - **Action:** Background on structure prediction

---

## Notes

### Strengths of this paper:
1. **Complete pipeline** - from PPI identification to validated inhibitor
2. **Rigorous validation** - multiple orthogonal assays (MST, BRET, MTase, antiviral)
3. **Mutagenesis confirmation** - proves hot spots are functionally important
4. **Large-scale screening** - 350M compounds, not just focused library
5. **Open methods** - computational code available on GitHub
6. **Clinically relevant** - shows antiviral activity, suggests therapeutic potential

### Limitations:
1. **Micromolar potency** - Compound 459 not optimized (Kd ~13 µM)
2. **Moderate inhibition** - 50% MTase inhibition, not complete block
3. **No crystal structure of compound 459 bound** - would validate binding mode
4. **Limited SAR** - only 15 compounds tested, no optimization shown
5. **In vitro/cellular gap** - antiviral IC50 (39.5 µM) higher than PPI IC50 (9.2 µM)
6. **No in vivo data** - mouse model would strengthen translational potential

### Why this matters for my work:
- **Proof of concept:** RTC PPI interfaces CAN be drugged
- **Validated target:** NSP10-NSP16 is druggable and functionally important
- **Benchmark:** Provides comparison point for my computational predictions
- **Roadmap:** Shows complete path from computation to validation
- **Hot spots:** Gives me specific residues to focus on (Lys93, Asp106)
- **Credibility:** Published in high-quality journal (MSB), validates approach

### Key takeaway:
**"If they can find a micromolar inhibitor for this challenging PPI interface using virtual screening, I can too - and maybe improve upon it by targeting multiple interfaces or identifying better scaffolds."**

---

**Last Updated:** January 28, 2025
