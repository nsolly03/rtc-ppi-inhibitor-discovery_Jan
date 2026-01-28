# Paper Analysis: Nanobodies against SARS-CoV-2 non-structural protein Nsp9 inhibit viral replication by targeting innate immunity

**Date Read:** January 28, 2025  
**Reviewer:** Olivier Nsekuye  
**Paper Number:** 02

---

## Citation

**Authors:** Tomas Venit, Jeremy Blavier, Sibusiso B. Maseko, Sam Shu, Lilia Espada, Christopher Breunig, Hans-Peter Holthoff, Sabrina C. Desbordes, Martin Lohse, Gennaro Esposito, Jean-Claude Twizere, Piergiorgio Percipalle  
**Title:** Nanobodies against SARS-CoV-2 non-structural protein Nsp9 inhibit viral replication by targeting innate immunity  
**Journal:** bioRxiv preprint (NOT YET PEER-REVIEWED)  
**Year:** 2023  
**DOI:** 10.1101/2023.10.12.561992  
**Posted:** October 13, 2023  

---

## Summary (3-5 sentences)

This preprint describes the development and validation of anti-NSP9 nanobodies delivered as mRNA encapsulated in lipid nanoparticles (LNPs) to inhibit SARS-CoV-2 replication. The lead nanobody, 2NSP23, was previously characterized by NMR and shown to stabilize non-functional tetrameric assemblies of NSP9, thereby preventing RTC assembly. When delivered as LNP-mRNA-2NSP23 to infected cells, the nanobody achieves EC50 of 1.8 µM (comparable to Remdesivir) and inhibits >90% replication of multiple SARS-CoV-2 variants including Wuhan, Delta, Mu, and Omicron (Beta slightly less at ~60%). RNA-seq analysis reveals that the nanobody rescues host cell gene expression programs disrupted by infection (particularly mitochondrial function) and remarkably activates innate immune response genes even in uninfected cells, suggesting a dual mechanism of action involving both direct RTC disruption and immune system priming.

---

## Relevance to Project

**Relevance Score:** 4/5 ⭐⭐⭐⭐

**Why relevant:**
- **Validates NSP9 as druggable RTC component** - proves targeting RTC proteins blocks viral replication
- **Different RTC target than NSP10-NSP16** - NSP9 is another essential component of the replication complex
- **Pan-coronavirus potential** - NSP9 highly conserved, similar to my targets
- **Mechanism: oligomerization disruption** - 2NSP23 stabilizes non-functional tetrameric NSP9
- **Comprehensive validation** - RNA-seq, qPCR, multiple viral variants tested
- **Prof. Twizere is co-author** - direct connection to my supervisor's lab!
- **Immune activation finding** - suggests multi-pronged antiviral mechanisms possible

**Limitations for my project:**
- **Preprint status** - NOT peer-reviewed yet, findings preliminary
- **Nanobody approach** - different from my small-molecule screening strategy
- **mRNA-LNP delivery** - not applicable to traditional drug development
- **Limited structural details** - exact NSP9 epitope/binding site not fully described in this paper

**How I will use this:**
1. **Background on RTC composition** - NSP9 role in viral replication complex
2. **Evidence for RTC as druggable target** - validates my overall approach
3. **Potential alternative/additional target** - NSP9 could complement NSP10-NSP16 targeting
4. **Validation strategies** - RNA-seq approach for assessing antiviral effects
5. **Pan-coronavirus potential** - conservation argument applies to my targets too
6. **Discussion point** - compare different RTC targeting strategies (NSP9 vs NSP10-NSP16)
7. **Ask Prof. Twizere** - can discuss NSP9 structure, epitopes, potential collaboration

---

## Key Findings

### Main Results

1. **LNP-mRNA delivery system validated:**
   - 2NSP23 nanobody mRNA successfully delivered via lipid nanoparticles
   - Transfection efficiency peaked at 1/50 dilution
   - mRNA translated into functional nanobody within 16h in HEK293T cells
   - Batch-to-batch variability observed (S53 > S55)

2. **Antiviral activity against SARS-CoV-2:**
   - **EC50 = 1.8 µM** (icSARS-CoV-2-mNeonGreen reporter assay)
   - **CC50 >30 µM** (no cytotoxicity)
   - Therapeutic index >16
   - Comparable efficacy to Remdesivir (FDA/EMA-approved antiviral)

3. **Broad-spectrum activity across variants:**
   - **Wuhan strain:** >90% inhibition at 0.4 µg/ml
   - **Delta (B1.617.X):** >90% inhibition
   - **Mu (B1.621):** >90% inhibition  
   - **Omicron (B1.1.529):** >90% inhibition
   - **Beta (UK variant):** ~60% inhibition (slightly more resistant)
   - Validated in both VeroE6 and HEK293-ACE2 cell lines

4. **Mechanism: NSP9 tetramerization:**
   - From previous work (Esposito et al. 2021): 2NSP23 stabilizes non-functional tetrameric NSP9
   - Functional NSP9 is monomer/dimer in RTC complex
   - Tetrameric form cannot participate in RTC assembly
   - Prevents viral RNA replication

5. **RNA-seq findings - Viral infection effects:**
   - SARS-CoV-2 reads: 25-50% of total reads in infected control cells
   - **2862 genes upregulated, 2668 genes downregulated** upon infection (without nanobody)
   - Major pathways affected:
     - **Mitochondrial function SUPPRESSED** (oxidative phosphorylation, ATP synthesis, electron transport)
     - **Transcriptional regulation DYSREGULATED** (RNA Pol II, chromatin organization)

6. **RNA-seq findings - Nanobody rescue:**
   - SARS-CoV-2 reads: **near 0%** in nanobody-treated infected cells
   - **Only 219 upregulated, 241 downregulated** in infected + nanobody vs non-infected + nanobody
   - **Rescues host cell transcriptome** to near-normal levels
   - Restores mitochondrial gene expression (COX5B, NDUFS8, MRPL9, MRPL24, etc.)

7. **Immune activation in uninfected cells (CRITICAL FINDING):**
   - **460 genes upregulated, 271 downregulated** by nanobody treatment in healthy cells
   - GO terms enriched: "Defense response to virus", "Innate immune response"
   - Key genes activated BEFORE infection:
     - **RSAD2** (radical SAM domain-containing 2)
     - **OAS1, OAS2, OAS3** (2'-5'-oligoadenylate synthetases)
     - **IFIT1, IFIT2, IFIT3** (interferon-induced proteins)
     - **ISG15, MX1** (interferon-stimulated genes)
   - Suggests **prophylactic potential** - primes cells against infection

### Structural Information

**PDB IDs mentioned:**
- **NO new structures** reported in this paper
- References previous NSP9 structures (not specified)
- Cites Esposito et al. 2021 for nanobody characterization

**NSP9 in RTC context:**
- Part of mini-RTC: NSP7-2×NSP8-NSP12-2×NSP13-**NSP9**
- Located in catalytic center of RNA-dependent RNA polymerase NSP12
- Functions as **monomer or homodimer** (NOT tetramer)
- **Dimerization interface residues** coincide with binding contacts in NSP12

**Resolution:** N/A - no crystal structure of NSP9-nanobody complex provided

**Chains/proteins:**
- NSP9: Small RNA-binding protein (~12 kDa)
- Functions in RTC assembly and RNA replication
- Strong tendency to oligomerize in vitro

**Ligands co-crystallized:** N/A

### Interface Analysis

**NOTE:** This paper does NOT provide detailed interface residues. It refers to previous work (Esposito et al. 2021) for epitope mapping.

**From this paper's description:**

**NSP9 functional states:**
- **Monomer/Dimer:** Functional form for RTC assembly
- **Tetramer:** Non-functional form stabilized by 2NSP23 nanobody
- **Dimerization interface = NSP12 binding site**

**Nanobody 2NSP23 binding:**
- Binds NSP9 and stabilizes tetrameric assembly
- **Epitope characterized by NMR** (Esposito et al. 2021)
- Binding prevents functional monomer/dimer formation
- Blocks NSP9 incorporation into RTC

**Buried surface area:** Not provided

**Key interactions:**
- Not specified in this paper
- Refers to Esposito et al. 2021 for molecular dynamics simulations
- Epitope mapping revealed "epitopes on wild-type NSP9 protein"

**Hot Spots:**

**NOT explicitly defined in this paper** - would need to read Esposito et al. 2021

**Inferred critical features:**
- Residues involved in **NSP9 dimerization**
- Residues at **NSP12 binding interface**
- Nanobody likely binds at/near dimerization interface

**Mutational data:**
- None provided in this paper
- Previous studies (cited): Mutations affecting dimerization impair viral propagation

**Conservation:**
- **NSP9 highly conserved across coronaviruses**
- **"Extremely low degree of mutagenicity"** (cited: Abbasian et al 2023)
- Makes it attractive for **pan-coronavirus targeting**
- Less prone to resistance mutations than Spike protein

### Druggability Assessment

**Pockets mentioned:**
- **NO specific pockets identified** (nanobody approach, not small molecule)
- Nanobody targets **protein-protein interface** (oligomerization)
- Suggests **dimerization interface is druggable**

**Druggability score (if given):**
- Not provided (not a small-molecule study)

**Known inhibitors:**

1. **2NSP23 Nanobody (this study):**
   - **EC50 = 1.8 µM** (antiviral in cells)
   - Delivered as LNP-encapsulated mRNA
   - Molecular weight: ~15 kDa (typical nanobody)
   - Mechanism: Stabilizes non-functional NSP9 tetramer
   - No cytotoxicity at tested concentrations
   - Works against multiple SARS-CoV-2 variants

2. **2NSP90 Nanobody (mentioned):**
   - Another anti-NSP9 nanobody generated
   - Also characterized (Esposito et al. 2021)
   - Less detail provided in this paper

3. **136 unique nanobodies generated total** against NSP9
   - 2NSP23 and 2NSP90 most promising
   - Specifically recognize NSP9 in COVID-19 patient saliva samples

**Comparison to approved drugs:**
- **Comparable to Remdesivir** (FDA/EMA-approved, EC50 in similar range)
- **Better than many Spike-targeting nanobodies** (more conserved target)

**Advantages of NSP9 as target:**
- **Conserved** across variants (unlike Spike)
- **Essential** for viral replication
- **Intracellular** - less subject to antibody escape mechanisms
- **Low mutation frequency**

---

## Figures/Tables of Interest

### Figure 1: LNP Delivery and Expression

**1A - LNP composition:**
- Ionizable lipid (C12-200), PEGylated lipid, phospholipid, cholesterol
- Nano-assembly microfluidic mixing technology

**1C - Dose-response of tdTomato expression:**
- Peak expression at 1/50 LNP dilution
- Used to optimize delivery

**1E-F - Nanobody expression:**
- Immunostaining shows successful intracellular expression
- Batch S53 > S55 in efficiency
- **Action:** mRNA-LNP delivery is feasible but variable

### Figure 2: Viral Inhibition Assays (CRITICAL)

**2A - Experimental pipeline:**
- Cells seeded → nanobody added → infection → readout at 16h
- Clear workflow for antiviral testing

**2B - Dose-response curve:**
- **EC50 = 1.8 µM** for LNP-mRNA-2NSP23
- **CC50 >30 µM** (no cytotoxicity)
- **KEY BENCHMARK:** Aim for similar or better EC50 in my work

**2C - Control (tdTomato):**
- No antiviral effect
- Confirms specificity of 2NSP23

**2D - Representative images:**
- Green = viral mNeonGreen
- Clear reduction in viral replication with 2NSP23
- **Visual proof** of antiviral activity

**2E - Dual cell line validation:**
- VeroE6 and HEK293-ACE2
- Dose-dependent inhibition in both
- **Action:** Consider testing in multiple cell lines

**2F - Variant inhibition (VERY IMPORTANT):**
- **>90% inhibition:** Wuhan, Delta, Mu, Omicron
- **~60% inhibition:** Beta variant
- **Insight:** NSP9 targeting works across variants (pan-coronavirus potential!)

### Figure 3: RNA-Seq Analysis - Viral Replication

**3A - Alignment to genomes:**
- Bar graph showing % SARS-CoV-2 vs human reads
- **Infected + tdTomato:** 25-50% viral reads
- **Infected + 2NSP23:** Near 0% viral reads
- **KEY FINDING:** Nanobody blocks viral RNA replication

**3B - Hierarchical clustering:**
- All samples cluster by infection status (with tdTomato)
- But nanobody-treated samples cluster together regardless of infection
- **Insight:** Nanobody suppresses transcriptional impact of infection

**3C - Distance matrix (UK variant):**
- Clear separation between infected/non-infected (tdTomato)
- Separation blurred with nanobody treatment

**3D - PCA plot:**
- Variance between infected/non-infected smaller with nanobody
- **Action:** Could use PCA in my validation analyses

**3E-H - MA plots (CRITICAL):**
- **3E:** UKdTOM infected vs non-infected: **2862 up, 2668 down**
- **3F:** UKaNSP infected vs non-infected: **219 up, 241 down**
- **3G:** Infected: UKaNSP vs UKdTOM: **2003 up, 2023 down**
- **3H:** Non-infected: NCaNSP vs NCdTOM: **460 up, 271 down**
- **Interpretation:** Nanobody rescues ~90% of infection-induced changes

### Figure 4: Gene Ontology - Infection Effects

**4A - GO terms (all DE genes):**
- **Biological Process:** Transcription, mitochondria, respiration
- **Cellular Component:** Nucleoplasm, nucleus, mitochondrion
- **KEGG Pathway:** Oxidative phosphorylation, Huntington disease, Parkinson disease

**4B - Upregulated genes:**
- Transcriptional regulation
- Chromatin organization
- **Action:** Virus hijacks transcription machinery

**4C - Downregulated genes:**
- **Mitochondrial function HEAVILY SUPPRESSED**
- Oxidative phosphorylation
- Electron transport chain
- **Key genes:** COX5B, NDUFS8, MRPL9, MRPL24

**4D - Venn diagram:**
- 5120 genes DE in UKdTOM vs NCdTOM
- 412 genes DE in UKaNSP vs NCaNSP
- 48 genes unique to nanobody treatment
- **Overlap shows nanobody rescues majority of infection effects**

**4E-G - Heatmaps:**
- **4E:** All DE genes - clear infection signature rescued by nanobody
- **4F:** Transcription-related genes - partially rescued
- **4G:** Mitochondrial genes - **FULLY RESCUED by nanobody**
- **Action:** Mitochondrial gene expression = good biomarker for antiviral efficacy

### Figure 5: Immune Activation (CRITICAL FINDING)

**5A - GO terms (nanobody in uninfected cells):**
- **Defense response to virus**
- **Innate immune response**
- **Negative regulation of viral genome replication**
- **Interferon signaling pathways**
- **KEY:** These are activated BEFORE infection!

**5B - Upregulated genes only:**
- Same immune pathways
- Confirms activation (not suppression)

**5C - Venn diagram:**
- 246 genes unique to NCaNSP vs NCdTOM
- 485 genes common between non-infected and infected nanobody effects
- 3541 genes unique to infected comparison
- **Insight:** Nanobody has specific immune-activating signature

**5D - Heatmap of all DE genes:**
- Clear upregulation pattern with nanobody treatment
- **Prophylactic potential** suggested

**5E - Venn diagram (Defense + Innate immunity):**
- 21 genes specific to "Defense response"
- 45 genes in intersection
- 29 genes specific to "Innate immune response"

**5F-H - Gene expression heatmaps (VERY INFORMATIVE):**
- **5F:** Intersection genes (RSAD2, OAS1, OAS2, IFI16, IFIH1, etc.)
- **5G:** Defense-specific genes (IFI44L, IFNB1, GBP1, IFIT1, ISG15, etc.)
- **5H:** Innate immunity-specific genes (APOL1, TRIM21, HLA-C, SERPING1, etc.)
- **Pattern:** Strong upregulation in nanobody-treated cells (both infected and non-infected)

**Key Genes Highlighted:**
- **RSAD2** (viperin) - broad antiviral activity
- **OAS1/2/3** - 2'-5'-oligoadenylate synthetases (interferon pathway)
- **IFIT1/2/3** - interferon-induced proteins (inhibit translation)
- **ISG15** - ISGylation (antiviral protein modification)
- **MX1** - interferon-induced GTPase (blocks viral replication)

---

## Methods Relevant to My Project

### Techniques used:

1. **LNP-mRNA delivery:**
   - NanoAssemblr microfluidic mixing
   - Particle size: 70.7 nm
   - Encapsulation of in vitro transcribed mRNA
   - **NOT applicable to my small-molecule work**

2. **Viral infection assays:**
   - **icSARS-CoV-2-mNeonGreen** (fluorescent reporter)
   - **icSARS-CoV-2-nanoLuc** (luminescent reporter)
   - Multiple variants tested (Wuhan, Beta, Delta, Mu, Omicron)
   - MOI 0.1 or 0.01
   - Readout at 16-24h post-infection
   - **Action:** These reporter assays could be used for my compound validation

3. **RT-qPCR for viral load:**
   - Target: SARS-CoV-2 **E gene**
   - Primers: E_Sarbeco_F1/R2
   - Probe: E_Sarbeco_P1 (FAM-BHQ)
   - **Action:** Standard method for quantifying viral replication

4. **RNA-Seq (Deep transcriptional profiling):**
   - TRI Reagent RNA extraction
   - NEBNext Ultra II RNA Library Prep
   - NextSeq 500/550 sequencing
   - **Analysis pipeline:**
     - FastQC quality control
     - Trimmomatic trimming
     - HISAT2 alignment (to human GRCh38.81)
     - Unmapped reads → aligned to SARS-CoV-2 (NC_045512.2)
     - HTseq-count for quantification
     - NASQAR portal for DE analysis
   - **Action:** Comprehensive validation method I could use

5. **Differential expression analysis:**
   - DESeq2 or similar
   - Cutoffs: log2(FC) ≥0.5 or ≤-0.5, adj p-value <0.05
   - GO enrichment: DAVID Bioinformatics
   - **Action:** Standard bioinformatics for RNA-seq

6. **Cell viability assays:**
   - Cell Titer-Glo Luminescent Assay (Promega)
   - Measures ATP (metabolically active cells)
   - **Action:** Essential for determining cytotoxicity

7. **Immunostaining:**
   - AF488-conjugated goat anti-alpaca antibody
   - For detecting expressed nanobody
   - **NOT applicable to small molecules**

8. **Live-cell imaging:**
   - Incucyte S3 system
   - Automated fluorescence monitoring
   - **Action:** Could adapt for compound screening

### Could I use these methods?

**YES - Directly applicable:**
- ✅ **Viral reporter assays** (mNeonGreen, nanoLuc)
- ✅ **RT-qPCR** for viral E gene
- ✅ **RNA-Seq** for comprehensive transcriptional profiling
- ✅ **Cell viability** assays (Cell Titer-Glo)
- ✅ **GO enrichment** analysis
- ✅ **Incucyte live-cell imaging** (if available)

**MAYBE - With modifications:**
- ⚠️ **Multiple cell lines** (VeroE6, HEK293-ACE2) - depends on BSL-3 access
- ⚠️ **Multiple viral variants** - depends on availability

**NO - Not applicable:**
- ❌ **LNP-mRNA delivery** (specific to nanobody approach)
- ❌ **Immunostaining** for nanobody detection

**Recommendation:**
- **Use viral reporter assays** for initial compound screening
- **RT-qPCR** for quantification
- **RNA-Seq** for in-depth mechanism studies (if funding allows)
- **Collaborate with Prof. Twizere** for access to assays/expertise

---

## Questions Raised

1. **What is the exact epitope of 2NSP23 on NSP9?**
   - Paper refers to Esposito et al. 2021 for NMR mapping
   - Need to read that paper for structural details
   - Critical for understanding if small molecules could target same site

2. **Is there a crystal structure of NSP9-2NSP23 complex?**
   - Not provided in this paper
   - Would be very useful for structure-based design
   - Could ask Prof. Twizere if unpublished

3. **How does 2NSP23 stabilize NSP9 tetramers?**
   - Mechanism mentioned but not detailed
   - Binding at dimerization interface?
   - Would help understand druggable sites

4. **Can small molecules mimic nanobody binding?**
   - Nanobodies are large (~15 kDa)
   - Small molecules typically <500 Da
   - Is the epitope accessible to small molecules?

5. **Why is Beta variant less susceptible (~60% vs >90%)?**
   - Does Beta have mutations in NSP9?
   - Could inform resistance mechanisms
   - Important for designing robust inhibitors

6. **What is the stoichiometry of 2NSP23:NSP9 binding?**
   - 1:1? 2:2? 4:4?
   - Affects how many NSP9 molecules are "trapped"
   - Relevant for potency considerations

7. **Does immune activation contribute to antiviral effect?**
   - Direct mechanism (RTC disruption) vs indirect (immune priming)
   - Which is more important?
   - Could small molecules also trigger immune response?

8. **Why do uninfected cells activate immune genes with nanobody?**
   - Off-target effect?
   - NSP9-like protein in human cells?
   - Foreign protein recognition?

9. **What is the half-life of the nanobody in cells?**
   - mRNA stability?
   - Protein stability?
   - Determines dosing frequency

10. **Could NSP9 and NSP10-NSP16 be targeted simultaneously?**
    - Combination therapy potential?
    - Synergistic effects?
    - Worth exploring in my project?

---

## Action Items

### Immediate (Week 2):
- [x] Read and annotate this paper completely
- [ ] **Read Esposito et al. 2021** for NSP9 epitope details (HIGH PRIORITY)
- [ ] **Search for NSP9 crystal structures** in PDB
- [ ] **Download NSP9 structure** for visualization
- [ ] **Check NSP9 conservation** across coronaviruses
- [ ] **Ask Prof. Twizere:**
  - NSP9 structural details
  - Unpublished epitope mapping data
  - Possibility of collaboration on validation
  - Whether NSP9 should be added to my target list

### Week 3-4:
- [ ] **Compare NSP9 vs NSP10-NSP16** as targets:
  - Conservation levels
  - Structural features
  - Druggability
  - Potential for combination therapy
- [ ] **Visualize NSP9 oligomerization interface** if structure available
- [ ] **Check if NSP9 interacts with NSP10-NSP16** in RTC
- [ ] **Document NSP9 role in RTC** for thesis background

### Week 5-7:
- [ ] **Decision point:** Add NSP9 as target or focus on NSP10-NSP16?
- [ ] If adding NSP9:
  - Identify dimerization interface residues
  - Define potential binding pockets
  - Set up docking for NSP9
- [ ] If not adding NSP9:
  - Keep as backup target
  - Cite as validation of RTC-targeting approach

### For Thesis:
- [ ] **Cite in Introduction:**
  - RTC components are druggable targets
  - NSP9 example of oligomerization-dependent function
  - Pan-coronavirus conservation
- [ ] **Cite in Background:**
  - RTC composition includes NSP9
  - Oligomerization states affect function
  - Intracellular targeting strategies
- [ ] **Discussion:**
  - Compare NSP9 vs NSP10-NSP16 targeting
  - Nanobody vs small molecule approaches
  - Multi-target combination potential
  - Immune activation as secondary mechanism

---

## Related Papers to Read

### PRIORITY 1 - NSP9 structural/functional:

1. **Esposito et al. (2021)** - **MUST READ NEXT**
   - NMR analysis of 2NSP23-NSP9 interaction
   - Epitope mapping
   - Molecular dynamics simulations
   - **Adv Biol (Weinh) 5:e2101113**
   - **Action:** This is the companion paper with structural details!

2. **Littler et al. (2020)** - NSP9 crystal structure
   - **iScience 23:101258**
   - **Action:** Get NSP9 structure for visualization

3. **Zhang et al. (2020)** - NSP9 multimerization
   - Structural basis for NSP9 oligomerization
   - **Mol Biomed 1:5**

### PRIORITY 2 - NSP9 function in replication:

4. **Yan et al. (2021)** - NSP9 in extended RTC
   - **Cell 184:184-193.e10**
   - Shows NSP9 in context of complete RTC
   - **Also needed for NSP12-NSP7-NSP8 (already on list)**

5. **Miknis et al. (2009)** - NSP9 dimerization essential
   - SARS-CoV-1 NSP9
   - **J Virol 83:3007-3018**

6. **Sutton et al. (2004)** - Original NSP9 structure and function
   - **Structure 12:341-353**

7. **Egloff et al. (2004)** - NSP9 RNA-binding
   - **Proc Natl Acad Sci USA 101:3792-3796**

### PRIORITY 3 - Host cell response:

8. **Zaffagni et al. (2022)** - NSP14 effects on transcriptome
   - Shows viral proteins affect host transcription
   - **Elife 11**
   - **Action:** Compare NSP9 vs NSP14 effects

9. **Blanco-Melo et al. (2020)** - SARS-CoV-2 host response
   - **Cell 181:1036-1045.e9**

### PRIORITY 4 - Conservation analysis:

10. **Abbasian et al. (2023)** - SARS-CoV-2 mutation landscape
    - Shows NSP9 conservation
    - **J Transl Med 21:152**

---

## Notes

### Strengths of this paper:
1. **Comprehensive validation** - multiple assays, cell lines, viral variants
2. **RNA-seq depth** - detailed transcriptional profiling
3. **Dual mechanism** - direct RTC disruption + immune activation
4. **Pan-coronavirus potential** - works across variants
5. **Clear experimental design** - well-controlled experiments
6. **Prof. Twizere involved** - potential collaboration opportunity

### Limitations:
1. **PREPRINT STATUS** - not yet peer-reviewed
2. **Nanobody approach** - different from small molecules
3. **Structural details lacking** - refers to other paper
4. **No crystal structure** of NSP9-nanobody complex
5. **mRNA-LNP delivery** - complex, expensive, not traditional drug
6. **Mechanism of immune activation unclear**
7. **Clinical translation uncertain** - mRNA delivery in vivo?

### Why this matters for my work:
- **Validates RTC targeting** - proves concept works
- **NSP9 as alternative** - backup if NSP10-NSP16 difficult
- **Conservation argument** - supports pan-coronavirus approach
- **Validation strategies** - RNA-seq, reporter assays applicable
- **Prof. Twizere connection** - potential collaboration
- **Multi-target potential** - could combine NSP9 + NSP10-NSP16

### Key takeaway:
**"NSP9 is a validated, druggable target in the RTC complex. The success of nanobodies suggests small molecules targeting the same oligomerization interfaces could also work. Combined with NSP10-NSP16 targeting, multi-component RTC disruption could provide robust pan-coronavirus antiviral activity."**

### Important consideration:
**This is PROF. TWIZERE'S work - I should:**
- Discuss with him about adding NSP9 to my targets
- Ask about unpublished structural data
- Explore collaboration for validation experiments
- Consider if multi-target approach (NSP9 + NSP10-NSP16) makes sense

---

**Last Updated:** January 28, 2025
