# RTC-PPI Inhibitor Discovery: A Research Journey

**Author:** Olivier Nsekuye  
**Institution:** GIGA-VIN Laboratory, University of Liège  
**Supervisor:** Prof. Jean-Claude Twizere  
**Funding:** ARES                                   
**Project Duration:** 2025-2028

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [The Scientific Challenge](#the-scientific-challenge)
3. [Research Strategy](#research-strategy)
4. [The Journey Begins](#the-journey-begins)
5. [Key Milestones](#key-milestones)
6. [Technical Decisions](#technical-decisions)
7. [Challenges and Solutions](#challenges-and-solutions)
8. [Daily Progress](#daily-progress)
9. [Lessons Learned](#lessons-learned)
10. [Looking Forward](#looking-forward)

---

## Project Overview

### The Big Picture

Coronaviruses have caused three major outbreaks in the 21st century: SARS (2003), MERS (2012), and COVID-19 (2019-present). These viruses share a common replication machinery - the replication-transcription complex (RTC) - yet no drugs specifically target the protein-protein interactions that hold this complex together.

**My Mission:** Discover novel small molecules that disrupt critical protein-protein interactions within the coronavirus RTC, potentially creating pan-coronavirus inhibitors effective against current and future viral threats.

### Why This Matters

**Scientific Impact:**
- First systematic targeting of coronavirus RTC protein-protein interactions
- Novel mechanism of action complementary to existing antivirals
- Potential for pan-coronavirus activity (broad spectrum)
- May prevent resistance development

**Public Health Impact:**
- Preparedness for future coronavirus pandemics
- New therapeutic options for COVID-19 and variants
- Foundation for rapid drug development in future outbreaks

**Personal Motivation:**
As a scientist from Rwanda, I witnessed how infectious diseases disproportionately affect low-resource settings. This project combines cutting-edge computational methods with practical drug development, aiming to create accessible therapies that could benefit communities worldwide, including my home country.

---

## The Scientific Challenge

### Understanding the Enemy: Coronavirus Replication

Coronaviruses are RNA viruses with the largest genomes among RNA viruses (~30 kb). Their replication depends on a complex molecular machine - the replication-transcription complex (RTC) - composed of multiple non-structural proteins (NSPs).

**Key RTC Components:**

1. **NSP12 (RNA-dependent RNA polymerase):** The core enzyme that synthesizes viral RNA
2. **NSP7 and NSP8:** Cofactors that stabilize NSP12 and enhance processivity
3. **NSP14 (Exonuclease):** Proofreading enzyme that maintains replication fidelity
4. **NSP10:** Cofactor that activates both NSP14 and NSP16
5. **NSP16 (2'-O-methyltransferase):** Modifies viral RNA to evade immune detection
6. **NSP13 (Helicase):** Unwinds RNA during replication

### The Innovation: Targeting Protein-Protein Interactions

**Traditional Approach:**
- Most antiviral drugs target enzymatic active sites
- Examples: Remdesivir (NSP12 active site), Paxlovid (NSP5 protease)
- Limitation: Viruses can mutate active sites → resistance

**Our Approach:**
- Target the **interfaces** where proteins bind each other
- Disrupting these interactions blocks RTC function
- Harder for virus to develop resistance (mutations destabilize complex)
- Novelty: No approved drugs target coronavirus PPIs

**Target Interfaces:**
1. NSP12-NSP7 interface (stabilizes polymerase)
2. NSP12-NSP8 interface (enhances processivity)
3. NSP7-NSP8 interface (cofactor interaction)
4. NSP10-NSP14 interface (activates proofreading)
5. NSP10-NSP16 interface (activates cap modification)

### Why This Is Challenging

**Scientific Challenges:**
1. PPI interfaces are large and flat (difficult to drug)
2. Need molecules larger than typical drugs (MW 350-700 Da)
3. Must screen billions of compounds efficiently
4. Validation requires specialized protein interaction assays

**Technical Challenges:**
1. Ultra-large-scale virtual screening (1.4 billion compounds)
2. Requires high-performance computing infrastructure
3. Need accurate protein structure preparation
4. Complex data analysis and hit selection

**Resource Challenges:**
1. Computational: Access to HPC clusters
2. Financial: Compound ordering costs (€500-1000 per compound)
3. Time: 4-year PhD timeline
4. Experimental: BSL-3 facility access for viral assays

---

## Research Strategy

### Overall Approach: AI-Guided Drug Discovery

**Phase 1: Computational Screening (Months 1-6)**
- Structure preparation and validation
- Virtual screening of 1.4 billion compounds
- Hit identification and clustering
- Selection of top 500 compounds per target

**Phase 2: Experimental Validation (Months 6-24)**
- Protein interaction assays (NanoBiT, GPCA, N2H)
- Dose-response curves (IC50 determination)
- Hit-to-lead optimization
- Structure-activity relationship (SAR) studies

**Phase 3: Viral Testing (Months 24-36)**
- BSL-3 testing with live SARS-CoV-2
- Testing against variants (Omicron, etc.)
- Testing pan-coronavirus activity (MERS, SARS)
- Resistance mutation mapping

**Phase 4: Optimization and Publication (Months 36-48)**
- Lead optimization guided by crystallography
- Mechanism of action studies
- Manuscript preparation
- Patent application (if applicable)

### Computational Methods

**Virtual Screening Platform:**
- **Software:** VirtualFlow (Gorgulla et al., Nature 2020)
- **Infrastructure:** CECI HPC cluster (Belgium)
- **Scale:** 1.4 billion compounds × 5 targets
- **Timeline:** 2-4 weeks per target
- **Cost:** Computational time only (~€5,000 equivalent)

**Compound Libraries:**
- **ZINC15:** 1.48 billion commercially available compounds
- **Enamine REAL:** 200 million synthesizable compounds
- **Filtering:** Drug-like properties, PPI-optimized
- **Format:** SMILES → 3D → PDBQT for docking

**Docking Protocol:**
- **Receptor preparation:** Protein structure optimization
- **Grid box definition:** Cover PPI interface (~25 Å cube)
- **Scoring function:** AutoDock Vina affinity
- **Output:** Binding poses and predicted affinity

### Experimental Validation Strategy

**Tier 1: Protein Interaction Assays** (Primary screen)
- NanoBiT (split luciferase) - high throughput
- GPCA (Gaussia luciferase) - quantitative
- N2H (mammalian two-hybrid) - cellular context
- Read-out: Luminescence decrease = PPI disruption
- Throughput: 96-well or 384-well format
- Timeline: 2-4 weeks per compound set

**Tier 2: Viral Replication Assays** (Secondary validation)
- BSL-3 live virus assays
- Multiple SARS-CoV-2 variants
- Cytotoxicity counter-screening
- Timeline: 4-6 weeks per compound set

**Tier 3: Resistance Studies** (Lead characterization)
- Serial passaging under drug pressure
- Sequencing resistant variants
- Mapping resistance mutations
- Structural interpretation

### Collaborations

**Computational:**
- Dr. Christoph Gorgulla (Harvard/St. Jude)
  - VirtualFlow optimization
  - Ultra-large-scale docking expertise

**Medicinal Chemistry:**
- Prof. Steven Ballet (VUB Brussels)
  - SAR studies
  - Chemical optimization
  - Analogue synthesis

**Structural Biology:**
- Prof. Das (KU Leuven REGA Institute)
  - Protein crystallography
  - Structural validation
  - Mechanism studies

**Mass Spectrometry:**
- Dr. Elisabetta Boeri Erba (Grenoble)
  - Native MS for PPI validation
  - Stoichiometry determination

---

## The Journey Begins

### January 26, 2025: Day 1 - Foundation Building

**Morning: Setting Up the Computational Environment**

The first challenge in any computational drug discovery project is establishing a robust, reproducible environment. Unlike wet-lab experiments, computational work requires careful software management to ensure reproducibility.

**Decision 1: Python Environment Manager**

I chose **Miniconda** over Anaconda because:
- Lighter footprint (~500 MB vs 5 GB)
- Faster environment creation
- Only installs what you need
- Standard in computational chemistry
- Better for eventual HPC deployment
```bash
# Installation command
/opt/miniconda3/bin/conda init zsh
```

**Learning Point:** On Mac, you must initialize conda for your shell (zsh) before it recognizes the `conda` command. This wasn't immediately obvious and cost 15 minutes of troubleshooting.

**Decision 2: Package Selection**

For computational drug discovery, I needed:

**Core Chemistry:**
- **RDKit:** Molecular manipulation, descriptor calculation
- **OpenBabel:** Format conversion
- **Meeko:** Ligand preparation for docking

**Structural Biology:**
- **BioPython:** PDB file handling, sequence analysis

**Data Science:**
- **Pandas, NumPy:** Data manipulation
- **SciPy, scikit-learn:** Statistical analysis, clustering

**Visualization:**
- **Matplotlib, Seaborn:** 2D plots
- **py3Dmol:** 3D molecular visualization in notebooks

**Development:**
- **JupyterLab:** Interactive analysis
- **VS Code:** Code development

**Challenge:** Some packages (like AutoDock Vina) aren't compiled for Apple Silicon (ARM64 architecture). 

**Solution:** For local work, I created a simple PDBQT converter. Production docking will run on CECI HPC where Vina is properly installed. This pragmatic decision allowed me to continue development without waiting for ARM64 support.

**Afternoon: Structure Acquisition**

**Target Selection Rationale:**

I'm focusing on five protein structures representing three critical RTC complexes:

1. **7DFG & 6XEZ:** NSP12-NSP7-NSP8 (RdRp holoenzyme)
   - Priority: **Highest**
   - Reason: Core replication machinery
   - Three potential interfaces to target

2. **6W4H:** NSP10-NSP16 (2'-O-MTase)
   - Priority: **High**
   - Reason: Immune evasion mechanism
   - Well-defined interface

3. **7EDI:** NSP10-NSP14 (ExoN proofreading)
   - Priority: **High**  
   - Reason: Resistance prevention
   - Could synergize with nucleoside analogues

4. **6W9C:** NSP13 (Helicase)
   - Priority: **Medium**
   - Reason: Alternative target
   - For comparison/control

**Script Development: 01_download_pdb_structures.py**

This script automates structure acquisition from the RCSB PDB database. Key features:
```python
# Download from RCSB
url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

# Clean structures (remove water, ions, ligands)
class ProteinCleaner(Select):
    def accept_residue(self, residue):
        return residue.id[0] == ' '  # Only standard residues

# Extract structure information
chains = list(structure.get_chains())
for chain in chains:
    print(f"Chain {chain.id}: {len(residues)} residues")
```

**Success:** Downloaded and cleaned all 5 structures in ~5 minutes.

**Evening: Critical Bug Discovery**

While creating the visualization notebook, I noticed something wrong with 7EDI:

**Expected:** Two chains (NSP14 + NSP10)  
**Observed:** Only one chain (NSP14 in purple)

**Investigation:**
```python
# Check what chains exist
structure = parser.get_structure('7EDI', 'data/targets/7EDI_clean.pdb')
for chain in structure.get_chains():
    print(f"Chain {chain.id}: {len(residues)} residues")

# Output showed:
# Chain A: 1024 residues (NSP14)
# Chain B: 139 residues (NSP10)  # ← This was missing!
```

**Root Cause:** In my receptor preparation script, I only extracted Chain A:
```python
'7EDI': {'chains': ['A']}  # Wrong - missing NSP10!
```

**Impact Analysis:**
- If left uncorrected, I would have screened against NSP14 alone
- The NSP10-NSP14 **interface** wouldn't exist in my structure
- Compounds would dock in wrong pockets
- Entire 7EDI screening would be invalid
- Months of wasted computational time
- Thousands of euros in wasted compound ordering

**Fix:**
```python
'7EDI': {'chains': ['A', 'B']}  # Correct - both proteins
```

**Lesson Learned:** Always verify biological assembly contains all functional components. This is exactly the type of error that could invalidate an entire PhD project if caught late.

**Verification Strategy Developed:**
```python
# Always run after extraction
for chain in structure.get_chains():
    residues = list(chain.get_residues())
    print(f"Chain {chain.id}: {len(residues)} residues")
```

**Reflection:** This caught error on Day 1 is a perfect example of why careful validation at each step is crucial. The 3D visualization wasn't just for pretty pictures - it revealed a fundamental structural error.

---

### Key Technical Decisions Made

**Decision 1: Local Development vs HPC**

**Local Mac (Development):**
- Structure preparation and validation
- Script development and testing
- Data analysis and visualization
- Manuscript figure generation

**CECI HPC (Production):**
- Ultra-large-scale docking (billions of molecules)
- Parallel processing (100+ cores)
- Long-running computations

**Rationale:** Separates development (iterative, interactive) from production (batch processing, high-throughput).

**Decision 2: PPI-Optimized vs Standard Filters**

Traditional drug-like filters (Lipinski's Rule of Five):
- MW ≤ 500 Da
- LogP ≤ 5
- HBD ≤ 5, HBA ≤ 10

**My PPI-optimized filters:**
- **MW: 350-700 Da** (PPIs need larger molecules)
- **LogP: 2-5** (moderate lipophilicity)
- **Rings ≥ 3** (rigidity for binding flat surfaces)
- **RotBonds ≤ 10** (not too flexible)

**Rationale:** PPI interfaces are large and flat, requiring molecules with greater surface area than typical active site inhibitors. This is based on literature analysis of successful PPI inhibitors.

**Decision 3: Version Control from Day 1**

Initialized Git repository immediately and created auto-backup script:
```bash
#!/bin/bash
git add .
git commit -m "Auto-backup: $(date)"
git push origin main
```

**Rationale:** 
- Protects against data loss
- Creates audit trail for thesis
- Enables collaboration
- Professional practice

**Outcome:** Already have complete history of development with 15+ commits documenting every major decision.

---

## Key Milestones

### Milestone 1: Environment Setup ✅ (January 26, 2025)

**Achieved:**
- Functional Python environment with all required packages
- VS Code configured for development
- Git version control operational
- Project structure established

**Deliverables:**
- `rtc-chem` conda environment
- 11 essential packages installed
- GitHub repository initialized
- Complete documentation

**Evidence:**
```bash
conda list  # Shows all 11 packages
git log --oneline  # Shows commit history
```

**Time:** 2 hours (planned: 3 hours) ✅

---

### Milestone 2: Structure Preparation ✅ (January 26, 2025)

**Achieved:**
- 5 PDB structures downloaded and cleaned
- 3 receptor complexes prepared for docking
- Critical bug discovered and fixed (7EDI chain issue)
- Interactive visualization working

**Deliverables:**
- 10 PDB files (5 original + 5 cleaned)
- 3 PDBQT receptors (7DFG, 6W4H, 7EDI)
- Jupyter notebook with 3D visualization
- Validation workflow established

**Quality Metrics:**
- All structures visualized ✅
- All chains verified ✅
- File sizes reasonable ✅
- Center of mass calculated ✅

**Time:** 4 hours (planned: 4 hours) ✅

---

### Milestone 3: Ligand Pipeline ✅ (January 26, 2025)

**Achieved:**
- Complete ligand preparation workflow
- Download instructions for ZINC15/Enamine
- Drug-like filtering implemented
- SMILES to PDBQT conversion working

**Deliverables:**
- 3 Python scripts (download, filter, convert)
- Demo compound library
- 3 prepared ligand PDBQT files
- Processing pipeline validated

**Quality Metrics:**
- Filtering logic correct ✅
- PDBQT format valid ✅
- Energy minimization working ✅
- Success rate: 100% (3/3 demo compounds) ✅

**Time:** 3 hours (planned: 4 hours) ⚡ (ahead of schedule)

---

### Milestone 4: Documentation System ✅ (January 26, 2025)

**Achieved:**
- Comprehensive setup guide
- Quick reference commands
- Troubleshooting log
- Work log template
- Project narrative (this document)

**Deliverables:**
- 5 documentation files
- Auto-backup script
- Daily log helper scripts
- Weekly summary template

**Quality Metrics:**
- Complete installation steps ✅
- All issues documented ✅
- Reproducible workflow ✅
- Thesis-ready format ✅

**Time:** 1.5 hours (planned: 2 hours) ⚡

---

### Upcoming Milestones (Next 3 Months)

**Milestone 5: Literature Review** (Week 2-3)
- [ ] Read 20 key papers on NSP structures
- [ ] Document known PPI inhibitors
- [ ] Identify druggable pockets from literature
- [ ] Create annotated bibliography

**Milestone 6: Binding Site Analysis** (Week 4-5)
- [ ] Run fpocket on all structures
- [ ] Manual inspection of interfaces
- [ ] Hot spot identification
- [ ] Select 2-3 pockets per target

**Milestone 7: Grid Box Definition** (Week 6)
- [ ] Calculate optimal grid centers
- [ ] Determine grid dimensions
- [ ] Validate with known inhibitors (if available)
- [ ] Create VirtualFlow configs

**Milestone 8: Test Docking** (Week 7-8)
- [ ] Local docking with 1,000 compounds
- [ ] Validate docking protocol
- [ ] Check for biases
- [ ] Optimize parameters

**Milestone 9: HPC Deployment** (Week 9-10)
- [ ] Transfer receptors to CECI
- [ ] Download full compound libraries
- [ ] Configure VirtualFlow
- [ ] Launch production runs

**Milestone 10: Full Screening** (Week 11-16)
- [ ] Screen 1.4B compounds × 5 targets
- [ ] Monitor progress daily
- [ ] Collect results
- [ ] Initial ranking

---

## Challenges and Solutions

### Challenge 1: Apple Silicon Compatibility

**Problem:** Many computational chemistry tools aren't compiled for ARM64 (Apple Silicon M1/M2/M3 chips).

**Specific Issues:**
- AutoDock Vina: No conda package for osx-arm64
- Some Python packages need Rosetta translation layer
- Performance penalties for x86_64 emulation

**Solutions Implemented:**
1. **Pragmatic approach:** Skip local Vina installation
2. **Workaround:** Created simple PDBQT converter for development
3. **Production solution:** Will use CECI HPC (Intel architecture)

**Lesson:** Don't let perfect be the enemy of good. Use workarounds for development, rely on proper infrastructure for production.

**Time Cost:** 1 hour (debugging installation attempts)

---

### Challenge 2: Meeko API Changes

**Problem:** Meeko's PDBQT writing API changed between versions.

**Error:** `write() argument must be str, not tuple`

**Investigation:**
```python
# Old API (didn't work):
pdbqt_string = PDBQTWriterLegacy.write_string(setup)

# New API (works):
pdbqt_string = preparator.write_pdbqt_string()
```

**Solution:** Updated to new API after consulting documentation.

**Lesson:** Software documentation may lag behind actual API. When documentation fails, experiment with the actual object methods.

**Time Cost:** 30 minutes

---

### Challenge 3: Missing Protein Components

**Problem:** 7EDI extraction only included NSP14, missing the NSP10 cofactor.

**Impact:** Would have invalidated entire 7EDI screening campaign.

**Detection Method:** 3D visualization revealed single protein instead of complex.

**Root Cause Analysis:**
```python
# Bug was here:
'7EDI': {'chains': ['A']}  # Only extracted NSP14

# Should have been:
'7EDI': {'chains': ['A', 'B']}  # NSP14 + NSP10
```

**Prevention Strategy:**
1. Always visualize structures after preparation
2. Cross-reference with PDB database
3. Verify expected number of chains
4. Check biological assembly information

**Lesson:** Visualization isn't just for presentations - it's a QC step.

**Time Saved:** Potentially months of wasted work

---

### Challenge 4: Git Remote Conflicts

**Problem:** Push rejected with "remote contains work you don't have locally"

**Cause:** Made changes on GitHub web interface while also working locally.

**Solution:**
```bash
git pull origin main --rebase  # Get remote changes
git push origin main           # Now push works
```

**Prevention:** Always `git pull` before starting work session.

**Lesson:** Version control requires discipline. Establish clear workflow.

**Time Cost:** 10 minutes

---

## Daily Progress

### January 26, 2025 - Day 1 Summary

**Hours Worked:** 8 hours

**Major Achievements:**
1. ✅ Complete computational environment (Miniconda, 11 packages, VS Code)
2. ✅ Downloaded and prepared 5 target structures
3. ✅ Created 7 Python scripts + 3 bash scripts
4. ✅ Built interactive visualization notebook
5. ✅ Established ligand preparation pipeline
6. ✅ Created comprehensive documentation system
7. ✅ Discovered and fixed critical 7EDI bug
8. ✅ Set up version control with auto-backup

**Deliverables Count:**
- Scripts: 10
- Documentation files: 5
- Notebooks: 1
- PDB structures: 10 (5 original + 5 clean)
- PDBQT receptors: 3
- PDBQT ligands: 3 (demo)
- Git commits: 15+

**Problems Solved:** 6 technical issues

**Key Learning:** Validation at every step prevents downstream disasters

**Mood:** Energized and productive

**Tomorrow's Priority:** Begin literature review on NSP structures

---

### [Date] - Day [X]

**Hours Worked:**

**Major Achievements:**
- 

**Deliverables:**
- 

**Problems Encountered:**
- 

**Solutions Found:**
- 

**Key Learning:**
- 

**Tomorrow's Priority:**
- 

---

## Lessons Learned

### Technical Lessons

**Lesson 1: Start with Structure Validation**

Don't assume PDB structures are ready to use. Always:
1. Visualize in 3D
2. Check for missing chains
3. Verify biological assembly
4. Remove water/ions/ligands
5. Check for unusual residues

**Lesson 2: Platform-Specific Issues Are Real**

Apple Silicon created unexpected compatibility issues:
- Not all scientific software supports ARM64
- Workarounds are acceptable for development
- Plan for production on standard architectures (x86_64)

**Lesson 3: API Documentation Can Be Outdated**

When documentation doesn't match behavior:
- Inspect actual object methods (`dir(object)`)
- Check GitHub issues for similar problems
- Experiment with small test cases
- Document what actually works

**Lesson 4: Visualization Catches Errors**

Creating the 3D visualization wasn't just for aesthetics:
- Revealed missing NSP10 chain
- Confirmed structure integrity
- Validated chain extraction
- Provides QC checkpoint

Investment: 1 hour  
Value: Prevented months of wasted work

---

### Project Management Lessons

**Lesson 5: Documentation from Day 1**

Starting with comprehensive documentation:
- Captures decisions while fresh
- Creates reproducible workflow
- Provides material for thesis Methods
- Helps collaborators understand approach
- Facilitates troubleshooting

**Lesson 6: Git Discipline Pays Off**

Version control isn't optional:
- Protects against data loss
- Documents decision history
- Enables experimentation
- Facilitates collaboration
- Provides audit trail

**Lesson 7: Pragmatic Over Perfect**

Don't let perfect be the enemy of progress:
- Workarounds are fine for development
- Focus on production quality where it matters
- Time boxing prevents over-optimization
- Shipping beats perfection

**Lesson 8: Early Detection Saves Time**

Bugs caught early are cheap to fix:
- 7EDI chain issue: 30 min to fix on Day 1
- If caught after screening: weeks/months lost
- Validation steps are investments, not overhead

---

### Scientific Lessons

**Lesson 9: PPI Drugging Is Different**

Protein-protein interfaces require different approach:
- Larger molecules (350-700 Da vs typical 300-500 Da)
- More rings (rigidity)
- Flat surface binding
- Can't apply standard drug rules blindly

**Lesson 10: Structure Quality Matters**

Garbage in, garbage out:
- Missing chains → wrong binding site
- Wrong protonation → wrong scoring
- Incorrect biological assembly → invalid results
- Time spent on preparation prevents downstream problems

---

## Looking Forward

### Week 2: Literature Deep Dive

**Reading List (Priority Order):**

**NSP12 (RdRp):**
1. Gao et al. (2020) Science - 7DFG structure
2. Hillen et al. (2020) Nature - Remdesivir mechanism
3. Wang et al. (2020) Cell - NSP7/NSP8 function

**NSP14 (ExoN):**
1. Lin et al. (2021) NAR - 7EDI structure
2. Ferron et al. (2018) Nature - Proofreading mechanism
3. Ma et al. (2015) PNAS - NSP10 activation

**NSP16 (MTase):**
1. Raman et al. (2020) - 6W4H structure
2. Decroly et al. (2011) PLOS Path - Cap methylation
3. Bouvet et al. (2010) PLOS Path - NSP10 role

**PPI Inhibitors (General):**
1. Scott et al. (2016) Nat Rev Drug Disc - PPI drugging strategies
2. Nero et al. (2014) Nat Rev Cancer - PPI as targets
3. Arkin et al. (2014) Nat Rev Drug Disc - Small molecules for PPIs

**Goal:** By end of Week 2, understand:
- Structure-function relationships
- Known hot spots at interfaces
- Existing inhibitor data (if any)
- Druggability assessment

---

### Weeks 3-4: Binding Site Analysis

**Computational Analysis:**
- Install and run fpocket
- Calculate interface metrics
- Identify hot spots
- Rank pockets by druggability

**Manual Analysis:**
- Visualize each interface in detail
- Map conserved residues
- Identify key interactions
- Compare across coronavirus strains

**Validation:**
- Compare predictions with literature
- Check against known inhibitors
- Discuss with Prof. Twizere
- Refine target selection

**Output:** 
- Ranked list of 2-3 pockets per target
- Pocket property tables
- Visualization figures for thesis
- Grid box coordinate proposals

---

### Weeks 5-6: Docking Protocol Development

**Test Set Creation:**
- Download 1,000 diverse compounds
- Include known active molecules (if available)
- Add decoy compounds
- Prepare all as PDBQT

**Parameter Optimization:**
- Test different grid box sizes
- Vary exhaustiveness settings
- Compare scoring functions
- Validate pose reproduction

**Quality Control:**
- Check pose clustering
- Analyze score distributions
- Identify potential biases
- Document optimal parameters

**Output:**
- Validated docking protocol
- VirtualFlow configuration files
- Benchmark results
- Standard operating procedure

---

### Months 2-3: Production Screening

**HPC Preparation:**
- Request CECI account (if not yet done)
- Transfer structures and configs
- Download full compound libraries
- Test job submission

**Screening Campaign:**
- Launch VirtualFlow jobs
- Monitor progress daily
- Collect results incrementally
- Backup data continuously

**Expected Output:**
- ~1.4 billion docking poses
- Affinity scores for all compounds
- Top 10,000 per target
- Clustered by chemical similarity

**Estimated Time:** 2-4 weeks per target

---

### Months 4-6: Hit Analysis and Selection

**Primary Analysis:**
- Rank by predicted affinity
- Cluster by chemical similarity
- Filter for favorable properties
- Remove frequent hitters

**Secondary Analysis:**
- Visual inspection of top 1,000
- Check binding poses for consistency
- Identify key interactions
- Group by scaffolds

**Selection Criteria:**
- Predicted affinity (< -8 kcal/mol)
- Ligand efficiency
- Chemical diversity
- Commercial availability
- Reasonable price (<€500/compound)

**Output:**
- Top 500 compounds per target
- Ordering list with supplier info
- Budget estimate (~€100,000)
- Backup compounds (next 500)

---

### Months 7-12: Experimental Validation

**Protein Production:**
- Clone NSP pairs into expression vectors
- Express and purify protein complexes
- Validate complex formation
- Optimize assay conditions

**PPI Assay Development:**
- Set up NanoBiT system
- Optimize signal-to-noise
- Determine dynamic range
- Validate with controls

**Primary Screening:**
- Test top 500 compounds
- 8-point dose-response curves
- Calculate IC50 values
- Counter-screen for cytotoxicity

**Expected Hit Rate:** 1-5% (5-25 validated hits per target)

---

### Year 2: Optimization and Mechanism

**Hit-to-Lead:**
- Order analogues of best hits
- Structure-activity relationships
- Improve potency and selectivity
- Optimize pharmacokinetics

**Mechanistic Studies:**
- Determine binding mode
- Co-crystallization (if possible)
- Native MS to validate PPI disruption
- Resistance mutation studies

**Viral Validation:**
- BSL-3 testing with live SARS-CoV-2
- Test against variants
- Combination with existing drugs
- Resistance selection experiments

**Output:**
- 2-5 lead compounds
- Complete SAR dataset
- Mechanism of action
- 1-2 first-author publications

---

### Years 3-4: Translation and Publication

**Lead Optimization:**
- Improve drug-like properties
- ADME profiling
- Toxicity assessment
- Formulation studies

**Broader Testing:**
- Other coronaviruses (MERS, SARS)
- Animal models (if applicable)
- Patient samples (if available)
- Resistance mapping

**Thesis Writing:**
- Introduction (coronavirus biology)
- Methods (computational pipeline)
- Results (screening and validation)
- Discussion (implications)
- Conclusions (future directions)

**Publications Target:**
- 3-4 first-author papers
- Nature Communications / J Med Chem level
- Computational methods paper
- Validation results paper
- Mechanism/SAR paper

**Patent:**
- File if commercially viable leads identified
- Before first publication
- University IP office handles process

---

### Long-Term Vision (Post-PhD)

**Scientific Impact:**
- Establish PPI targeting as viable strategy
- Provide tools for next pandemic
- Create chemical probes for RTC biology
- Contribute to drug development pipeline

**Career Development:**
- Post-doc in drug discovery (industry or academia)
- Continue in antiviral research
- Maintain Rwanda-Belgium connection
- Contribute to African scientific capacity

**Societal Impact:**
- Accessible treatments for low-resource settings
- Pandemic preparedness
- Scientific diplomacy
- Inspire next generation of African scientists

---

## Reflection and Motivation

### Why This Project Matters to Me

Growing up in Rwanda, I witnessed how infectious diseases disproportionately affect communities with limited resources. The COVID-19 pandemic amplified global health inequities, with Africa initially lacking access to vaccines and treatments.

This PhD project represents an opportunity to:
1. Develop expertise in cutting-edge computational methods
2. Contribute to global pandemic preparedness
3. Create potentially affordable therapeutic options
4. Build scientific bridges between Rwanda and Europe
5. Inspire future African scientists in computational biology

### The Power of Computational Approaches

What excites me about this project is the democratizing potential of computational drug discovery:
- Requires computers, not expensive lab equipment
- Results are reproducible and shareable
- Can be done anywhere with internet
- Levels playing field for resource-limited settings
- Enables contribution from developing countries

My work demonstrates that a scientist from Rwanda, trained in Belgium, using open-source software and public databases, can compete at the frontiers of drug discovery.

### Daily Motivation

On difficult days, I remind myself:
- Every bug fixed brings me closer to functional drugs
- Each script written could help fight the next pandemic
- My success path for future African computational scientists
- Science is a marathon, not a sprint
- Perseverance matters more than perfection

### Gratitude

I'm grateful for:
- ARES funding this work
- Prof. Twizere's mentorship and trust
- GIGA-VIN lab colleagues' support
- CECI HPC infrastructure access
- International collaborators' expertise
- Family's understanding and encouragement

---

## Appendices

### A. Project Timeline
```
Month 1-2: Environment setup, structure preparation ✅
Month 3-4: Literature review, binding site analysis
Month 5-6: Protocol development, test docking
Month 7-8: HPC deployment, production screening
Month 9-12: Hit analysis and selection
Year 2: Experimental validation, hit-to-lead
Year 3: Optimization, viral testing, publications
Year 4: Thesis writing, defense preparation
```

### B. Budget Estimate

**Computational Resources:** €5,000 (HPC time)
**Compound Ordering:** €100,000 (500 compounds × 2 rounds)
**Assay Reagents:** €20,000 (proteins, kits, consumables)
**BSL-3 Testing:** €30,000 (facility access, testing)
**Travel:** €10,000 (conferences, collaborations)
**Publication:** €5,000 (open access fees)

**Total:** ~€170,000 (covered by FRIA fellowship + lab resources)

### C. Risk Mitigation

**Risk 1: Low hit rate**
- Mitigation: Screen multiple targets, 500 compounds each
- Backup: Lower threshold, screen more compounds

**Risk 2: Assay development challenges**
- Mitigation: Use established NanoBiT platform
- Backup: Alternative assays (GPCA, N2H, HTRF)

**Risk 3: Computational resource limitations**
- Mitigation: CECI access pre-arranged
- Backup: Partner with other HPC facilities

**Risk 4: BSL-3 access issues**
- Mitigation: REGA Institute collaboration established
- Backup: Surrogate virus models, pseudovirus assays

---

**This narrative will be updated weekly as the project progresses.**

**Last Updated:** January 26, 2025  
**Next Update:** February 2, 2025
