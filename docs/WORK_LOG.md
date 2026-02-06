# PhD Research Work Log

**Project:** RTC-PPI Inhibitor Discovery  
**Author:** Olivier Nsekuye  
**Institution:** GIGA-VIN, University of Liège  
**Supervisor:** Prof. Jean-Claude Twizere  
**Started:** January 26, 2025

---

## How to Use This Log

**Daily entries should include:**
- Date and time spent
- Tasks completed
- Problems encountered
- Solutions found
- Ideas for next steps
- References consulted

**Weekly summaries should include:**
- Major accomplishments
- Challenges faced
- Lessons learned
- Plans for next week

---

## Week 1: January 26 - February 1, 2025

### Sunday, January 26, 2025 (8 hours)

**Goals:**
- Set up complete computational pipeline for virtual screening
- Configure development environment
- Download and prepare target structures

**Completed:**
- ✅ Installed Miniconda and created rtc-chem environment
- ✅ Installed all required packages (RDKit, Meeko, BioPython, etc.)
- ✅ Configured VS Code with Python and Jupyter extensions
- ✅ Set up Git version control and GitHub repository
- ✅ Downloaded 5 PDB structures (7DFG, 6XEZ, 7EDI, 6W4H, 6W9C)
- ✅ Created receptor preparation pipeline
- ✅ Fixed 7EDI structure (added missing NSP10 chain)
- ✅ Created interactive 3D visualization notebook
- ✅ Created ligand preparation pipeline (download, filter, convert)
- ✅ Created comprehensive documentation

**Scripts Created:**
1. 01_download_pdb_structures.py - Download and clean PDB files
2. 02_prepare_receptors.py - Initial receptor prep (had issues)
3. 02b_prepare_receptors_simple.py - Working receptor prep
4. 03_download_zinc_libraries.py - ZINC15 download instructions
5. 04_filter_ligands.py - Apply drug-like filters
6. 05_prepare_ligands_for_docking.py - Convert to PDBQT
7. git_backup.sh - Auto-backup to GitHub

**Notebooks Created:**
1. 01_visualize_structures.ipynb - 3D structure visualization

**Documentation Created:**
1. docs/SETUP_GUIDE.md - Complete installation guide
2. docs/QUICK_REFERENCE.md - Command cheat sheet
3. docs/TROUBLESHOOTING.md - Issues and solutions log

**Problems Encountered:**
1. Conda not recognized after installation
   - Solution: Ran /opt/miniconda3/bin/conda init zsh and restarted terminal
   
2. AutoDock Vina not available for Apple Silicon
   - Solution: Will use CECI HPC for docking
   
3. Meeko required prody package
   - Solution: pip install prody
   
4. 7EDI visualization only showed NSP14, missing NSP10
   - Solution: Updated script to extract both chains A and B
   - Lesson: Always verify biological assembly contains all functional components
   
5. Meeko PDBQT writing API error
   - Solution: Changed to preparator.write_pdbqt_string()

**Key Decisions:**
- Use simple PDBQT conversion for local work
- Focus on PPI-optimized filters (MW 350-700 Da)
- Defer binding pocket identification until after literature review

**Files Generated:**
- 10 PDB structures (5 original + 5 cleaned)
- 3 receptor PDBQT files (7DFG, 6W4H, 7EDI)
- 3 demo ligand PDBQT files
- Complete project structure with documentation

**Next Steps:**
- Read key papers on NSP protein structures
- Study PPI interfaces in detail using notebook
- Request Enamine REAL academic license
- Schedule meeting with Prof. Twizere
- Begin literature review on NSP inhibitors

**Time Breakdown:**
- Environment setup: 2 hours
- Script development: 4 hours
- Documentation: 1.5 hours
- Troubleshooting: 0.5 hours

**References Consulted:**
- RDKit documentation: https://www.rdkit.org/docs/
- Meeko documentation: https://github.com/forlilab/Meeko
- ZINC15 database: https://zinc15.docking.org/
- BioPython tutorial: https://biopython.org/wiki/Documentation

**Notes:**
- Apple Silicon (M1/M2) has limited support for some chemistry tools
- VirtualFlow on HPC will handle production docking
- Local Mac is for development, analysis, and visualization
- GitHub: https://github.com/nsolly03/rtc-ppi-inhibitor-discovery_Jan

**Mood/Reflection:**
Excellent first day! Set up complete pipeline from structure download to ligand preparation. Encountered and solved several technical issues. Documentation is comprehensive and will be valuable for thesis. Ready to move to literature review phase.

---

### Monday, January 27, 2025

**Goals:**
- 

**Completed:**
- 

**Problems Encountered:**
- 

**Next Steps:**
- 

**Time Spent:**

**References:**
- 

**Notes:**
- 

---

## Week 2: February 2 - February 8, 2025

### [Date]

**Goals:**
- 

**Completed:**
- 

---

## Weekly Summary Template

### Week [X]: [Date Range]

**Major Accomplishments:**
- 

**Challenges Faced:**
- 

**Lessons Learned:**
- 

**Progress Toward Milestones:**
- [ ] Environment setup
- [ ] Structure preparation
- [ ] Ligand preparation
- [ ] Binding site identification
- [ ] Grid box definition
- [ ] Test docking
- [ ] HPC deployment
- [ ] Full screening

**Plans for Next Week:**
- 

**Hours Worked:** X hours

---

## Monthly Milestones

### January 2025
- [x] Environment setup complete
- [x] Receptor preparation pipeline
- [x] Ligand preparation pipeline
- [ ] Literature review
- [ ] Binding pocket identification

### February 2025
- [ ] Complete binding site analysis
- [ ] Define docking grid boxes
- [ ] Test docking locally
- [ ] Deploy to CECI HPC

### March 2025
- [ ] Run full virtual screening
- [ ] Analyze results
- [ ] Select top 500 compounds

---

## Ideas and Future Directions

### Research Ideas
- Investigate allosteric sites in addition to PPI interfaces
- Consider covalent inhibitors for NSP14 exonuclease
- Explore peptidomimetic scaffolds for NSP7-NSP8 interface

### Technical Improvements
- Automate daily backup with cron job
- Create visualization script for docking results
- Develop ML model to predict binding affinity

### Collaborations
- Dr. Christoph Gorgulla (Harvard) - VirtualFlow optimization
- Prof. Steven Ballet (VUB) - Medicinal chemistry
- Prof. Das (KU Leuven REGA) - Structural validation

---

## Important Contacts

**Supervisor:**
- Prof. Jean-Claude Twizere

**Collaborators:**
- Dr. Christoph Gorgulla (Harvard/St. Jude) - VirtualFlow
- Prof. Steven Ballet (VUB Brussels) - Medicinal chemistry
- Prof. Das (KU Leuven REGA) - Structural biology

**IT Support:**
- CECI HPC helpdesk

---

## Resources

**Key Papers:**
- [ ] Gao et al. (2020) - NSP12 structure (7DFG)
- [ ] Lin et al. (2021) - NSP10-NSP14 structure (7EDI)
- [ ] Raman et al. (2020) - NSP10-NSP16 structure (6W4H)

**Tools:**
- ZINC15: https://zinc15.docking.org/
- Enamine REAL: https://enamine.net/
- VirtualFlow: https://virtual-flow.org/
- PDB: https://www.rcsb.org/

---

**Last Updated:** January 26, 2025
## January 28, 2025 - Literature Review Session 1

### Session 4: Paper 01 Analysis (Trepte et al. 2024)

**Time:** [Your session time]

**Activities:**
- Read and comprehensively analyzed Paper 01: Trepte et al. 2024
- Created detailed annotation file (Paper_01_Trepte_2024_NSP10-NSP16_Inhibitor.md)
- Updated LITERATURE_REVIEW.md with NSP10-NSP16 findings
- Updated summary tables (Interface Residues, Hot Spots, Inhibitors)
- Created PAPERS_TRACKING.md system
- Committed to GitHub

**Key Findings:**
- NSP10-NSP16 interface is VALIDATED as druggable
- Critical hot spots: NSP10 Lys93 ↔ NSP16 Asp106 (salt bridge)
- Compound 459: Kd 12.97 µM, IC50 9.2 µM (PPI disruption), IC50 39.5 µM (antiviral)
- Virtual screening of ~350M compounds successful
- Docking box: 75.647 × 16.822 × 17.631 Å centered on Lys93
- Top hits: ~-8.5 kcal/mol docking scores

**Action Items Generated:**
- Download 6W4H structure
- Visualize Lys93-Asp106 interaction
- Extract coordinates for grid box definition
- Install PDBePISA
- Plan docking box centered on Lys93
- Consider Enamine REAL library access

**Files Created/Modified:**
- `docs/literature_review/notes/Paper_01_Trepte_2024_NSP10-NSP16_Inhibitor.md` (NEW)
- `docs/literature_review/LITERATURE_REVIEW.md` (UPDATED)
- `docs/literature_review/PAPERS_TRACKING.md` (NEW)
- `docs/WORK_LOG.md` (THIS FILE)

**Next Steps:**
- Continue with Paper 02 or proceed with immediate action items
- Download 6W4H structure
- Begin visualization of interface

**Status:** ✅ Paper 01 complete and documented

## January 28, 2025 - Literature Review Session 2

### Session 5: Paper 02 Analysis (Venit et al. 2023 - PREPRINT)

**Time:** [Your session time]

**Activities:**
- Read and comprehensively analyzed Paper 02: Venit et al. 2023 (bioRxiv preprint)
- Created detailed annotation file (Paper_02_Venit_2023_NSP9_Nanobody_Inhibitor.md)
- Updated LITERATURE_REVIEW.md with NSP9 findings
- Added NSP9 section to literature review
- Updated summary tables (Known Inhibitors)
- Updated PAPERS_TRACKING.md
- Committed to GitHub

**Key Findings:**
- NSP9 is VALIDATED druggable RTC target (EC50 1.8 µM)
- 2NSP23 nanobody stabilizes non-functional tetrameric NSP9
- Prevents RTC assembly → blocks viral replication
- Works on multiple variants (>90% inhibition Wuhan, Delta, Mu, Omicron)
- RNA-seq: Rescues host cell transcriptome, activates innate immunity
- Dual mechanism: Direct (RTC) + Indirect (immune activation)
- **Prof. Twizere is co-author** - collaboration opportunity!

**Important Notes:**
- **PREPRINT** - not yet peer-reviewed
- Different approach (nanobody vs small molecule)
- Structural details in companion paper (Esposito et al. 2021)
- NSP9 highly conserved (low mutation rate)
- Potential alternative/complementary target to NSP10-NSP16

**Action Items Generated:**
- **URGENT:** Read Esposito et al. 2021 for NSP9 epitope details
- Download NSP9 crystal structure
- Ask Prof. Twizere about:
  - NSP9 structural data
  - Adding NSP9 to target list
  - Collaboration for validation
  - Multi-target approach feasibility
- Decision point: NSP9 as primary/backup/complementary target?

**Files Created/Modified:**
- `docs/literature_review/notes/Paper_02_Venit_2023_NSP9_Nanobody_Inhibitor.md` (NEW)
- `docs/literature_review/LITERATURE_REVIEW.md` (UPDATED - added NSP9 section)
- `docs/literature_review/PAPERS_TRACKING.md` (UPDATED)
- `docs/WORK_LOG.md` (THIS FILE)

**Next Steps:**
- **Option A:** Continue with Paper 03
- **Option B:** Read Esposito et al. 2021 (companion paper for NSP9)
- **Option C:** Contact Prof. Twizere about NSP9
- **Option D:** Start working on immediate action items from Paper 01

**Status:** ✅ Paper 02 complete and documented (with preprint caveat)

## January 28, 2025 - Session 6: 6W4H Structure Analysis

**MAJOR BREAKTHROUGH:** Discovered residue numbering discrepancy

**Activities:**
1. Downloaded 6W4H structure from PDB
2. Analyzed residue numbering system
3. Discovered critical discrepancy between paper and structure
4. Identified actual hot spot residues
5. Calculated grid box coordinates
6. Generated analysis files

**Critical Finding:**
- Paper (Trepte et al. 2024) refers to K93-D106
- 6W4H PDB uses polyprotein numbering: K76-D107
- Distance: 5.15 Å (salt bridge likely)
- Grid box center: (75.883, 11.641, 10.087)

**Files Created:**
- `scripts/analyze_6W4H_corrected.py` - Main analysis
- `scripts/check_6W4H_residues.py` - Residue verification
- `scripts/find_hotspots.py` - Hot spot discovery
- `data/analysis_results/6W4H_analysis.json` - Data
- `data/analysis_results/6W4H_vina_config.txt` - Docking config
- `data/analysis_results/6W4H_interface_residues.csv` - Interface list
- `docs/6W4H_ANALYSIS_SUMMARY.md` - Complete documentation

**Interface Composition:**
- Total: 22 residues
- NSP10: 16 residues
- NSP16: 6 residues
- Features: K-D salt bridge, multiple cysteines, hydrophobic core

**Next Steps:**
- Week 3-4: Install and run fpocket
- Validate K76-D107 site is in identified pocket
- Prepare structure for docking (Week 5)

**Status:** ✅ Week 2 objectives COMPLETE - ahead of schedule!

### Important Note: Residue Numbering

**Discrepancy discovered and resolved:**
- Paper (Trepte et al.) uses individual protein numbering: K93-D106
- 6W4H PDB uses polyprotein numbering: K76-D107
- **SAME FUNCTIONAL HOT SPOT** - just different labeling system
- See `docs/RESIDUE_NUMBERING_EXPLANATION.md` for full details

**For thesis:** Always mention both numbering systems and explain the conversion.

---

## January 29, 2025 - Session 7: PepPrCLIP Discovery & Strategic Planning

**MAJOR DEVELOPMENTS:**

### New Opportunities Identified

**1. PepPrCLIP Paper (Bhat et al., Science Advances, Jan 2025)**
- Sequence-based peptide binder design using ESM-2 + CLIP
- Validates our sequence-based approach
- Works on disordered proteins (no structure needed!)
- Potential complementary approach to small molecules

**2. NVIDIA DGX Spark Possibility**
- Personal AI supercomputer ($3,999)
- 1 PFLOP FP4 performance, 128GB memory
- Would eliminate HPC dependencies
- Enable local virtual screening + PepPrCLIP experiments
- Prof. Twizere considering for the project

**3. LLM Compression Techniques**
- Prof mentioned importance of pruning/distillation/quantization
- Relevant for custom model training (Month 3+)
- Essential for optimizing ESM-2 inference
- Key skill for computational biology

### Strategic Discussions

**Hybrid Approach Proposed:**
- **Primary:** Small molecule virtual screening (proven, druggable)
- **Secondary:** PepPrCLIP peptide design (novel, validation)
- **Rationale:** Dual modalities strengthen thesis, multiple papers possible
- **Timeline:** Compatible with 6-month PhD scope

**Computational Strategy:**
- **Option A:** JURECA HPC (original plan, queue times)
- **Option B:** DGX Spark (local, instant, full control)
- **Decision:** Pending Prof. Twizere feedback

### Technical Insights

**PepPrCLIP vs Our Approach:**
- **Similarity:** Both target protein-protein interfaces
- **Similarity:** Both use ESM-2 embeddings
- **Difference:** They generate peptides, we screen molecules
- **Difference:** They need no structure, we use 6W4H
- **Complementarity:** Could validate K76-D107 from both angles

**DGX Spark Benefits:**
- No HPC queue times (5x faster iteration)
- Run PepPrCLIP locally
- Unlimited AlphaFold-Multimer runs
- Train custom coronavirus models
- Complete computational independence

### Documentation Updates

**Files Created:**
- `docs/PEPPRCLIP_RELEVANCE.md` - Analysis of paper relevance
- `docs/RESIDUE_NUMBERING_EXPLANATION.md` - K93/K76 discrepancy resolved
- `docs/6W4H_CHARGED_CLUSTER_ANALYSIS.md` - Comprehensive Lys-Asp analysis
- `notebooks/03_comprehensive_lys_asp_analysis.ipynb` - Full systematic analysis

**Papers Added to Queue:**
- Paper #29: Bhat et al. (2025) PepPrCLIP - Science Advances

### Next Steps (Pending Prof Feedback)

**Immediate (Week 3-4):**
- [ ] Clarify project scope with Prof. Twizere
- [ ] Confirm computational resources (HPC vs DGX Spark)
- [ ] Decide on hybrid approach timing
- [ ] Continue fpocket installation and analysis

**If DGX Spark Approved:**
- [ ] Learn DGX Spark setup and optimization
- [ ] Plan expanded scope (3 targets, dual modalities)
- [ ] Prepare for custom model training
- [ ] Study LLM compression techniques

**If Staying with Original Plan:**
- [ ] Continue Week 3-4: fpocket pocket identification
- [ ] Note PepPrCLIP as future direction
- [ ] Proceed with JURECA HPC timeline

### Key Decisions Made

1. ✅ **Keep small molecules as primary approach** (proven, druggable, safe)
2. ✅ **Consider PepPrCLIP as validation** (novel, comprehensive)
3. ✅ **Wait for Prof guidance** before major pivots
4. ✅ **Stay flexible** on computational resources
5. ✅ **Document everything** for thesis

### Status Summary

**Completed:**
- ✅ Week 1-2 objectives (structure analysis)
- ✅ 6W4H interface mapped
- ✅ K76-D107 hot spot validated
- ✅ Comprehensive Lys-Asp analysis (144 pairs)
- ✅ Charged triad discovered (K76-K78-D107)
- ✅ Grid box parameters defined
- ✅ Residue numbering discrepancy resolved

**In Progress:**
- 🔄 Strategic planning with Prof. Twizere
- 🔄 Computational resource determination
- 🔄 Project scope finalization

**Upcoming:**
- ⏳ Week 3-4: fpocket pocket identification
- ⏳ Literature review (PepPrCLIP + compression)
- ⏳ Project plan finalization

**Overall Progress:** ✅ ON TRACK (ahead of schedule!)

---

## January 29, 2025 - Session 8: Documentation Update

**Activities:**
- Updated WORK_LOG.md with recent developments
- Added PepPrCLIP to literature review
- Created strategic planning documents
- Consolidated all recent analyses

**Status:** Documentation current, ready for Week 3-4 work

---

## February 6, 2026 - Week 3, Day 3

### 7DFG Interface Discovery Analysis

**Completed:**
- ✅ Created complete discovery-based analysis notebook: `05_7DFG_interface_discovery.ipynb`
- ✅ Implemented systematic charged interaction screening (Lys-Asp, Arg-Asp, etc.)
- ✅ Built automated chain detection for NSP12, NSP7, NSP8
- ✅ Integrated py3Dmol visualizations
- ✅ Set up complete workflow: discovery → mapping → grid boxes → export

**Notebook Structure (10 Sections):**
1. Setup and imports
2. Structure parsing with automatic chain detection
3. 3D visualization (py3Dmol)
4A. NSP12-NSP7 hot spot discovery
4B. NSP12-NSP8 hot spot discovery
5A. NSP12-NSP7 visualization
5B. NSP12-NSP8 visualization
6A. NSP12-NSP7 interface mapping (10 Å)
6B. NSP12-NSP8 interface mapping (10 Å)
7. Docking grid box definition (25×25×25 Ų)
8. Summary table
9. Export results (CSV + JSON)
10. Conclusions

**Methodology:**
- Discovery-based approach (no prior knowledge required)
- Distance cutoffs: <5 Å for hot spots, <10 Å for interfaces
- Identifies ALL charged interactions, ranks by distance
- Exports docking-ready grid box coordinates

**Next Steps:**
- Run 05_7DFG_interface_discovery.ipynb to discover hot spots
- Analyze results and compare with 6W4H findings
- Prepare for fpocket pocket identification
- Continue with 6XEZ analysis (alternative RdRp structure)

**Files Created:**
- `notebooks/05_7DFG_interface_discovery.ipynb` (complete analysis)
- Scripts: `add_sections_5_10_fixed.py`, `add_sections_7_10.py`

---

## February 6, 2026 - Week 3, Day 3

### 7DFG Interface Discovery Analysis

**Completed:**
- ✅ Created complete discovery-based analysis notebook: `05_7DFG_interface_discovery.ipynb`
- ✅ Implemented systematic charged interaction screening (Lys-Asp, Arg-Asp, etc.)
- ✅ Built automated chain detection for NSP12, NSP7, NSP8
- ✅ Integrated py3Dmol visualizations
- ✅ Set up complete workflow: discovery → mapping → grid boxes → export

**Notebook Structure (10 Sections):**
1. Setup and imports
2. Structure parsing with automatic chain detection
3. 3D visualization (py3Dmol)
4A. NSP12-NSP7 hot spot discovery
4B. NSP12-NSP8 hot spot discovery
5A. NSP12-NSP7 visualization
5B. NSP12-NSP8 visualization
6A. NSP12-NSP7 interface mapping (10 Å)
6B. NSP12-NSP8 interface mapping (10 Å)
7. Docking grid box definition (25×25×25 Ų)
8. Summary table
9. Export results (CSV + JSON)
10. Conclusions

**Methodology:**
- Discovery-based approach (no prior knowledge required)
- Distance cutoffs: <5 Å for hot spots, <10 Å for interfaces
- Identifies ALL charged interactions, ranks by distance
- Exports docking-ready grid box coordinates

**Next Steps:**
- Run 05_7DFG_interface_discovery.ipynb to discover hot spots
- Analyze results and compare with 6W4H findings
- Prepare for fpocket pocket identification
- Continue with 6XEZ analysis (alternative RdRp structure)

**Files Created:**
- `notebooks/05_7DFG_interface_discovery.ipynb` (complete analysis)
- Scripts: `add_sections_5_10_fixed.py`, `add_sections_7_10.py`

LOG_EOF

echo "✓ Work log updated"
```

---

## Step 3: What to Expect from the Analysis

When you run the notebook, you'll see:

**Section 2 Output:**
```
Chain A: 908 residues → NSP12 (RdRp)
Chain C: 69 residues → NSP7 (cofactor)
Chain B or G: ~106-117 residues → NSP8 (cofactor)
```

**Section 4A Output (NSP12-NSP7):**
```
✓ DISCOVERED X charged interactions!

Charged Interactions (sorted by distance):
NSP12_Res  NSP12_PDB  NSP7_Res  NSP7_PDB  Distance  Type
...

🔥 STRONGEST HOT SPOT (NSP12-NSP7):
  LYSXXX (NSP12) ←→ ASPYYY (NSP7)
  Distance: X.XX Å
  Type: Salt bridge / Ionic interaction
```

**Section 7 Output (Grid Boxes):**
```
FOR AUTODOCK VINA CONFIG FILES:
# NSP12-NSP7
center_x = XXX.XXX
center_y = XXX.XXX
center_z = XXX.XXX
size_x = 25.0
size_y = 25.0
size_z = 25.0
