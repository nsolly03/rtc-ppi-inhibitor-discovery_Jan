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
