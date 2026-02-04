# rtc-ppi-inhibitor-discovery_Jan

# Computational Drug Discovery for Coronavirus Protein-Protein Interactions

**PhD Project** | Rwanda Biomedical Centre & University of Liège  
**Student:** Olivier Nsekuye (MSc Epi, MSc Big Data Analytics)  
**Supervisor:** Prof. Jean-Claude Twizere  
**Timeline:** 6 months (January - June 2025)

---

## 🎯 Project Status

**Current Phase:** Week 3-4 - Pocket Identification  
**Overall Progress:** ✅ AHEAD OF SCHEDULE  
**Last Updated:** January 29, 2025

### Completed Milestones

- ✅ **Week 1-2:** Literature review & structure analysis
  - 6W4H structure downloaded and analyzed
  - K76-D107 hot spot identified and validated
  - Charged triad discovered (K76-K78-D107)
  - Comprehensive Lys-Asp analysis (144 pairs)
  - Grid box coordinates calculated

- ✅ **Documentation:** Comprehensive project setup
  - Git repository initialized
  - All analysis scripts created
  - Jupyter notebooks functional
  - Literature review framework established

### Current Work

- 🔄 **Strategic Planning:**
  - PepPrCLIP integration assessment
  - DGX Spark computational strategy
  - Hybrid approach evaluation
  
- ⏳ **Week 3-4:** Pocket identification (fpocket)
- ⏳ **Literature Review:** 2/29 papers complete

### Recent Developments

**🆕 January 29, 2025:**
- Discovered PepPrCLIP paper (Bhat et al., Science Advances 2025)
- Evaluating hybrid approach (small molecules + peptides)
- Considering NVIDIA DGX Spark for computational independence
- Learning LLM compression techniques for future model optimization

---

## 📊 Project Overview

### Primary Objective
Design small molecule inhibitors targeting critical coronavirus protein-protein interactions using computational methods.

### Potential Secondary Objective
Validate interface druggability using peptide binder design (PepPrCLIP) as complementary approach.

### Targets
1. **NSP10-NSP16** (2'-O-methyltransferase) - PRIMARY
2. **NSP10-NSP14** (exoribonuclease) - SECONDARY
3. **NSP12-NSP7-NSP8** (RNA polymerase) - TERTIARY

---

## 🔬 Methodology

### Current Approach
**Small Molecule Virtual Screening** (Established)
- Structure-based drug design
- AutoDock Vina docking
- Grid-based pocket targeting
- Clinically relevant modality

### Under Consideration
**PepPrCLIP Peptide Design** (Novel)
- Sequence-based binder generation
- ESM-2 + CLIP architecture
- No structural information required
- Validation of interface druggability

### Computational Resources
**Options:**
- JURECA HPC (current plan)
- NVIDIA DGX Spark (under consideration)

---

## 📁 Repository Structure
```
Project/
├── data/
│   ├── structures/pdb/          # PDB files (6W4H, etc.)
│   └── analysis_results/        # Analysis outputs
├── notebooks/
│   ├── 01_visualize_structures.ipynb
│   ├── 02_analyze_6W4H_interface.ipynb
│   └── 03_comprehensive_lys_asp_analysis.ipynb
├── scripts/
│   ├── analyze_6W4H_corrected.py
│   ├── find_all_lys_asp_pairs.py
│   └── check_6W4H_residues.py
├── docs/
│   ├── WORK_LOG.md                      # Daily work log
│   ├── PROJECT_TIMELINE.md              # Detailed timeline
│   ├── literature_review/               # Paper summaries
│   ├── 6W4H_ANALYSIS_SUMMARY.md
│   ├── 6W4H_CHARGED_CLUSTER_ANALYSIS.md
│   ├── RESIDUE_NUMBERING_EXPLANATION.md
│   └── PEPPRCLIP_RELEVANCE.md
└── README.md                             # This file
```

---

## 🎓 Key Findings

### 6W4H NSP10-NSP16 Interface

**Hot Spot Validated:**
- NSP10 K76 (PDB 4346) - NSP16 D107 (PDB 6904)
- Distance: 5.15 Å (salt bridge likely)
- Charged triad: K76-K78-D107
- Interface: 22 residues (16 NSP10, 6 NSP16)

**Grid Box Defined:**
- Center: (75.883, 11.641, 10.087)
- Size: 25 × 25 × 25 Å³
- Ready for docking

**Validation:**
- #1 strongest Lys-Asp interaction (of 144 pairs)
- No competing hot spots found
- Paper (Trepte et al. 2024) confirms druggability

---

## 📚 References

**Key Papers:**
1. Trepte et al. (2024) - NSP10-NSP16 hot spots validated
2. Bhat et al. (2025) - PepPrCLIP sequence-based design
3. Rosas-Lemus et al. (2020) - 6W4H structure

**Progress:** 2/29 papers reviewed

---

## 🚀 Next Steps

**Immediate (Week 3-4):**
- [ ] Install fpocket
- [ ] Identify binding pockets
- [ ] Validate K76-D107 pocket
- [ ] Finalize computational strategy

**Strategic:**
- [ ] Clarify project scope with supervisor
- [ ] Confirm computational resources
- [ ] Decide on hybrid approach
- [ ] Continue literature review

---

**Last Updated:** January 29, 2025  
**Status:** ✅ ON TRACK (documentation current)
