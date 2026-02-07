# Notebook Run Checklist

Use this before running ANY analysis notebook to avoid errors.

---

## Pre-Run Checklist

### 1. Environment Setup ✓
- [ ] Conda environment activated: `conda activate rtc-chem`
- [ ] In correct directory: `cd ~/Desktop/Botanique/Project`
- [ ] Jupyter / VS Code ready

### 2. File Verification ✓
- [ ] PDB file exists: `ls data/structures/pdb/[STRUCTURE].pdb`
- [ ] Results directory exists: `ls data/analysis_results/`
- [ ] No path issues (check paths in Section 1)

### 3. Dependencies Check ✓
- [ ] py3Dmol installed: `pip list | grep py3Dmol`
- [ ] BioPython installed: `pip list | grep biopython`
- [ ] All imports work (run Section 1 first)

### 4. Kernel Selection ✓
- [ ] VS Code: Select **rtc-chem** kernel (top right)
- [ ] Jupyter: Kernel → Change Kernel → rtc-chem

---

## Running the Notebook

### Sequential Execution
1. **Run Section 1** (Setup) → Verify no errors  
2. **Run Section 2** (Structure parsing) → Check chain assignments  
3. **Run Section 3** (Visualization) → Verify 3D structure loads  
4. **Run Sections 4+** (Analysis) → Watch for hot spots  

---

## Common Issues & Fixes

**Issue:** PDB file not found  
**Fix:** Check path in Section 1  
Verify: `ls data/structures/pdb/7DFG.pdb`

**Issue:** ModuleNotFoundError: py3Dmol  
**Fix:** `pip install py3Dmol --break-system-packages`

**Issue:** Chain X not found  
**Fix:** Check Section 2 output  
Update `chain_assignments` if needed

**Issue:** Empty DataFrame in Section 4  
**Meaning:** No charged interactions < 5 Å  
**Action:** Verify chains are correct

---

## Post-Run Actions

### 1. Review Results ✓
- [ ] Count hot spots
- [ ] Identify strongest interaction
- [ ] Check K76 involvement (NSP10)
- [ ] Review interface size

### 2. Save Outputs ✓
- [ ] CSV files in `data/analysis_results/`
- [ ] JSON metadata saved
- [ ] Summary table exported

### 3. Document Findings ✓
- [ ] Update `docs/RESULTS_SUMMARY.md`
- [ ] Update `docs/WORK_LOG.md`
- [ ] Note key findings

---

## Structure-Specific Notes

### 7DFG (NSP12-NSP7-NSP8)
- Two interfaces: NSP12–NSP7, NSP12–NSP8
- Longer runtime due to NSP12 size

### 7N0C (NSP10-NSP14)
- Single interface
- Check K76 conservation vs 6W4H

### 6ZSL (NSP13 Helicase)
- Homodimer interface
- Expect symmetric interactions

---

## ✅ Checklist Complete
