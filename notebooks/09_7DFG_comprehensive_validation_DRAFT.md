# 7DFG Comprehensive Validation - Planning

Based on discovery results:
- Interface 1: GLU431-LYS2 (NSP12-NSP7)
- Interface 2: ARG331-ASP112 (NSP12-NSP8)

## Analysis Plan

### Part A: NSP12-NSP7 Exhaustive Analysis

**Find ALL charged residues:**
- NSP12 Chain A: ALL LYS + ARG (~100-150 residues expected)
- NSP7 Chain C: ALL ASP + GLU (~10-20 residues expected)
- Total combinations: ~1000-3000 pairs

**Key Questions:**
1. Is GLU431-LYS2 really #1?
2. Are there other GLU431-LYS interactions (cluster)?
3. Any competing hot spots?

### Part B: NSP12-NSP8 Exhaustive Analysis

**Find ALL charged residues:**
- NSP12 Chain A: ALL LYS + ARG (same as above)
- NSP8 Chain B: ALL ASP + GLU (~15-25 residues expected)
- NSP8 Chain G: ALL ASP + GLU (~15-25 residues expected)
- Total combinations: ~2000-4000 pairs (two NSP8 chains!)

**Key Questions:**
1. Is ARG331-ASP112 really #1?
2. Are there other ARG331-ASP interactions (cluster)?
3. Do chains B and G have different hot spots?

### Outputs Needed

For each interface:
- [ ] Complete distance matrix CSV
- [ ] Distance distribution histogram
- [ ] Heatmap (residue × residue)
- [ ] Top 10 interactions table
- [ ] Charged cluster identification
- [ ] 3D visualization of clusters
- [ ] Strategic implications

