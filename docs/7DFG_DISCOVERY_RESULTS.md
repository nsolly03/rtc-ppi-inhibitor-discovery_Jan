# 7DFG Discovery Results - February 7, 2026

**Status:** ✅ Discovery notebook completed successfully  
**Notebook:** `05_7DFG_interface_discovery.ipynb`  
**Structure:** NSP12-NSP7-NSP8 RdRp Holoenzyme

---

## Hot Spots Discovered

### Interface 1: NSP12-NSP7
**Primary Hot Spot:** GLU431 (NSP12) ←→ LYS2 (NSP7)
- **Type:** Likely charged interaction
- **Grid Box Center:** (136.880, 101.085, 152.896)
- **Grid Box Size:** 25 × 25 × 25 Å³

**Notes:**
- GLU (glutamate) in NSP12 interacting with LYS (lysine) in NSP7
- Negative-Positive interaction (classic salt bridge)

### Interface 2: NSP12-NSP8
**Primary Hot Spot:** ARG331 (NSP12) ←→ ASP112 (NSP8)
- **Type:** Likely charged interaction  
- **Grid Box Center:** (133.006, 168.424, 145.106)
- **Grid Box Size:** 25 × 25 × 25 Å³

**Notes:**
- ARG (arginine) in NSP12 interacting with ASP (aspartate) in NSP8
- Positive-Negative interaction (classic salt bridge)

---

## Key Questions for Comprehensive Validation

### 1. Are these really the strongest interactions?
- Need to check ALL charged pairs
- Could there be hidden hot spots?

### 2. How many total hot spots were found?
- ❓ How many charged pairs < 5 Å for NSP12-NSP7?
- ❓ How many charged pairs < 5 Å for NSP12-NSP8?

### 3. Are there charged clusters/triads?
- Like 6W4H had K76-K78-D107 triad
- Could GLU431 interact with multiple lysines?
- Could ARG331 interact with multiple aspartates?

### 4. Distance details?
- ❓ What's the exact distance for GLU431-LYS2?
- ❓ What's the exact distance for ARG331-ASP112?
- ❓ Gap to next closest interaction?

---

## Next Steps: Comprehensive Validation

Need to create `09_7DFG_comprehensive_validation.ipynb` to:

### Interface 1: NSP12-NSP7
1. Find ALL positive residues in NSP12 (LYS + ARG)
2. Find ALL negative residues in NSP7 (ASP + GLU)
3. Calculate ALL pairwise distances
4. Rank all interactions
5. Validate GLU431-LYS2 as #1
6. Check for charged clusters

### Interface 2: NSP12-NSP8
1. Find ALL positive residues in NSP12 (LYS + ARG)
2. Find ALL negative residues in NSP8 (ASP + GLU)
3. Calculate ALL pairwise distances (chains B AND G)
4. Rank all interactions
5. Validate ARG331-ASP112 as #1
6. Check for charged clusters

---

## Comparison with 6W4H

| Metric | 6W4H | 7DFG NSP12-NSP7 | 7DFG NSP12-NSP8 |
|--------|------|-----------------|-----------------|
| Primary Hot Spot | K76-D107 | GLU431-LYS2 | ARG331-ASP112 |
| Distance | 5.15 Å | ❓ | ❓ |
| Rank | #1 of 144 | ❓ | ❓ |
| Charged Triad? | Yes (K76-K78-D107) | ❓ | ❓ |
| Interface Size | 22 residues | ❓ | ❓ |

---

## Action Items

- [ ] Check notebook output for total hot spots found
- [ ] Note exact distances for primary hot spots
- [ ] Check if secondary hot spots were found
- [ ] Create comprehensive validation notebook
- [ ] Run comprehensive analysis
- [ ] Compare results with 6W4H methodology

