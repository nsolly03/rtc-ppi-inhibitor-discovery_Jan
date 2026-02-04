# Paper #29: PepPrCLIP - De novo Peptide Binder Design

**Full Title:** De novo design of peptide binders to conformationally diverse targets with contrastive language modeling

**Authors:** Bhat, S., Palepu, K., Hong, L., et al.

**Journal:** Science Advances, Vol. 11, Issue 4, eadr8638

**Published:** January 22, 2025

**DOI:** 10.1126/sciadv.adr8638

---

## Executive Summary

**Key Innovation:** Sequence-only peptide binder design using ESM-2 embeddings and CLIP-based screening. No structural information required.

**Relevance to Our Project:** ⭐⭐⭐⭐⭐ (CRITICAL)
- Validates sequence-based approaches
- Uses same ESM-2 foundation model
- Targets protein-protein interfaces
- Works on disordered proteins

**Impact:** Opens door to hybrid approach (small molecules + peptides)

---

## Key Findings

### 1. Method: PepPrCLIP Pipeline

**Two-Step Process:**

**Generation:**
- Sample natural peptides from PDB training set
- Embed using ESM-2-650M
- Apply Gaussian noise to embeddings
- Decode to novel peptide sequences
- Generate ~1000 peptides per minute

**Screening:**
- CLIP-based discriminator (contrastive learning)
- Trained on ~11,597 peptide-protein pairs
- Binary accuracy: 95.4%
- Top 10% accuracy: 82%
- Ranks peptides by binding probability

### 2. Validation Results

**In Silico (AlphaFold-Multimer):**
- Hit rate: 17.7% (peptides with ipTM ≥ test set)
- RFDiffusion: 30.1% (but requires structure!)
- 22% achieve ipTM ≥ 0.7 (stable binding)
- Competitive without using structural information

**Experimental (UltraID enzyme):**
- 19/20 PepPrCLIP peptides showed inhibition
- 4/20 showed >75% inhibition
- Outperformed RFDiffusion on average (P=0.0364)
- Demonstrates functional binding + inhibition

**β-catenin (transcription factor):**
- 4/6 peptides reduced transcriptional activity
- β-cat_PpC_3 and β-cat_PpC_4: >50% degradation
- Binding affinity: KD = 150-200 nM (mid-nanomolar!)
- Proteasome-mediated degradation confirmed

**SS18-SSX1 (disordered fusion oncoprotein):**
- Successfully designed binders for highly disordered target
- SS_PpC_4: >40% reduction in fusion protein levels
- Proves method works on "undruggable" targets

### 3. Key Technical Details

**Model Architecture:**
- Encoders: ESM-2-650M embeddings
- CLIP: Contrastive learning (image-caption analogy)
- Training: Noisy dataset (11,597 pairs)
- Testing: Strict dataset (1,002 pairs)
- No structural information used

**Generation Parameters:**
- Noise scaling factor (k): 5-22
- Source: Natural PDB peptides
- Output: 15-25 amino acid peptides
- Speed: ~1000 peptides/min on T4 TPU

**Performance Metrics:**
- Binary accuracy: 95.4%
- Top 1 accuracy: 53%
- Top 10% accuracy: 82%
- CLIP score distribution: Clear separation

---

## Relevance to NSP10-NSP16 Project

### Direct Applications

**1. Interface Validation:**
- Could design peptides targeting K76-D107
- Validate interface as druggable from orthogonal angle
- Compare peptide vs small molecule efficacy

**2. Complementary Approach:**
- Small molecules: Primary (proven, druggable)
- PepPrCLIP peptides: Secondary (validation)
- Dual modalities strengthen thesis

**3. Disordered Target Capability:**
- NSP10-NSP16 has some disorder
- Future targets (fusion oncoproteins) highly disordered
- Method generalizes beyond structured proteins

### Strategic Considerations

**Advantages for Our Project:**
- ✅ Uses ESM-2 (same foundation model)
- ✅ Sequence-only (no structure bottleneck)
- ✅ Fast generation (minutes, not hours)
- ✅ Experimentally validated approach
- ✅ Works on PPIs (our exact use case)

**Challenges:**
- ❌ Peptides less druggable than small molecules
- ❌ Requires DGX Spark or similar compute
- ❌ New method (less established)
- ❌ Delivery challenges (cell permeability)

**Hybrid Strategy:**
```
Phase 1: Small molecule screening (proven)
Phase 2: PepPrCLIP validation (novel)
Result: Comprehensive dual-modality study
Output: 2-3 papers instead of 1
```

---

## Technical Comparisons

### PepPrCLIP vs RFDiffusion

| Feature | PepPrCLIP | RFDiffusion |
|---------|-----------|-------------|
| Input | Sequence only | 3D structure |
| Speed | ~1000 peptides/min | ~0.5 peptides/min |
| Hit rate | 17.7% | 30.1% |
| Training data | PDB peptide-protein | PDB structures |
| Disorder tolerance | High | Low |
| Compute | T4 TPU sufficient | A100 GPU needed |

**Verdict:** RFDiffusion better for structured targets, PepPrCLIP better for disordered

### PepPrCLIP vs Our Approach

| Aspect | PepPrCLIP | Our Virtual Screening |
|--------|-----------|----------------------|
| Target type | PPI interface | Binding pocket |
| Modality | Peptides (15-25 aa) | Small molecules |
| Design | Generate new | Screen existing |
| Druggability | Low (injection) | High (oral) |
| Specificity | Very high | Moderate |
| Speed | Fast (min) | Fast (hours) |
| Clinical path | Difficult | Established |

**Verdict:** Complementary, not competitive

---

## Key Methods to Adopt

### 1. ESM-2 Embedding Usage
```python
# Their approach
from transformers import AutoModel

model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D")
embeddings = model(sequences)

# We already use this for analysis
# Could use for custom model training
```

### 2. CLIP Architecture for Screening
```python
# Contrastive learning for binding prediction
# Could train custom CLIP for:
# - NSP10-NSP16 specific prediction
# - Small molecule binding prediction
# - Interface scoring
```

### 3. AlphaFold-Multimer Validation
```python
# They validate with ipTM scores
# We should do same for:
# - Top docking hits
# - Predicted complexes
# - Binding pose validation
```

---

## Integration into Our Project

### Immediate Actions

**Week 3-4:**
- ✅ Add to literature review (this document)
- ✅ Cite in introduction (validates approach)
- ✅ Note as complementary method

**Month 3-4 (If resources allow):**
- 🔄 Install PepPrCLIP (if DGX Spark available)
- 🔄 Generate NSP10-NSP16 targeting peptides
- 🔄 Compare with small molecule hits

**Thesis Integration:**
- Methods: Acknowledge PepPrCLIP as alternative
- Discussion: Compare modalities
- Future work: Peptide-based disruption

### Potential Experiments

**If Pursuing Hybrid Approach:**
```
Experiment 1: Small Molecule Screening
├── AutoDock Vina
├── Target: K76-D107 pocket
└── Output: Top 10 small molecule hits

Experiment 2: PepPrCLIP Generation
├── Input: NSP10-NSP16 sequences
├── Generate: 10,000 peptides
├── Screen: CLIP model
└── Output: Top 10 peptide hits

Experiment 3: Comparative Validation
├── AlphaFold-Multimer for both
├── Compare ipTM scores
├── Experimental testing (if possible)
└── Analysis: Which modality more effective?

Result: Comprehensive dual-modality thesis
```

---

## Citations and References

**Key Papers to Cite:**

1. **This paper (PepPrCLIP):**
   - For: Sequence-based binder design
   - For: ESM-2 usage validation
   - For: PPI targeting methods

2. **ESM-2 (Lin et al., 2023):**
   - For: Protein language model foundation
   - For: Embedding generation
   - For: Evolutionary-scale predictions

3. **CLIP (Radford et al., 2021):**
   - For: Contrastive learning architecture
   - For: Multi-modal embeddings

4. **AlphaFold-Multimer (Evans et al., 2022):**
   - For: Complex prediction
   - For: Validation methodology
   - For: ipTM scoring

**Citation Format:**
```
Bhat, S., et al. (2025). De novo design of peptide binders to 
conformationally diverse targets with contrastive language modeling. 
Science Advances, 11(4), eadr8638. doi:10.1126/sciadv.adr8638
```

---

## Questions for Further Investigation

1. **Sequence requirements:**
   - What sequence features make good PepPrCLIP targets?
   - Does interface size affect success rate?
   - Optimal peptide length for our target?

2. **Comparison with docking:**
   - Can we benchmark both approaches on same target?
   - Which gives better experimental validation?
   - Cost-benefit analysis?

3. **Custom training:**
   - Could we train coronavirus-specific CLIP?
   - Would fine-tuning improve performance?
   - Data requirements?

4. **Experimental validation:**
   - Cell permeability of designed peptides?
   - Stability in biological conditions?
   - Comparison with small molecules?

---

## Conclusions

**Impact on Project:**
- ⭐⭐⭐⭐⭐ CRITICAL - Changes strategic options
- Validates sequence-based approach
- Opens dual-modality possibility
- Strengthens thesis potential

**Recommendation:**
- Primary: Continue small molecule screening (proven)
- Secondary: Consider PepPrCLIP validation (novel)
- Strategy: Hybrid approach if resources permit

**Status:** 
- Added to project consideration
- Awaiting Prof. Twizere guidance
- Ready to implement if approved

---

**Reviewed:** January 29, 2025  
**Reviewer:** Olivier Nsekuye  
**Status:** Active consideration for project integration
