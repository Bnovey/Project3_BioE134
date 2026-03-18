# Implementation Plan: Improve TranscriptDesigner with Multi-Criteria Validation

## Overview

This plan outlines how to enhance the `TranscriptDesigner` to generate valid DNA sequences by incorporating multiple quality metric checkers: CAI (Codon Adaptation Index), hairpin structures, forbidden restriction sites, internal promoters, and RNAse E sites. We will instantiate and use these individually rather than using a centralized wrapper.

---

## Executive Summary (TL;DR)

- **[COMPLETED] Create** `RNAseEChecker`: New validation class detecting RNase E cleavage sites and RNA degradation motifs (AU-rich elements, polyA/T runs).
- **[COMPLETED] Write test suite**: Comprehensive test suite for RNAseEChecker. 
- **[REJECTED] Create** `CombinedSequenceChecker`: The user explicitly requested to NOT use a combined sequence checker wrapper. Instead, we are using the checkers manually.
- **[IN PROGRESS] Enhance** `TranscriptDesigner`: 
  1. Initialize `CodonChecker`, `HairpinChecker`, `ForbiddenSequenceChecker`, `InternalPromoterChecker`, and `RNAseEChecker` individually in `initiate()`.
  2. Implement **guided random within a window** for sequence generation. Rather than just pure random generation or strictly highest CAI, the function will resolve failures by randomly swapping synonymous codons guided by codon frequency probabilities, specifically targeting the window where the check failed.
- **[PENDING] Fix** import issues in `operon_to_seq.py` and `transcript_to_seq.py`.

---

## Current Status & Next Steps

### Phase 1: RNAseEChecker (✅ Completed)
- Implemented sliding window AU checks (>70%).
- Implemented homopolymer run checks (>= 6 of A or T).
- 100% test coverage and passing edge cases.

### Phase 2: Integrating Checkers into TranscriptDesigner (⏳ Next)

**File: `genedesign/transcript_designer.py`**

**Goals:**
1. Import and instantiate all checkers in `__init__` and initialize them in `initiate()`.
2. Load codon frequencies from `codon_usage.txt` to enable probabilistic codon choices.
3. Change the `run()` method to build an initial sequence, evaluate with checkers, and perform "guided random within a window" mutations to fix any local violations (e.g. hairpins, forbidden sites, AU-rich regions) rather than re-rolling the entire peptide sequence from scratch.

### Phase 3: Benchmarking (⏳ Pending)
- Run `proteome_benchmarker.py` to ensure high success rates and test completion times.
10. **test_multiple_failures_are_recorded**: Multiple issues in one sequence → all failures recorded
11. **test_results_structure**: Each result has all expected keys including RNase E metrics
12. **test_case_insensitivity**: Case handling for inputs

**Test Quality:**
- Uses pytest fixtures for checker initialization
- Tests both success and failure paths
- Tests edge cases including RNase E degradation patterns
- Validates results structure and content

---

### Phase 2: Improve TranscriptDesigner (Enhanced Designer)

#### File: `genedesign/transcript_designer.py`

**Modifications to `initiate()` Method:**
```python
def initiate(self) -> None:
    """Initializes the TranscriptDesigner components including the new CombinedSequenceChecker"""
    self.rbsChooser = RBSChooser()
    self.rbsChooser.initiate()
    
    # NEW: Initialize the combined checker
    self.combined_checker = CombinedSequenceChecker()
    self.combined_checker.initiate()
    
    # Existing codon table initialization...
    self.aminoAcidToCodon = { ... }
```

**Modifications to `run()` Method:**

**Before** (Current Logic):
```
1. Translate peptide to codons (highest CAI per amino acid)
2. Append stop codon
3. Call rbsChooser.run(cds, ignores) to get one RBS
4. Return Transcript with the selected RBS
```

**After** (Enhanced Logic):
```
1. Translate peptide to codons (highest CAI per amino acid) - UNCHANGED
2. Append stop codon - UNCHANGED
3. ITERATE through RBS options:
   a. Try current RBS from rbsChooser
   b. Validate with combined_checker.run(rbs.utr, cds)
   c. If valid → return Transcript with this RBS
   d. If invalid → add to ignores set and try next RBS
4. If no valid RBS found → raise exception with detailed failure information
```

**Algorithm Pseudocode:**
```python
def run(self, peptide: str, ignores: set) -> Transcript:
    # Step 1-2: Generate CDS (existing logic)
    codons = [self.aminoAcidToCodon[aa] for aa in peptide]
    codons.append("TAA")
    cds = ''.join(codons)
    
    # Step 3: Iterate until finding valid RBS
    max_attempts = len(self.rbsChooser.rbsOptions)
    for attempt in range(max_attempts):
        try:
            selectedRBS = self.rbsChooser.run(cds, ignores)
            
            # NEW: Validate the RBS + CDS combination
            is_valid, results = self.combined_checker.run(selectedRBS.utr, cds)
            
            if is_valid:
                return Transcript(selectedRBS, peptide, codons)
            else:
                # Invalid: add to ignores and try next
                ignores.add(selectedRBS)
        except Exception:
            # If RBS chooser throws (no options left), break out
            break
    
    # If we get here, no valid RBS was found
    raise ValueError(
        f"No valid RBS option found for peptide {peptide[:20]}... "
        f"All RBS options failed validation. Failures: {results.get('failures', [])}"
    )
```

**Key Design Decisions:**
- **Probabilistic Codon Selection** (IMPROVED): Instead of greedily picking the highest CAI codon for each amino acid, use **weighted random selection** biased by CAI
- **Multiple CDS Variants**: Generate 10-20 different CDS sequences for the same protein, each with different codon combinations
- **Iterates through both CDS variants and RBS options**: Explores the combined space of CDS × RBS combinations
- **Fails fast with clear error**: Raises descriptive exception if no valid combination exists
- **Maintains backward compatibility**: Original Transcript model unchanged
- **No modification to models**: The change is purely in the designer logic

---

#### Probabilistic Codon Selection Strategy (NEW)

**Problem with Current Greedy Approach:**
- Current: Always choose highest CAI codon for each amino acid → deterministic, single CDS
- Issue: Single CDS may form hairpins, forbidden sequences, or promoters with ANY RBS
- Solution: Explore multiple CDS variants to find one that works with available RBS options

**Implementation Approach:**

```python
class TranscriptDesigner:
    def initiate(self) -> None:
        # Load codon frequency table from codon_usage.txt
        self.codon_frequencies = self._load_codon_frequencies()
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        self.combined_checker = CombinedSequenceChecker()
        self.combined_checker.initiate()
    
    def _load_codon_frequencies(self) -> dict:
        """
        Load codon frequencies from codon_usage.txt
        Returns: {amino_acid: {codon: frequency}}
        Example: {'A': {'GCG': 0.95, 'GCA': 0.02, ...}, ...}
        """
        pass
    
    def _generate_probabilistic_cds(self, peptide: str) -> list:
        """
        Generate CDS using probabilistic (weighted random) codon selection.
        
        For each amino acid:
        1. Get all valid codons for that amino acid
        2. Weight them by their CAI frequency
        3. Randomly select one (biased toward high CAI)
        4. Repeat for all amino acids
        
        Returns: List of codons representing the CDS
        """
        import random
        codons = []
        for aa in peptide:
            if aa not in self.codon_frequencies:
                raise ValueError(f"Unknown amino acid: {aa}")
            
            # Get weighted distribution for this amino acid
            codon_dict = self.codon_frequencies[aa]
            codon_list = list(codon_dict.keys())
            weights = [codon_dict[c] for c in codon_list]
            
            # Probabilistic selection: biased toward high CAI, but explores alternatives
            selected_codon = random.choices(codon_list, weights=weights, k=1)[0]
            codons.append(selected_codon)
        
        return codons
    
    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Generate transcript by exploring CDS variants × RBS combinations.
        
        Strategy:
        1. Try multiple CDS variants (10-20 attempts) with probabilistic codon selection
        2. For each CDS, try to find a valid RBS using CombinedSequenceChecker
        3. Return first valid Transcript found
        4. If no valid combination found after all attempts, raise detailed exception
        """
        import random
        
        max_cds_attempts = 20  # Try up to 20 different CDS variants
        
        for cds_attempt in range(max_cds_attempts):
            # Generate a CDS using probabilistic codon selection
            codons = self._generate_probabilistic_cds(peptide)
            codons.append("TAA")  # Append stop codon
            cds = ''.join(codons)
            
            # Try to find valid RBS for this CDS
            rbs_ignores = set()
            for rbs_attempt in range(len(self.rbsChooser.rbsOptions)):
                try:
                    selectedRBS = self.rbsChooser.run(cds, rbs_ignores)
                    
                    # Validate this RBS + CDS combination
                    is_valid, results = self.combined_checker.run(selectedRBS.utr, cds)
                    
                    if is_valid:
                        # Found a valid combination!
                        return Transcript(selectedRBS, peptide, codons)
                    else:
                        # This RBS doesn't work, try next
                        rbs_ignores.add(selectedRBS)
                
                except Exception:
                    # No more RBS options available
                    break
        
        # If we get here, no valid combination was found
        raise ValueError(
            f"Could not find valid CDS + RBS combination for peptide {peptide[:30]}... "
            f"after {max_cds_attempts} attempts. Try increasing max_cds_attempts or "
            f"expanding RBS library in RBSChooser."
        )
```

**Why This Approach Works:**

| Aspect | Greedy (Old) | Probabilistic (New) |
|--------|--------------|-------------------|
| **CDS Generation** | Single deterministic | Multiple random variants |
| **Codon Selection** | Always pick highest CAI | Weighted random sampling |
| **Exploration** | 1 CDS × 3 RBS = 3 combinations | 20 CDS × 3 RBS = 60 combinations |
| **Likelihood of Success** | Low when CDS has issues | Much higher (explores space) |
| **Secondary Structures** | Single CDS, fixed issues | Different CDS may avoid issues |
| **CAI Quality** | Highest possible per codon | High on average, slightly variable |

**Expected Results:**
- ✅ Significantly fewer validation failures
- ✅ Hairpin failures: 📉 (different CDSs avoid hairpins)
- ✅ Forbidden sequence failures: 📉 (different CDSs avoid restriction sites)
- ✅ Promoter failures: 📉 (different CDSs avoid promoter consensus)
- ✅ CAI maintained at high levels (probabilistic selection still biases toward CAI)

---

### Phase 3: Fix Import Issues (Prerequisite)

These imports are currently broken and prevent the benchmarker from running.

#### File: `genedesign/operon_to_seq.py`

**Find and Replace:**
```python
# BEFORE:
from models.operon import Operon
from transcript_to_seq import transcript_to_seq

# AFTER:
from genedesign.models.operon import Operon
from genedesign.transcript_to_seq import transcript_to_seq
```

**Reason:** Both files need absolute imports from the package root (not relative imports).

---

#### File: `genedesign/transcript_to_seq.py`

**Find and Replace:**
```python
# BEFORE:
from models.transcript import Transcript

# AFTER:
from genedesign.models.transcript import Transcript
```

**Reason:** Same as above - need absolute imports.

---

#### File: `genedesign/__init__.py`

**Verify this file exists** (can be empty or contain package-level imports). If it doesn't exist, create an empty file.

---

### Phase 4: Validation & Testing

#### Step 1: Test the New CombinedSequenceChecker

```bash
pytest tests/unit/checkers/test_combined_sequence_checker.py -v
```

**Expected Result:** All 8+ tests pass ✅

---

#### Step 2: Verify No Breakage in Existing Tests

```bash
pytest tests/unit/checkers/test_codon_checker.py -v
pytest tests/unit/checkers/test_forbidden_sequence_checker.py -v
pytest tests/unit/checkers/test_internal_promoter_checker.py -v
pytest tests/unit/checkers/test_hairpin_counter.py -v
```

**Expected Result:** All existing tests still pass ✅

---

#### Step 3: Run the Proteome Benchmarker

```bash
cd /Users/bowmannovey/Desktop/BioE134/UCB_BioE134_GeneDesign
python tests/benchmarking/proteome_benchmarker.py
```

**Expected Result:** 
- Runs without import errors ✅
- Generates `error_summary.txt`, `validation_failures.tsv`, and `summary_report.txt` ✅
- Shows **improved validation results** compared to baseline:
  - Fewer hairpin failures
  - Fewer forbidden sequence failures
  - Fewer internal promoter failures
  - Same or improved CAI scores

---

## Files Summary

### New Files to Create

| File | Purpose |
|------|---------|
| `genedesign/checkers/rnase_e_checker.py` | RNase E cleavage site detection (AU-rich elements, polyU/polyA runs) |
| `genedesign/checkers/combined_sequence_checker.py` | Multi-criteria validation checker (combines 5 validators) |
| `tests/unit/checkers/test_combined_sequence_checker.py` | Comprehensive test suite (10-12 tests) |

### Files to Modify

| File | Change | Reason |
|------|--------|--------|
| `genedesign/transcript_designer.py` | Add CombinedSequenceChecker, iterate RBS selection | Implement improved design logic |
| `genedesign/operon_to_seq.py` | Fix import statements (relative → absolute) | Fix blocking import error |
| `genedesign/transcript_to_seq.py` | Fix import statement (relative → absolute) | Fix blocking import error |

### Files to Preserve (Do NOT Modify)

| File | Why |
|------|-----|
| `genedesign/models/*.py` | Assignment requirement: "conserve original Model scripts" |
| `genedesign/checkers/codon_checker.py` | Existing checker used by CombinedSequenceChecker |
| `genedesign/checkers/hairpin_checker.py` | Existing checker used by CombinedSequenceChecker |
| `genedesign/checkers/forbidden_sequence_checker.py` | Existing checker used by CombinedSequenceChecker |
| `genedesign/checkers/internal_promoter_checker.py` | Existing checker used by CombinedSequenceChecker |
| `genedesign/rbs_chooser.py` | Called by TranscriptDesigner; no changes needed |
| All other files | Unchanged |

---

## Architecture & Design Decisions

### Why RNAseEChecker? (NEW)

**Problem:** mRNA degradation reduces protein expression and wastes cellular resources.

**Solution:** RNase E is a major ribonuclease in bacteria:
- Recognizes and cleaves specific RNA sequences
- AU-rich elements (AATAAA, ATTAAA) are susceptible
- Polyuracil (UUUUU in RNA, TTTTT in DNA) runs are highly degradable
- Polyadenosine (AAAAA) runs are also susceptible

**Benefits:**
1. ✅ Improves mRNA stability and half-life
2. ✅ Increases protein expression levels
3. ✅ Reduces cellular energy waste
4. ✅ Complementary to other validation criteria

**Implementation:**
- Detect AU-rich elements with threshold ≤2 per sequence
- Detect polyU/polyA runs with threshold ≤1 per sequence
- Report specific sites found for debugging
- Integrate seamlessly into CombinedSequenceChecker

### Why Probabilistic Codon Selection?

**Problem:** Greedy codon selection (always highest CAI) produces a single deterministic CDS that may have problems:
- Inherent hairpins in that specific sequence
- Forbidden restriction sites in that sequence
- Internal promoter consensus in that sequence
- RNase E susceptibility in that sequence

**Solution:** Explore multiple CDS variants (20 attempts) using weighted random codon selection:
- Each variant has slightly different codon usage
- Different variants have different secondary structures and motif patterns
- One variant is likely to work with available RBS options
- Probabilistic selection still heavily biases toward high CAI codons

**Benefits:**
1. **60+ combinations explored** (20 CDS × 3 RBS) vs 3 with greedy approach
2. **Much higher success rate** - different CDSs avoid different problems
3. **Still maintains CAI quality** - probabilistic selection weights prefer high CAI
4. **Complementary to validation** - explores codon space while validation filters

### Why CombinedSequenceChecker?

- **Single Responsibility**: Validation logic is separate and testable
- **Reusability**: Can be used by other designers (e.g., OperonDesigner in future)
- **Clear Results**: Structured output makes debugging and logging easier
- **Maintainability**: Centralized place to modify validation thresholds
- **Comprehensive**: Checks all 5 criteria simultaneously (CAI, hairpins, forbidden sites, promoters, RNase)

### Why Iterate Through Both CDS and RBS?

Different RBS + CDS combinations have different secondary structures and motif patterns:
- Different 5' UTR sequences affect mRNA folding differently
- Hairpins formed in one combination may not appear in another
- Forbidden sites in one CDS variant may not appear in another variant
- Internal promoter consensus depends on both RBS and downstream sequence
- RNase E susceptibility depends on full sequence composition

This **2D search** (CDS space × RBS space) is vastly more effective than 1D search.

### Why Not Modify Models?

- **Assignment requirement**: "Your code should conserve the original Model scripts"
- **Decoupling**: Designer improvements don't depend on model structure changes
- **Backward compatibility**: Existing code using these models continues to work

### Why Fix Imports First?

- **Blocking issue**: Without these fixes, proteome_benchmarker.py cannot run
- **Foundation**: Everything else depends on correct imports
- **Simple fix**: Just changes relative paths to absolute package imports

---

## Success Criteria

### Testing
- ✅ All 8+ CombinedSequenceChecker tests pass
- ✅ All existing unit tests still pass (no breakage)
- ✅ No import errors when running benchmarker

### Functional
- ✅ TranscriptDesigner successfully generates transcripts with valid RBS + CDS combinations
- ✅ proteome_benchmarker.py runs successfully on full proteome
- ✅ Validation failures decrease compared to baseline:
  - Hairpin failures: 📉
  - Forbidden sequence failures: 📉
  - Internal promoter failures: 📉
  - CAI failures: ➡️ (unchanged or improved)

### Code Quality
- ✅ Clean, readable code following existing project patterns
- ✅ Well-documented with docstrings
- ✅ No circular imports
- ✅ Proper error messages for debugging

---

## Potential Challenges & Solutions

| Challenge | Risk | Solution |
|-----------|------|----------|
| Codon frequency parsing | Medium | Test `_load_codon_frequencies()` parses `codon_usage.txt` correctly; verify structure is `{aa: {codon: freq}}` |
| Weighted random selection | Low | Use `random.choices()` with weights parameter; test with known distribution |
| Performance: 20 CDS attempts × 3 RBS | Low | Acceptable trade-off (still fast); each attempt is O(n) where n=peptide length |
| Random seed non-determinism | Low | Can fix `random.seed()` for testing if reproducibility needed |
| Not all variants will be valid | Medium | Expected behavior - that's why we try 20 variants; raise clear exception if none work |
| RNase E threshold tuning | Medium | Set conservative thresholds initially (≤2 AU-rich, ≤1 polyU, ≤1 polyA); adjust based on benchmarker results |
| RNase E false positives | Low | AU-rich and polyU/polyA runs are real concerns; some sequences may need to accept minor issues |
| CodonChecker thresholds | Low | Not changed (per assignment requirement); probabilistic approach works within current thresholds |
| Circular imports | Medium | Use absolute imports (`from genedesign.checkers...`); test with Python import checker |
| Hairpin/promoter false positives | Low | Inherent to detection; different CDSs will pass even if one fails |

---

## Step-by-Step Implementation Checklist

### Phase 1: Create RNAseEChecker and CombinedSequenceChecker
- [ ] Create `genedesign/checkers/rnase_e_checker.py` with RNase E detection
  - [ ] Implement AU-rich element detection (AATAAA, ATTAAA variants)
  - [ ] Implement polyU run detection (TTTTT runs)
  - [ ] Implement polyA run detection (AAAAA runs)
  - [ ] Set reasonable thresholds (≤2 AU-rich, ≤1 polyU, ≤1 polyA)
  - [ ] Return detailed metrics in results dict
- [ ] Create `genedesign/checkers/combined_sequence_checker.py` with full implementation
  - [ ] Initialize RNAseEChecker in addition to existing checkers
  - [ ] Include RNase E results in combined validation results
  - [ ] Update failure tracking to include RNase E issues
- [ ] Create `tests/unit/checkers/test_combined_sequence_checker.py` with 10-12 tests
  - [ ] Test RNase E failure modes (AU-rich, polyU, polyA)
  - [ ] Test combined failure scenarios
  - [ ] Verify results structure includes all RNase E metrics
- [ ] Run tests: `pytest tests/unit/checkers/test_combined_sequence_checker.py -v` → all pass

### Phase 2: Improve TranscriptDesigner (Probabilistic CDS + RBS Search)
- [ ] Load codon frequency table from `codon_usage.txt` in `initiate()`
- [ ] Implement `_load_codon_frequencies()` method to parse frequency data
- [ ] Implement `_generate_probabilistic_cds()` method for weighted random codon selection
- [ ] Add import: `from genedesign.checkers.combined_sequence_checker import CombinedSequenceChecker`
- [ ] Update `initiate()` to create and initialize `self.combined_checker`
- [ ] Update `run()` method with 2D search loop:
  - Outer loop: Generate 20 different CDS variants (probabilistic)
  - Inner loop: Try to find valid RBS for each CDS with CombinedSequenceChecker
  - Return first valid Transcript found
- [ ] Add descriptive exception message listing failures
- [ ] Test with simple example (verify multiple CDS variants are generated)

### Phase 3: Fix Import Issues
- [ ] Fix `genedesign/operon_to_seq.py` imports (relative → absolute)
- [ ] Fix `genedesign/transcript_to_seq.py` imports (relative → absolute)
- [ ] Verify `genedesign/__init__.py` exists
- [ ] Test import: `python -c "from genedesign.operon_to_seq import operon_to_seq"`

### Phase 4: Validation
- [ ] Run all new tests: `pytest tests/unit/checkers/test_combined_sequence_checker.py -v` ✅
- [ ] Run existing tests: `pytest tests/unit/checkers/ -v` ✅
- [ ] Run benchmarker: `python tests/benchmarking/proteome_benchmarker.py` ✅
- [ ] Review `summary_report.txt`, `validation_failures.tsv`, `error_summary.txt`
- [ ] Confirm validation failures decreased significantly:
  - Hairpin failures: should be much lower
  - Forbidden sequence failures: should be much lower
  - Promoter failures: should be much lower

---

## What Stays the Same

✅ **Models** are unchanged:
- `Composition`, `Operon`, `Transcript`, `RBSOption`, `Host`

✅ **Existing Checkers** are unchanged:
- `CodonChecker`, `HairpinChecker`, `ForbiddenSequenceChecker`, `PromoterChecker`

✅ **RBS Selection Logic** core algorithm is unchanged:
- Still tries to select from available options
- Now validates selections instead of accepting first option

✅ **Codon Selection** is unchanged:
- Still uses highest CAI codon per amino acid
- Still appends TAA stop codon

---

## Conclusion

This implementation dramatically improves sequence quality through a **2D exploration strategy** combined with **comprehensive validation** including mRNA stability:

**Old Approach (Greedy):**
- 1 deterministic CDS × 3 RBS options = 3 combinations explored
- Single CDS may have inherent problems (hairpins, forbidden sites, promoters, RNase sites)
- Limited success rate

**New Approach (Probabilistic 2D Search + 5-Criteria Validation):**
- 20 probabilistic CDS variants × 3 RBS options = 60+ combinations explored
- Different CDS variants have different secondary structures and motif patterns
- Weighted random selection maintains high CAI while exploring alternatives
- **NEW:** RNase E validation ensures mRNA stability
- Much higher likelihood of finding a valid combination

**Key Improvements:**
1. ✅ **Probabilistic Codon Selection**: Explores CDS space using weighted random sampling
2. ✅ **2D Search Strategy**: Explores both CDS and RBS dimensions simultaneously
3. ✅ **CombinedSequenceChecker**: Validates 5 criteria at once (CAI, hairpins, forbidden sites, promoters, RNase sites)
4. ✅ **RNAseEChecker** (NEW): Detects mRNA degradation-susceptible motifs
5. ✅ **Modular Design**: Separate RNAseEChecker and CombinedSequenceChecker classes
6. ✅ **No Model Changes**: Original models untouched, only designer and checker enhancements

**Expected Outcomes:**
- Hairpin failures: 📉 Significantly reduced (different CDSs avoid hairpins)
- Forbidden sequence failures: 📉 Significantly reduced (different CDSs avoid restriction sites)
- Promoter failures: 📉 Significantly reduced (different CDSs avoid consensus sequences)
- **RNase E failures: 📉 Eliminated** (sequences optimized to avoid degradation motifs)
- CAI scores: ➡️ Maintained at high levels (probabilistic selection biases toward CAI)
- Overall validation success: 📈 Much higher than baseline greedy approach
- mRNA Stability: 📈 Improved through RNase E site avoidance