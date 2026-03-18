import random
from typing import List, Dict, Tuple
def fix_failure_window(
    peptide: str,
    codons: List[str],
    full_seq_start_nt: int,
    full_seq_end_nt: int,
    utr_len: int,
    codon_frequencies: Dict[str, Dict[str, float]],
) -> Tuple[List[str], bool]:
    """
    Mutate codons in the CDS window overlapping [full_seq_start_nt, full_seq_end_nt)
    on the full sequence (RBS+CDS), using codon_frequencies.
    Returns: (new_codons, changed_flag)
    """
    cds_start_nt = max(0, full_seq_start_nt - utr_len)
    cds_end_nt = max(0, full_seq_end_nt - utr_len)
    # failure fully in UTR → no change
    if cds_end_nt <= 0:
        return codons, False
    start_codon = cds_start_nt // 3
    end_codon = min(len(peptide), (cds_end_nt + 2) // 3)
    new_codons = codons[:]
    for i in range(start_codon, end_codon):
        aa = peptide[i]
        codon_dict = codon_frequencies[aa]
        codon_list = list(codon_dict.keys())
        weights = [codon_dict[c] for c in codon_list]
        new_codons[i] = random.choices(codon_list, weights=weights, k=1)[0]
    return new_codons, True