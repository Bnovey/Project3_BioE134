from typing import Dict
def load_amino_acid_codon_frequencies() -> Dict[str, Dict[str, float]]:
    """
    Returns: {amino_acid: {codon: relative_frequency}}
    """
    freq_table = {}  # type: Dict[str, Dict[str, float]]
    with open("genedesign/data/codon_usage.txt", "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            codon = parts[0].strip()
            aa = parts[1].strip()
            try:
                rel_freq = float(parts[2])
            except ValueError:
                continue
            if aa not in freq_table:
                freq_table[aa] = {}
            freq_table[aa][codon] = rel_freq
    return freq_table