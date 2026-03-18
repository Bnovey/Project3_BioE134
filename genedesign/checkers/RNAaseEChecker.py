class RNAseEChecker:
    """
    Source: Copilot chat
    Detects RNase E susceptibility using two checks:
    - High AU content in a sliding window
    - Poly-A or poly-T runs
    """

    def __init__(self, window_size: int = 10, au_threshold: float = 0.70, run_length: int = 6):
        """
        Args:
            window_size: Size of the sliding window for AU check (default 10)
            au_threshold: AU content fraction to fail on (default 0.70 = 70%)
            run_length: Minimum consecutive A or T run length to fail on (default 6)
        """
        self.window_size = window_size
        self.au_threshold = au_threshold
        self.run_length = run_length

    def run(self, dna_seq: str) -> tuple:
        """
        Check sequence for AU-rich windows and poly-A/T runs.

        Args:
            dna_seq: DNA sequence (string)

        Returns:
            tuple[bool, str]: (is_valid, failure_message)
            - is_valid: True if both pass
            - failure_message: description of first failure, or None if valid
        """
        seq = dna_seq.upper()

        current_base = ''
        current_run = 0
        for i, base in enumerate(seq):
            if base == current_base and base in ('A', 'T'):
                current_run += 1
                if current_run >= self.run_length:
                    start = i - current_run + 1
                    return False, f"Poly-{base} run of {current_run} bases at position {start}"
            else:
                current_base = base
                current_run = 1

        for i in range(len(seq) - self.window_size + 1):
            window = seq[i:i + self.window_size]
            au_fraction = (window.count('A') + window.count('T')) / self.window_size
            # Only fail when AU content exceeds (not equals) the threshold.
            if au_fraction > self.au_threshold:
                return False, f"AU content {au_fraction:.0%} exceeds {self.au_threshold:.0%} threshold at position {i}"

        return True, None


if __name__ == "__main__":
    checker = RNAseEChecker()

    tests = [
        ("ATGCGCATGCGCATGCGC", True,  "Balanced GC-rich sequence"),
        ("ATTTTTTATTTTTTATTT", False, "PolyT homopolymer run"),
        ("AAAAAACGTAAAAAACGT", False, "PolyA homopolymer run"),
        ("ATATATATATCGTACGTA", False, "High AU content window"),
        ("ATGCGTACGTAGCTAGCT", True,  "Normal coding sequence"),
    ]

    for seq, expected, label in tests:
        valid, msg = checker.run(seq)
        print(valid, msg)
