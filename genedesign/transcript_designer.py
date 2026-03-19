import random
from pathlib import Path

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.RNAaseEChecker import RNAseEChecker


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence using the
    GeneOptimizer sliding-window algorithm, then selects an RBS.

    Algorithm (per the GeneOptimizer paper, PMID 21189842):
      1. Break the sequence into 3-aa windows (9 nt each), scan left to right.
      2. For each window, enumerate all synonymous codon combinations and score
         each in context of: upstream preamble (3 committed codons) + current
         3 codons + next 6 aa's (downstream lookahead from seed).
      3. Commit the best-scoring 3 codons and slide forward.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenSequenceChecker = None
        self.promoterChecker = None
        self.rnaseEchecker = None
        self.codonFrequencies = {}
        self.codonChoices = {}
        self.commitCodonsPerStep = 3
        self.downstreamLookaheadCodons = 6
        self.preambleCodons = 8  # 24 nt upstream context — keeps window ≥ 50 nt for hairpin checker
        self.maxTrialsPerWindow = 100

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()
        self.forbiddenSequenceChecker = ForbiddenSequenceChecker()
        self.forbiddenSequenceChecker.initiate()
        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()
        self.codonFrequencies = self._load_amino_acid_codon_frequencies()
        self.rnaseEchecker = RNAseEChecker()
        self.codonChoices = {
            aa: (tuple(freqs.keys()), tuple(freqs.values()))
            for aa, freqs in self.codonFrequencies.items()
        }


        # Highest-CAI codon for each amino acid (E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

    def _initialize_codons_highest_cai(self, peptide: str) -> list[str]:
        return [self.aminoAcidToCodon[aa] for aa in peptide]

    def _load_amino_acid_codon_frequencies(self) -> dict[str, dict[str, float]]:
        """
        Returns a codon frequency table in the form:
        {amino_acid: {codon: relative_frequency}}
        """
        codon_usage_path = Path(__file__).resolve().parent / "data" / "codon_usage.txt"
        freq_table: dict[str, dict[str, float]] = {}

        with codon_usage_path.open("r") as handle:
            for line in handle:
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

    def _score_window(self, codons: list[str], cutoff: float | None = None) -> float:
        """
        Scores a codon window. Lower is better (0.0 = fully passing).
        """
        local_dna = ''.join(codons)
        score = 0.0

        passes, _, _, _ = self.codonChecker.run(codons)
        if not passes:
            score += 1.0

        passes, _ = self.rnaseEchecker.run(local_dna)
        if not passes:
            score += 1.0

        passes, _ = self.forbiddenSequenceChecker.run(local_dna)
        if not passes:
            score += 1.0

        passes, _ = self.promoterChecker.run(local_dna)
        if not passes:
            score += 1.0


        passes, _ = hairpin_checker(local_dna)
        if not passes:
            score += 1.0

        return score

    def _optimize_tiled_windows(self, peptide: str, seed_codons: list[str], utr: str) -> list[str]:
        """
        GeneOptimizer-style rolling window optimization:
          - Prepends UTR (split into 3-nt chunks) so early windows can look back into it.
          - Iterates over CDS codons only; preamble reaches into UTR chunks when near the start.
          - Commits the best 3-codon block per window and advances by 3.
        """
        # Split UTR into 3-nt chunks (trim front so length is a multiple of 3)
        utr_trimmed = utr.upper()[len(utr) % 3:]
        utr_chunks = [utr_trimmed[i:i+3] for i in range(0, len(utr_trimmed), 3)]
        utr_offset = len(utr_chunks)

        # Full working sequence: UTR chunks + CDS seed codons
        final_codons = utr_chunks + seed_codons.copy()
        n = len(peptide)
        step = self.commitCodonsPerStep

        for i in range(utr_offset, utr_offset + n, step):
            commit_end = min(utr_offset + n, i + step)
            lookahead_end = min(utr_offset + n, commit_end + self.downstreamLookaheadCodons)
            preamble_start = max(0, i - self.preambleCodons)

            preamble = final_codons[preamble_start:i]
            lookahead = final_codons[commit_end:lookahead_end]
            amino_commit = peptide[i - utr_offset:commit_end - utr_offset]

            best_commit = list(final_codons[i:commit_end])
            best_score = self._score_window(preamble + best_commit + lookahead)

            for _ in range(self.maxTrialsPerWindow):
                if best_score == 0.0:
                    break
                trial_commit = []
                for aa in amino_commit:
                    choices, weights = self.codonChoices[aa]
                    trial_commit.append(random.choices(choices, weights=weights, k=1)[0])

                trial_score = self._score_window(
                    preamble + trial_commit + lookahead,
                    cutoff=best_score
                )
                if trial_score < best_score:
                    best_score = trial_score
                    best_commit = trial_commit

            final_codons[i:commit_end] = best_commit

        # Return only the CDS part (strip UTR chunks)
        return final_codons[utr_offset:]

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript with the selected RBS and codon list.
        """
        seed_codons = self._initialize_codons_highest_cai(peptide)
        initialRBS = self.rbsChooser.run(seed_codons, ignores)
        utr = initialRBS.utr

        cds = self._optimize_tiled_windows(peptide, seed_codons, utr)
        cds.append("TAA")

        selectedRBS = self.rbsChooser.run(cds, ignores)
        return Transcript(selectedRBS, peptide, cds)


def main() -> None:
    """Run a small example transcript design and print key outputs."""
    peptide = "MNPSDVFQIIEGHTKLMRDSIPLIASENLTSLSVRRCYVSDLGHRYAEGRVGERFYEGCKYVDQIESMAIELTRKIFEAEHANVQPISGVVANLAAFFALTNVGDTIMSISVPCGGHISHDRVSAAGLRGLRVIHYPFNSEEMSVDVDETRKVAERERPKLFILGSTLILFRQPVKEIREIADEIGAYVMYDASHVLGLIAGKAFQNPLKEGADVMTGSTHKTFFGPQRAIIASRKELAEKVDRAVFPGVVSNHHLNTLAGYVVAAMEMLEFGEDYAKQVVRNAKALAEELYSLGYKVLGEKRGFTETHQVAVDVREFGGGERVAKVLENAGIILNKNLLPWDSLEKTANPSGIRIGVQEVTRIGMKEEEMRAIAEIMDAAIKEKKSVDELRNEVKELKERFNVIKYSFDESEAYHFPDLR"

    designer = TranscriptDesigner()
    designer.initiate()
    transcript = designer.run(peptide, set())

    cds = "".join(transcript.codons)
    full_seq = transcript.rbs.utr.upper() + cds

    print("Peptide:", peptide)
    print("Selected RBS:", transcript.rbs.gene_name)
    print("UTR length:", len(transcript.rbs.utr))
    print("CDS length:", len(cds))
    print("Total transcript length:", len(full_seq))

if __name__ == "__main__":
    main()
