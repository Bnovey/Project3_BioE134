import random
import csv
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
    Reverse translates a protein sequence into a DNA sequence using a
    sliding-window algorithm with guided random codon sampling.

    For each 3-aa window, sample a limited number of RSCU-weighted codon
    combinations and keep the best one in context.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.codonChecker = None
        self.forbiddenSequenceChecker = None
        self.promoterChecker = None
        self.rnaseEChecker = None
        self.aaToSynonymousCodons = {}
        self.aaToWeights = {}

        self.commitCodonsPerStep = 3
        self.downstreamLookaheadCodons = 6
        self.preambleCodons = 8
        # Needs to be 51 for hairpin checker
        # Main speed/quality knob
        self.samplesPerWindow = 20

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.forbiddenSequenceChecker = ForbiddenSequenceChecker()
        self.forbiddenSequenceChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.rnaseEChecker = RNAseEChecker()

        # Highest-CAI seed codon for each amino acid
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

        codon_usage_path = Path(__file__).resolve().parent / "data" / "codon_usage.txt"
        aa_codon_freq = {}

        with open(codon_usage_path, "r") as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 3:
                    continue

                codon = row[0].strip().upper()
                aa = row[1].strip()

                try:
                    freq = float(row[2].strip())
                except ValueError:
                    continue

                if aa == '*':
                    continue

                aa_codon_freq.setdefault(aa, []).append((codon, freq))

        for aa, codon_freqs in aa_codon_freq.items():
            self.aaToSynonymousCodons[aa] = [cf[0] for cf in codon_freqs]
            self.aaToWeights[aa] = [cf[1] for cf in codon_freqs]

    def _initialize_codons_highest_cai(self, peptide: str) -> list[str]:
        return [self.aminoAcidToCodon[aa] for aa in peptide]

    def _sample_commit(self, window_aas: str) -> list[str]:
        return [
            random.choices(
                self.aaToSynonymousCodons[aa],
                weights=self.aaToWeights[aa],
                k=1
            )[0].upper()
            for aa in window_aas
        ]

    def _score(self, coding_codons: list[str], dna_context_codons: list[str]) -> int:
        coding_codons = [c.upper() for c in coding_codons]
        dna = "".join(c.upper() for c in dna_context_codons)

        score = 0

        if not self.codonChecker.run(coding_codons)[0]:
            score += 1
        if not self.forbiddenSequenceChecker.run(dna)[0]:
            score += 1
        if not self.promoterChecker.run(dna)[0]:
            score += 1
        if not self.rnaseEChecker.run(dna)[0]:
            score += 1
        if not hairpin_checker(dna)[0]:
            score += 1

        return score

    def _optimize(self, peptide: str, seed_codons: list[str], utr: str) -> list[str]:
        n = len(peptide)
        step = self.commitCodonsPerStep
        final_codons = [c.upper() for c in seed_codons]
        utr = utr.upper()

        for i in range(0, n, step):
            commit_end = min(n, i + step)
            lookahead_end = min(n, commit_end + self.downstreamLookaheadCodons)

            preamble_start = max(0, i - self.preambleCodons)
            preamble = final_codons[preamble_start:i]

            if i < self.preambleCodons:
                missing = self.preambleCodons - len(preamble)
                utr_tail = utr[-3 * missing:] if missing > 0 else ""
                utr_chunks = [utr_tail[j:j + 3] for j in range(0, len(utr_tail), 3)]
                preamble = utr_chunks + preamble

            current_commit = final_codons[i:commit_end]
            lookahead = final_codons[commit_end:lookahead_end]
            window_aas = peptide[i:commit_end]

            best_commit = current_commit[:]
            best_score = self._score(
                coding_codons=current_commit + lookahead,
                dna_context_codons=preamble + current_commit + lookahead
            )

            for _ in range(self.samplesPerWindow):
                trial_commit = self._sample_commit(window_aas)

                trial_score = self._score(
                    coding_codons=trial_commit + lookahead,
                    dna_context_codons=preamble + trial_commit + lookahead
                )

                if trial_score < best_score:
                    best_score = trial_score
                    best_commit = trial_commit

                    if best_score == 0:
                        break

            final_codons[i:commit_end] = best_commit

        return final_codons

    def run(self, peptide: str, ignores: set) -> Transcript:
        # deterministic high-CAI seed is fine and fast
        seed_codons = self._initialize_codons_highest_cai(peptide)

        selected_rbs = self.rbsChooser.run(seed_codons, ignores)
        utr = selected_rbs.utr.upper()

        cds = self._optimize(peptide, seed_codons, utr)
        cds.append("TAA")

        return Transcript(selected_rbs, peptide, cds)


def main() -> None:
    peptide = (
        "MNPSDVFQIIEGHTKLMRDSIPLIASENLTSLSVRRCYVSDLGHRYAEGRVGERFYEGCKYVDQIESMAIELTRK"
        "IFEAEHANVQPISGVVANLAAFFALTNVGDTIMSISVPCGGHISHDRVSAAGLRGLRVIHYPFNSEEMSVDVDETR"
    )

    designer = TranscriptDesigner()
    designer.initiate()
    transcript = designer.run(peptide, set())

    cds = "".join(transcript.codons)
    full_seq = transcript.rbs.utr.upper() + cds

    print("Selected RBS:", transcript.rbs.gene_name)
    print("CDS length:", len(cds))
    print("Total transcript length:", len(full_seq))


if __name__ == "__main__":
    main()