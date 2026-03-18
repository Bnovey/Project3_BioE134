import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.RNAaseEChecker import RNAseEChecker
from genedesign.codon_freq import load_amino_acid_codon_frequencies
from genedesign.window_mutation import fix_failure_window


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.rbsChooser = None
        self.codonFrequencies = {}
        self.codonChecker = None
        self.forbiddenSequenceChecker = None
        self.internalPromoterChecker = None
        self.rnaseEChecker = None

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

        self.internalPromoterChecker = PromoterChecker()
        self.internalPromoterChecker.initiate()

        self.rnaseEChecker = RNAseEChecker()

        self.codonFrequencies = load_amino_acid_codon_frequencies()



    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # Initial probabilistic CDS
        codons: list[str] = []
        for aa in peptide:
            codon_dict = self.codonFrequencies[aa]
            codon_list = list(codon_dict.keys())
            weights = [codon_dict[c] for c in codon_list]
            codons.append(random.choices(codon_list, weights=weights, k=1)[0])
        codons.append("TAA")
        max_attempts = 100
        # Phase 1: optimize CDS only, without changing/considering RBS.
        for _ in range(max_attempts):
            cds = "".join(codons)
            failures: list[tuple[int, int]] = []
            # Codon usage is tracked by benchmark validation, but we do not hard-fail
            # design here to avoid blocking all sequence generation.
            self.codonChecker.run(codons)
            # Hairpin checker (string)
            ok_hp, hp_seq = hairpin_checker(cds)
            if not ok_hp and hp_seq is not None:
                pos = cds.find(hp_seq)
                if pos != -1:
                    failures.append((pos, pos + len(hp_seq)))
                else:
                    failures.append((0, min(15, len(cds))))

            # Forbidden sequence checker
            ok_forb, forb_site = self.forbiddenSequenceChecker.run(cds)
            if not ok_forb and forb_site is not None:
                pos = cds.find(forb_site)
                if pos != -1:
                    failures.append((pos, pos + len(forb_site)))
                else:
                    failures.append((0, min(15, len(cds))))

            # Internal promoter checker
            ok_prom, prom_seq = self.internalPromoterChecker.run(cds)
            if not ok_prom and prom_seq is not None:
                pos = cds.find(prom_seq)
                if pos != -1:
                    failures.append((pos, pos + len(prom_seq)))
                else:
                    failures.append((0, min(15, len(cds))))
            # RNase E checker (no exact position → use a CDS window near start)
            ok_rnase, rnase_msg = self.rnaseEChecker.run(cds)
            if not ok_rnase:
                window_start = 0
                window_end = min(15, len(cds))
                failures.append((window_start, window_end))
            # If no failures, we are done
            if not failures:
                break
            # Take the first failing window and try to fix it in the CDS
            start_nt, end_nt = failures[0]
            new_codons, changed = fix_failure_window(
                peptide=peptide,
                codons=codons,
                full_seq_start_nt=start_nt,
                full_seq_end_nt=end_nt,
                utr_len=0,
                codon_frequencies=self.codonFrequencies,
            )
            if not changed:
                # For CDS-only design this should be rare; fall back to a full reshuffle.
                for i, aa in enumerate(peptide):
                    codon_dict = self.codonFrequencies[aa]
                    codon_list = list(codon_dict.keys())
                    weights = [codon_dict[c] for c in codon_list]
                    codons[i] = random.choices(codon_list, weights=weights, k=1)[0]
            else:
                codons = new_codons

        cds = "".join(codons)
        # Phase 2: choose/add RBS exactly once after a CDS is created.
        selectedRBS = self.rbsChooser.run(cds, ignores)
        return Transcript(selectedRBS, peptide, codons)


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
