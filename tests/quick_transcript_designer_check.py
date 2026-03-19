from tests.benchmarking.proteome_benchmarker import parse_fasta
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.RNAaseEChecker import RNAseEChecker
from genedesign.seq_utils.Translate import Translate


def main() -> None:
    proteome = parse_fasta('tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta')
    sample = list(proteome.items())[:10]

    designer = TranscriptDesigner(); designer.initiate()
    forbidden_checker = ForbiddenSequenceChecker(); forbidden_checker.initiate()
    promoter_checker = PromoterChecker(); promoter_checker.initiate()
    codon_checker = CodonChecker(); codon_checker.initiate()
    rnase_checker = RNAseEChecker()
    translator = Translate(); translator.initiate()

    summary = {k: 0 for k in ['translation', 'forbidden', 'promoter', 'hairpin', 'rnaseE', 'codon', 'overall']}

    for gene, protein in sample:
        transcript = designer.run(protein, set())
        cds = ''.join(transcript.codons)
        seq = transcript.rbs.utr.upper() + cds

        translation_ok = (
            len(cds) % 3 == 0
            and translator.run(cds) == protein
            and cds.startswith(("ATG", "GTG", "TTG"))
            and cds.endswith(("TAA", "TAG", "TGA"))
        )
        forbidden_ok, _ = forbidden_checker.run(seq)
        promoter_ok, _ = promoter_checker.run(seq)
        hairpin_ok, _ = hairpin_checker(seq)
        rnase_ok, _ = rnase_checker.run(seq)
        codon_ok, diversity, rare_count, cai = codon_checker.run(transcript.codons)

        overall_ok = all([translation_ok, forbidden_ok, promoter_ok, hairpin_ok, rnase_ok, codon_ok])

        summary['translation'] += int(translation_ok)
        summary['forbidden'] += int(forbidden_ok)
        summary['promoter'] += int(promoter_ok)
        summary['hairpin'] += int(hairpin_ok)
        summary['rnaseE'] += int(rnase_ok)
        summary['codon'] += int(codon_ok)
        summary['overall'] += int(overall_ok)

        print(
            f"{gene}\toverall={overall_ok}\ttranslation={translation_ok}\tforbidden={forbidden_ok}"
            f"\tpromoter={promoter_ok}\thairpin={hairpin_ok}\trnaseE={rnase_ok}\tcodon={codon_ok}"
            f"\tdiv={diversity:.3f}\trare={rare_count}\tcai={cai:.3f}"
        )

    print("\nSUMMARY")
    print(f"tested={len(sample)}")
    for key in ['translation', 'forbidden', 'promoter', 'hairpin', 'rnaseE', 'codon', 'overall']:
        print(f"{key}_pass={summary[key]}/{len(sample)}")


if __name__ == '__main__':
    main()
