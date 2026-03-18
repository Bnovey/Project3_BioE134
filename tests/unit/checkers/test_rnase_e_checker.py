import pytest
from genedesign.checkers.RNAaseEChecker import RNAseEChecker


@pytest.fixture
def checker():
    return RNAseEChecker()


VALID_SEQUENCES = [
    "ATGCGCATGCGCATGCGC",   # GC-rich repeat
    "ATGCGTACGTAGCTAGCT",   # balanced coding
    "GCTAGCTAGCTAGCTAGC",   # no AT-rich stretches
    "CGCGCTATGCGCTATGCG",   # mixed GC/AT, no runs
    "ATGCTAGCATGCTAGCAT",   # alternating GC/AT balanced
    "CGTACGTAGCTAGCTAGC",   # spread AT content
    "TAGCTAGCTAGCTAGCTA",   # uniform distribution
    "ATGCGCTAGCATGCGTAG",   # random-looking but safe
    "GCATGCATGCATGCATGC",   # GC-AT alternation
    "ATGCGATCGATCGATCGA",   # regular ATGC pattern
    "CGATCGATCGATCGATCG",   # CG repeat
    "TACGTACGTACGTACGTA",   # TACG repeat
    "GCTAGCTAGCATCGATCG",   # mixed safe sequence
    "ATGCGTTAGCATCGTAGC",   # coding-like
    "CGTAGCTAGCATGCATGC",   # safe AT spread
    "ACGTACGTACGTACGTAC",   # ACGT repeat
    "TAGCATGCTAGCATGCTA",   # interleaved
    "GTACGTACGTACGTACGT",   # GTAC repeat
    "CATGCATGCATGCATGCA",   # CATG repeat
    "ATGCATGCATGCATGCAT",   # ATGC repeat
    "CGCGATCGCGATCGCGAT",   # CG-heavy safe
    "GCATCGATCGATGCATCG",   # balanced mixed
    "TAGCGATCGATCGTAGCG",   # spread runs under 6
    "ATGCGTAGCATCGTAGCT",   # typical gene fragment
    "GCTAGCATGCGCTAGCAT",   # no high AU windows
]


def test_valid_sequences(checker):
    for seq in VALID_SEQUENCES:
        valid, msg = checker.run(seq)
        assert valid is True, f"Expected valid but got failure for '{seq}': {msg}"


def test_poly_t_run_fails(checker):
    """A run of 6+ T's should fail the poly-A/T check."""
    valid, msg = checker.run("ATGCTTTTTTATGCATGC")
    assert valid is False
    assert "Poly-T" in msg


def test_poly_a_run_fails(checker):
    """A run of 6+ A's should fail the poly-A/T check."""
    valid, msg = checker.run("ATGCAAAAAATGCATGCA")
    assert valid is False
    assert "Poly-A" in msg

def test_high_au_window_fails(checker):
    """A window with >75% AT content should fail."""
    valid, msg = checker.run("ATATATATATCGTACGTA")
    assert valid is False
    assert "AU content" in msg

def test_exactly_70_percent_au_passes(checker):
    """Exactly 70% AT (not exceeding threshold) should pass."""
    valid, msg = checker.run("GCATTAATGC")
    assert valid is True


def test_strict_au_threshold_50_percent_fails():
    """Strict 50% AU threshold should fail sequences that pass the default 75%."""
    strict = RNAseEChecker(au_threshold=0.50)
    # "ATGCATGCAT" = 6 AT / 10 = 60% AU — passes default but fails strict 50%
    valid, msg = strict.run("ATGCATGCATGCATGCAT")
    assert valid is False
    assert "AU content" in msg


def test_large_window_size_misses_local_peak():
    """A larger window (20) averages out a local AU-rich region that a smaller window catches."""
    seq = "GCGCGCATATATATATGCGC"
    small_window = RNAseEChecker(window_size=6, au_threshold=0.75)
    large_window = RNAseEChecker(window_size=20, au_threshold=0.75)
    valid_small, _ = small_window.run(seq)
    valid_large, _ = large_window.run(seq)
    assert valid_small is False   # small window catches the hot spot
    assert valid_large is True    # large window averages it out


def test_custom_run_length_stricter():
    """Custom run_length=4 should catch a run of 4 A's."""
    strict_checker = RNAseEChecker(run_length=4)
    valid, msg = strict_checker.run("GCAAAACGCATGCATG")
    assert valid is False
    assert "Poly-A" in msg


def test_custom_run_length_relaxed():
    """Custom run_length=8 should pass a run of 6 T's that would normally fail."""
    relaxed = RNAseEChecker(run_length=8)
    valid, msg = relaxed.run("GCTTTTTTCGCGCGCGCCGCGC")
    print(msg)
    assert valid is True

