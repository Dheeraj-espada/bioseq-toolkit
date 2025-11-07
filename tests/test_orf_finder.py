"""
tests/test_orf_finder.py
Unit tests for ORF-finding and analysis functions.

Covers:
- SequenceRecord.find_orfs()
- SequenceRecord.longest_orf()
- Wrapper functions find_orfs_in_seq() / find_longest_orf()
- analyze_record() integration
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bioseq import analyzer as ea


def test_find_orfs_simple():
    """Basic ORF: ATG...TAA in the forward strand"""
    seq = Seq("AAATGAAATTTTAA")  # ORF starts at index 2, ends at 13 (TAA)
    orfs = ea.find_orfs_in_seq(seq, allow_partial=False)

    # one clear ORF should be found
    assert len(orfs) == 1
    orf = orfs[0]

    # Check frame and expected protein length
    assert orf["frame"] == 2  # ATG starts at index 2 (frame 2)
    assert orf["aa_len"] == 3  # 9 nt -> 3 amino acids
    assert orf["prot_seq"] == "MKF"  # translated protein sequence


def test_find_orfs_partial():
    """If allow_partial=True, detect ORF with no stop codon at the end"""
    seq = Seq("CCCATGAAAAAA")  # ATG at index 3, no stop codon
    orfs_no = ea.find_orfs_in_seq(seq, allow_partial=False)
    orfs_yes = ea.find_orfs_in_seq(seq, allow_partial=True)

    # Without partial -> no ORF
    assert len(orfs_no) == 0
    # With partial -> one ORF
    assert len(orfs_yes) == 1

    orf = orfs_yes[0]
    assert orf["partial"] is True
    assert orf["aa_len"] >= 1


def test_longest_orf_forward_and_reverse():
    """Find longest ORF on both strands."""
    seq = Seq("AAATGAAATTTTAA")
    longest = ea.find_longest_orf(seq, both_strands=False, allow_partial=False)
    assert longest is not None
    assert longest["strand"] == "+"
    assert longest["aa_len"] == 3

    # Reverse complement should also yield a valid ORF when both_strands=True
    longest_both = ea.find_longest_orf(seq, both_strands=True, allow_partial=False)
    assert longest_both is not None
    assert longest_both["strand"] in {"+", "-"}


def test_sequence_record_methods():
    """Directly test SequenceRecord methods for ORFs."""
    rec = ea.SequenceRecord(id="s1", seq="AAATGAAATTTTAA")
    orfs = rec.find_orfs()
    assert len(orfs) == 1
    longest = rec.longest_orf()
    assert longest["aa_len"] == 3
    assert longest["strand"] == "+"


def test_analyze_record_integration():
    """Integration test: analyze_record should use SequenceRecord under the hood."""
    rec = SeqRecord(Seq("AAATGAAATTTTAA"), id="test1", description="example record")
    res = ea.analyze_record(rec, use_orf=True, both_strands=False, min_orf_aa=0, allow_partial=False)

    # check key fields
    assert res["id"] == "test1"
    assert res["orf_strand"] == "+"
    assert res["prot_len"] == 3
    assert res["orf_start"] != ""
    assert res["gc_percent"] > 0
    assert "rev_comp_first50" in res


