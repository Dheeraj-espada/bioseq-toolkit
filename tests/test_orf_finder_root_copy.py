# tests/test_orf_finder.py
import io
from Bio.Seq import Seq
import enhanced_seq_analyzer_cli as ea

def test_find_orfs_simple():
    seq = Seq("AAATGAAATTTTAA")  # ATG...TAA -> one ORF
    orfs = ea.find_orfs_in_seq(seq, allow_partial=False)
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf["frame"] == 2  # ATG starts at index 2 (0-based)
    assert orf["aa_len"] == 3  # translated aa length
    assert orf["prot_seq"].startswith("MKF") or orf["prot_seq"]  # common check

def test_find_orfs_allow_partial():
    seq = Seq("CCCATGAAAAAA")  # ATG but no stop -> partial allowed
    orfs_no = ea.find_orfs_in_seq(seq, allow_partial=False)
    assert len(orfs_no) == 0
    orfs_yes = ea.find_orfs_in_seq(seq, allow_partial=True)
    assert len(orfs_yes) >= 1
    partial = orfs_yes[0]
    assert partial["partial"] is True
    assert partial["aa_len"] >= 1

def test_find_longest_orf_both_strands():
    # forward ORF
    seq_fwd = Seq("AAATGAAATTTTAA")
    best = ea.find_longest_orf(seq_fwd, both_strands=False, allow_partial=False)
    assert best is not None
    assert best["strand"] == "+"
    assert best["aa_len"] == 3

    # reverse ORF: design sequence whose reverse complement has an ATG..stop
    forward_orf = Seq("AAATGAAATTTTAA")            # contains ATG...TAA
    seq_rev = forward_orf.reverse_complement()     # reverse complement will contain that ORF
    best_rev = ea.find_longest_orf(seq_rev, both_strands=True, allow_partial=False)
    assert best_rev is not None
    assert best_rev["strand"] in ("+","-")

def test_analyze_record_translation_and_fields():
    # Build a mock record using Bio.SeqIO style object -> use SeqIO to create a record
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("AAATGAAATTTTAA"), id="test1", description="test1")
    res = ea.analyze_record(rec, use_orf=True, both_strands=False, min_orf_aa=0, allow_partial=False)
    # should find ORF
    assert res["orf_start"] != ""
    assert res["prot_len"] == 3
    assert res["orf_strand"] == "+"
