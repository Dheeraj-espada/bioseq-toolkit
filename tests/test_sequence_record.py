# tests/test_sequence_record.py
import pytest
from bioseq.analyzer import SequenceRecord
from Bio.Seq import Seq

def test_gc_content_basic():
    # For "ATGCGC": G=2, C=2 => GC = 4/6 * 100 = 66.666...
    rec = SequenceRecord(id="r1", seq="ATGCGC")
    assert rec.length() == 6
    assert pytest.approx(rec.gc_content(), rel=1e-6) == (4/6 * 100)

def test_gc_content_lowercase_and_newlines():
    # input normalization: lowercase + newlines should be handled
    raw = "atgc\ngc"
    rec = SequenceRecord(id="r2", seq=raw)
    assert rec.length() == 6
    assert pytest.approx(rec.gc_content(), rel=1e-6) == (4/6 * 100)

def test_reverse_complement_simple():
    rec = SequenceRecord(id="r3", seq="ATGC")
    # reverse complement of ATGC -> GCAT
    rc = rec.reverse_complement()
    # rc is a Bio.Seq.Seq object; cast to str for comparison
    assert str(rc) == "GCAT"

def test_reverse_complement_palindrome():
    # palindromic sequence should equal its reverse complement
    seq = "GAATTC"  # EcoRI site; RC == original
    rec = SequenceRecord(id="r4", seq=seq)
    assert str(rec.reverse_complement()) == seq
