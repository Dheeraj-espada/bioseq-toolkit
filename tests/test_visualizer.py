# tests/test_visualizer.py
import os
from pathlib import Path
from bioseq import visualizer as vz

def sample_records():
    return [
        {"id": "r1", "length": 1000, "gc_percent": 50.0, "prot_len": 100},
        {"id": "r2", "length": 500, "gc_percent": 40.0, "prot_len": 0},
        {"id": "r3", "length": 1500, "gc_percent": 60.0, "prot_len": 200},
        {"id": "r4", "length": 200, "gc_percent": 30.0, "prot_len": 0},
    ]

def test_plot_length_histogram(tmp_path):
    out = tmp_path / "length_hist.png"
    records = sample_records()
    path = vz.plot_length_histogram(records, out)
    assert Path(path).exists()
    assert Path(path).stat().st_size > 0

def test_plot_gc_boxplot(tmp_path):
    out = tmp_path / "gc_box.png"
    records = sample_records()
    path = vz.plot_gc_boxplot(records, out)
    assert Path(path).exists()
    assert Path(path).stat().st_size > 0

def test_records_from_csv(tmp_path):
    # create a tiny CSV similar to analyzer output
    csv_path = tmp_path / "test.csv"
    csv_path.write_text("id,description,length,gc_percent,rev_comp_first50,prot_len,prot_seq_first50,prot_mw,prot_pI,aa_counts,orf_frame,orf_strand,orf_start,orf_end,orf_partial\n"
                        "r1,,1000,50.0,,,,,,,,,,,\n"
                        "r2,,500,40.0,,,,,,,,,,,\n")
    recs = vz.records_from_csv(csv_path)
    assert isinstance(recs, list)
    assert len(recs) == 2
    assert recs[0]["length"] == 1000
    assert recs[1]["gc_percent"] == 40.0

