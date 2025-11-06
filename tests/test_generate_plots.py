# tests/test_generate_plots.py
from pathlib import Path
from bioseq import analyzer as ea

def test_generate_plots_from_csv(tmp_path):
    # create a tiny CSV similar to analyzer output
    csv_path = tmp_path / "small.csv"
    csv_path.write_text(
        "id,description,length,gc_percent,rev_comp_first50,prot_len,prot_seq_first50,prot_mw,prot_pI,aa_counts,orf_frame,orf_strand,orf_start,orf_end,orf_partial\n"
        "r1,,1000,50.0,,,,,,,,,,,\n"
        "r2,,500,40.0,,,,,,,,,,,\n"
    )
    outdir = tmp_path / "plots_out"
    res = ea.generate_plots_from_csv(csv_path, out_dir=outdir)

    # check that returned paths exist and files are non-empty
    length_png = Path(res["length_histogram"])
    gc_png = Path(res["gc_boxplot"])
    assert length_png.exists()
    assert length_png.stat().st_size > 0
    assert gc_png.exists()
    assert gc_png.stat().st_size > 0

