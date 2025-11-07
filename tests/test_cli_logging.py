# tests/test_cli_logging.py
import subprocess
import sys
from pathlib import Path

def test_cli_quiet_and_logfile(tmp_path):
    # create a tiny FASTA file so the CLI proceeds (argparse won't exit early)
    fasta = tmp_path / "small.fa"
    fasta.write_text(">r1\nATGAAATAA\n")

    log_file = tmp_path / "log.txt"
    out_csv = tmp_path / "out.csv"

    # run the analyzer providing an actual input file so main() continues and configures logging
    res = subprocess.run(
        [
            sys.executable,
            "-m",
            "bioseq.analyzer",
            str(fasta),
            "--quiet",
            "--logfile",
            str(log_file),
            "--out_csv",
            str(out_csv),
        ],
        capture_output=True,
        text=True,
    )
    # should exit normally
    assert res.returncode == 0

    # logfile should be created and non-empty
    assert log_file.exists()
    assert log_file.stat().st_size > 0
