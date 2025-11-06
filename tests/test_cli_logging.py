# tests/test_cli_logging.py
import subprocess, sys
from pathlib import Path

def test_cli_quiet_and_logfile(tmp_path):
    log_file = tmp_path / "log.txt"
    res = subprocess.run(
        [sys.executable, "-m", "bioseq.analyzer", "--version", "--quiet", "--logfile", str(log_file)],
        capture_output=True, text=True
    )
    assert res.returncode == 0
    # logfile should be created
    assert log_file.exists()
    # should contain at least one line
    assert log_file.read_text().strip() != ""
