# tests/test_cli_entrypoint.py
import subprocess
import sys

def test_cli_version():
    # call module entrypoint (works without console script)
    res = subprocess.run([sys.executable, "-m", "bioseq.analyzer", "--version"], capture_output=True, text=True)
    assert res.returncode == 0
    assert res.stdout.strip() != ""  # prints version to stdout
