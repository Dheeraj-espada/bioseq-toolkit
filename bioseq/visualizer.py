# bioseq/visualizer.py
"""
Simple visualization utilities for BioSeq Toolkit.

Functions accept a list of record-summary dicts (as produced by analyze_record -> CSV).
Each plotting function writes a PNG file and returns the path.

Design choices:
- Keep plotting functions pure-ish: they accept data and an output path.
- Use matplotlib (no seaborn) and avoid specifying colors so visuals are flexible.
- Tests will check that files are created and non-empty.
"""
from typing import List, Dict, Union
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Optional

def _ensure_plot_dir(base_output: Optional[str] = None, subdir: str = "example_summary_plots") -> str:
    """
    Ensure plots directory exists and return full path.
    Priority:
      1. base_output argument
      2. BIOSEQ_OUTPUT_DIR environment variable
      3. default "results"
    Final path returned: <base>/<subdir>
    """
    base = base_output or __import__("os").environ.get("BIOSEQ_OUTPUT_DIR") or "results"
    plot_dir = __import__("os").path.join(base, subdir)
    __import__("os").makedirs(plot_dir, exist_ok=True)
    return plot_dir



def save_plot(fig, filename: str, base_output: Optional[str] = None):
    """
    Save a matplotlib figure.

    Rules:
      - Absolute path -> save exactly there (create parents).
      - Relative path:
          * if it begins with the plots-subdir (example_summary_plots), strip that
            and place the remainder under the controlled plots dir returned by
            _ensure_plot_dir(base_output) so we do NOT duplicate the subdir.
          * otherwise, place the relative path under the controlled plots dir,
            preserving subdirectories.
    Returns the full path to the saved file.
    """
    import os
    subdir = "example_summary_plots"

    # absolute path: respect it exactly
    if os.path.isabs(filename):
        full = os.path.abspath(filename)
        parent = os.path.dirname(full)
        if parent:
            os.makedirs(parent, exist_ok=True)
        fig.savefig(full, bbox_inches="tight")
        return full

    # relative path: normalize and sanitize
    rel = os.path.normpath(filename)

    # avoid writing outside of plots dir (no leading ..)
    if rel.startswith(os.pardir + os.sep):
        # fallback: use basename
        rel = os.path.basename(rel)

    # if the relative path already begins with the subdir, strip that component
    # so we don't create results/example_summary_plots/example_summary_plots/...
    comps = rel.split(os.path.sep)
    if comps and comps[0] == subdir:
        comps = comps[1:]
        rel = os.path.join(*comps) if comps else os.path.basename(filename)

    out_dir = _ensure_plot_dir(base_output, subdir=subdir)
    full = os.path.join(out_dir, rel)
    parent = os.path.dirname(full)
    if parent:
        os.makedirs(parent, exist_ok=True)
    fig.savefig(full, bbox_inches="tight")
    return full

def _ensure_dir(path: Union[str, Path]):
    p = Path(path)
    if not p.parent.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
    return p


def records_from_csv(csv_path: Union[str, Path]) -> List[Dict]:
    """Read CSV from analyzer and return list of dicts with parsed numeric fields."""
    records = []
    with open(str(csv_path), newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            # parse numeric fields we care about; keep others as-is
            parsed = dict(row)
            # safe numeric parsing
            try:
                parsed["length"] = int(row.get("length", 0))
            except Exception:
                parsed["length"] = 0
            try:
                parsed["gc_percent"] = float(row.get("gc_percent", 0.0))
            except Exception:
                parsed["gc_percent"] = 0.0
            try:
                parsed["prot_len"] = int(row.get("prot_len", 0) or 0)
            except Exception:
                parsed["prot_len"] = 0
            records.append(parsed)
    return records


def plot_length_histogram(records: List[Dict], out_path: Union[str, Path], bins: int = 30) -> Path:
    """
    Plot histogram of sequence lengths and save to out_path (PNG).
    records: list of dicts containing "length" key.
    """
    out_path = _ensure_dir(out_path)
    lengths = [r.get("length", 0) for r in records if r.get("length", 0) is not None]
    if not lengths:
        # create an empty figure to keep behavior deterministic for tests
        fig = plt.figure()
        fig.suptitle("No length data")
        out_path = save_plot(fig, str(out_path))
        plt.close(fig)
        return out_path

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(lengths, bins=bins)
    ax.set_xlabel("Sequence length (nt)")
    ax.set_ylabel("Count")
    ax.set_title("Sequence Length Distribution")
    fig.tight_layout()
    out_path = save_plot(fig, str(out_path))
    plt.close(fig)
    return out_path


def plot_gc_boxplot(records: List[Dict], out_path: Union[str, Path]) -> Path:
    """
    Plot GC% boxplot grouped by presence/absence of protein (prot_len > 0).
    If prot_len is missing, all values are put into a single box.
    """
    out_path = _ensure_dir(out_path)
    gc_values_with_prot = []
    gc_values_without_prot = []
    has_prot = False
    for r in records:
        gc = r.get("gc_percent", None)
        if gc is None:
            continue
        prot_len = r.get("prot_len", None)
        if prot_len is None:
            # fallback to single group
            gc_values_with_prot.append(gc)
        else:
            has_prot = True
            if prot_len > 0:
                gc_values_with_prot.append(gc)
            else:
                gc_values_without_prot.append(gc)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    if has_prot:
        data = []
        labels = []
        if gc_values_with_prot:
            data.append(gc_values_with_prot)
            labels.append("with_protein")
        if gc_values_without_prot:
            data.append(gc_values_without_prot)
            labels.append("without_protein")
        if not data:
            ax.text(0.5, 0.5, "No GC data", ha="center")
        else:
            ax.boxplot(data)
            ax.set_xticklabels(labels)
    else:
        # single boxplot
        if not gc_values_with_prot:
            ax.text(0.5, 0.5, "No GC data", ha="center")
        else:
            ax.boxplot([gc_values_with_prot])
            ax.set_xticklabels(["gc_percent"])

    ax.set_ylabel("GC%")
    ax.set_title("GC% distribution")
    fig.tight_layout()
    out_path = save_plot(fig, str(out_path))
    plt.close(fig)
    return out_path
