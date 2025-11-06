#!/usr/bin/env python3
"""
bioseq/analyzer.py
Enhanced FASTA analyzer with ORF finder, partial-ORF support, ORF FASTA output,
and optional plotting integration (now uses bioseq.visualizer in-process).

This file was refactored to:
- introduce SequenceRecord dataclass (OOP)
- add find_orfs() and longest_orf() methods on SequenceRecord
- keep backward-compatible wrapper functions for existing CLI/tests
- integrate plotting with bioseq.visualizer via generate_plots_from_csv()
"""
import argparse
import csv
import logging
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# relative import of visualizer module
from . import visualizer as _visualizer
from . import __version__
# constants
STOPS = {"TAA", "TAG", "TGA"}
START = "ATG"


# ---------------------------
# OOP: SequenceRecord class
# ---------------------------
@dataclass
class SequenceRecord:
    """
    Lightweight wrapper around a nucleotide sequence.
    Stores a Bio.Seq.Seq object in self.seq and provides common methods:
    length, gc_content, reverse_complement, translate, find_orfs, longest_orf.
    """
    id: str
    seq: str  # input can be string or Seq-like; will be normalized to Bio.Seq.Seq

    def __post_init__(self):
        # Normalize: remove newlines, uppercase, create Seq object
        seq_str = str(self.seq).replace("\n", "").upper()
        if not seq_str:
            raise ValueError("Empty sequence provided")
        self.seq: Seq = Seq(seq_str)

    def length(self) -> int:
        """Return sequence length (nt)."""
        return len(self.seq)

    def gc_content(self) -> float:
        """Return GC percentage (0-100). Uses simple counting (fast)."""
        n = self.length()
        if n == 0:
            return 0.0
        g = self.seq.count("G")
        c = self.seq.count("C")
        return (g + c) / n * 100

    def reverse_complement(self) -> Seq:
        """Return reverse complement as a Seq object."""
        return self.seq.reverse_complement()

    def translate(self, to_stop: bool = True, table: int = 1) -> str:
        """
        Translate nucleotide seq to amino acids.
        to_stop=True will stop at first stop; table chooses codon table.
        Trim trailing partial codon to avoid warnings.
        """
        trimmed_len = (len(self.seq) // 3) * 3
        trimmed_seq = self.seq[:trimmed_len]
        return str(trimmed_seq.translate(to_stop=to_stop, table=table))

    def to_summary_dict(self) -> Dict[str, object]:
        """Return a small dict summary helpful for DataFrame/CSV output."""
        return {
            "id": self.id,
            "length": self.length(),
            "gc_percent": round(self.gc_content(), 2),
            "rev_comp_first50": str(self.reverse_complement())[:50]
        }

    def find_orfs(self, allow_partial: bool = False) -> List[Dict]:
        """
        Scan this sequence for ORFs that start with ATG and end with a stop codon.
        Returns a list of ORF dicts with keys:
          frame (0-2), start_nt (0-based), end_nt (end exclusive, 0-based), aa_len, prot_seq, partial (bool)
        """
        seq_str = str(self.seq)
        n = len(seq_str)
        orfs: List[Dict] = []
        for frame in range(3):
            i = frame
            while i <= n - 3:
                codon = seq_str[i:i+3]
                if codon == START:
                    j = i + 3
                    found_stop = False
                    while j <= n - 3:
                        stop = seq_str[j:j+3]
                        if stop in STOPS:
                            orf_nt_start = i
                            orf_nt_end = j + 3  # end exclusive index
                            prot = Seq(seq_str[orf_nt_start:orf_nt_end]).translate(to_stop=True)
                            orfs.append({
                                "frame": frame,
                                "start_nt": orf_nt_start,
                                "end_nt": orf_nt_end,
                                "aa_len": len(prot),
                                "prot_seq": str(prot),
                                "partial": False
                            })
                            found_stop = True
                            break
                        j += 3
                    if not found_stop and allow_partial:
                        orf_nt_start = i
                        orf_nt_end = n
                        trimmed_len = ((orf_nt_end - orf_nt_start) // 3) * 3
                        if trimmed_len >= 3:
                            prot = Seq(seq_str[orf_nt_start:orf_nt_start + trimmed_len]).translate(to_stop=True)
                            orfs.append({
                                "frame": frame,
                                "start_nt": orf_nt_start,
                                "end_nt": orf_nt_start + trimmed_len,
                                "aa_len": len(prot),
                                "prot_seq": str(prot),
                                "partial": True
                            })
                    i += 3
                else:
                    i += 3
        return orfs

    def longest_orf(self, both_strands: bool = False, allow_partial: bool = False) -> Optional[Dict]:
        """
        Return the single longest ORF (optionally searching both strands).
        Returns a dict with mapped 1-based coordinates for start/end and strand.
        """
        candidates: List[Dict] = []
        n = len(self.seq)
        # forward strand candidates
        for orf in self.find_orfs(allow_partial=allow_partial):
            candidates.append({
                "strand": "+",
                "frame": orf["frame"],
                "start_nt_1based": orf["start_nt"] + 1,
                "end_nt_1based": orf["end_nt"],
                "aa_len": orf["aa_len"],
                "prot_seq": orf["prot_seq"],
                "partial": orf.get("partial", False)
            })

        # reverse strand candidates (map coordinates back to original)
        if both_strands:
            rc = self.reverse_complement()
            rc_rec = SequenceRecord(id=f"{self.id}_rc", seq=str(rc))
            for orf in rc_rec.find_orfs(allow_partial=allow_partial):
                rc_start = orf["start_nt"]
                rc_end_excl = orf["end_nt"]
                # mapping rc coords back to original (1-based)
                orig_start_1 = n - (rc_end_excl - 1)
                orig_end_1 = n - rc_start
                candidates.append({
                    "strand": "-",
                    "frame": orf["frame"],
                    "start_nt_1based": orig_start_1,
                    "end_nt_1based": orig_end_1,
                    "aa_len": orf["aa_len"],
                    "prot_seq": orf["prot_seq"],
                    "partial": orf.get("partial", False)
                })

        if not candidates:
            return None

        # choose longest by aa length, break ties by nucleotide span
        candidates.sort(key=lambda x: (x["aa_len"], x["end_nt_1based"] - x["start_nt_1based"]), reverse=True)
        best = candidates[0]
        return best


# ---------------------------
# Backwards-compatible wrapper functions
# ---------------------------
def find_orfs_in_seq(seq, allow_partial=False):
    """Wrapper: keep original function signature but use SequenceRecord implementation."""
    return SequenceRecord(id="tmp", seq=str(seq)).find_orfs(allow_partial=allow_partial)


def find_longest_orf(seq, both_strands=False, allow_partial=False):
    """Wrapper: keep original signature but use SequenceRecord.longest_orf."""
    return SequenceRecord(id="tmp", seq=str(seq)).longest_orf(both_strands=both_strands, allow_partial=allow_partial)


# ---------------------------
# analyze_record now uses SequenceRecord internally
# ---------------------------
def analyze_record(rec, use_orf=False, both_strands=False, min_orf_aa=0, allow_partial=False, codon_table=1):
    """
    Accepts a Biopython SeqRecord and returns the summary dict.
    """
    seq_rec = SequenceRecord(id=rec.id, seq=str(rec.seq))

    seq_len = seq_rec.length()
    gc = seq_rec.gc_content()

    orf_info = None
    prot_seq = ""
    prot_len = 0
    prot_mw = ""
    prot_pI = ""
    aa_counts = ""

    if use_orf:
        orf_info = find_longest_orf(seq_rec.seq, both_strands=both_strands, allow_partial=allow_partial)
        if orf_info and orf_info["aa_len"] >= min_orf_aa:
            prot_seq = orf_info["prot_seq"]
            prot_len = orf_info["aa_len"]
            try:
                pa = ProteinAnalysis(prot_seq)
                prot_mw = round(pa.molecular_weight(), 4)
                prot_pI = round(pa.isoelectric_point(), 4)
                aa_counts = ";".join(f"{k}:{v}" for k, v in pa.count_amino_acids().items())
            except Exception as e:
                logging.warning(f"ProteinAnalysis failed for {rec.id}: {e}")
        else:
            orf_info = None
    else:
        # translate trimmed sequence (no ORF detection)
        try:
            prot = Seq(str(seq_rec.seq)[: (seq_len // 3) * 3]).translate(to_stop=True, table=codon_table)
            prot_seq = str(prot)
            prot_len = len(prot_seq)
            if prot_len > 0:
                pa = ProteinAnalysis(prot_seq)
                prot_mw = round(pa.molecular_weight(), 4)
                prot_pI = round(pa.isoelectric_point(), 4)
                aa_counts = ";".join(f"{k}:{v}" for k, v in pa.count_amino_acids().items())
        except Exception as e:
            logging.warning(f"Translation failed for {rec.id}: {e}")

    rev_comp_first50 = str(seq_rec.reverse_complement())[:50]

    if orf_info:
        prot_seq_first50 = prot_seq[:50]
        prot_len_field = prot_len
        prot_mw_field = prot_mw
        prot_pI_field = prot_pI
        aa_counts_field = aa_counts
        orf_frame = orf_info["frame"]
        orf_strand = orf_info["strand"]
        orf_start = orf_info["start_nt_1based"]
        orf_end = orf_info["end_nt_1based"]
        orf_partial = orf_info.get("partial", False)
    else:
        prot_seq_first50 = ""
        prot_len_field = prot_len
        prot_mw_field = prot_mw
        prot_pI_field = prot_pI
        aa_counts_field = aa_counts
        orf_frame = ""
        orf_strand = ""
        orf_start = ""
        orf_end = ""
        orf_partial = ""

    return {
        "id": rec.id,
        "description": rec.description,
        "length": seq_len,
        "gc_percent": round(gc, 2),
        "rev_comp_first50": rev_comp_first50,
        "prot_len": prot_len_field,
        "prot_seq_first50": prot_seq_first50,
        "prot_mw": prot_mw_field,
        "prot_pI": prot_pI_field,
        "aa_counts": aa_counts_field,
        "orf_frame": orf_frame,
        "orf_strand": orf_strand,
        "orf_start": orf_start,
        "orf_end": orf_end,
        "orf_partial": orf_partial
    }


# ---------------------------
# Plot generation helper (integrated visualizer)
# ---------------------------
def generate_plots_from_csv(csv_path, out_dir=None):
    """
    Read CSV produced by analyzer and generate standard plots (length histogram, GC boxplot).
    Returns a dict with paths of generated files.
    """
    csv_path = Path(csv_path)
    if out_dir is None:
        out_dir = csv_path.with_suffix("").stem + "_plots"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = _visualizer.records_from_csv(csv_path)

    length_png = out_dir / "length_histogram.png"
    gc_png = out_dir / "gc_boxplot.png"

    _visualizer.plot_length_histogram(records, length_png)
    _visualizer.plot_gc_boxplot(records, gc_png)

    return {"length_histogram": str(length_png), "gc_boxplot": str(gc_png)}


# ---------------------------
# CLI (kept unchanged except plotting integration)
# ---------------------------
def main():
    parser = argparse.ArgumentParser(description="Enhanced FASTA analyzer with ORF finder and optional plotting")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("in_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("-o", "--out_csv", type=Path, default=Path("enhanced_summary.csv"), help="Output CSV file")
    parser.add_argument("--orf", action="store_true", help="Find longest ORF (ATG..stop) across frames")
    parser.add_argument("--both-strands", action="store_true", help="When using --orf, search both forward and reverse strands")
    parser.add_argument("--min-orf-aa", type=int, default=0, help="Minimum ORF length in amino acids (default 0)")
    parser.add_argument("--allow-partial", action="store_true", help="Allow ATG..end as partial ORF if no in-frame stop")
    parser.add_argument("--write-orf-fasta", action="store_true", help="Write nucleotide and amino-acid FASTA for found ORFs (prefix from out_csv)")
    parser.add_argument("--plot-stats", action="store_true", help="Generate GC% and protein-length plots from CSV after analysis")
    parser.add_argument("--to-stop", action="store_true", help="(Legacy) translate to first stop (used when --orf not set)")
    parser.add_argument("--codon-table", type=int, default=1, help="Codon table number (default 1)")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum nucleotide length to include (default 0)")
    parser.add_argument("--min-prot-len", type=int, default=0, help="Minimum protein length (aa) to include (default 0, used when --orf not set)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else (
        logging.WARNING if args.quiet else logging.INFO
    )
    logging.basicConfig(
        level=log_level,
        format="%(levelname)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

# Optional logfile handler
    if args.logfile:
        fh = logging.FileHandler(args.logfile)
        fh.setLevel(logging.DEBUG)
    	fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    	logging.getLogger().addHandler(fh)

# Keep this check right after logging is configured
    if not args.in_fasta.exists():
    	logging.error(f"Input FASTA {args.in_fasta} not found.")
    	raise SystemExit(1)

    # If writing ORF FASTA, load entire FASTA into memory (id -> SeqRecord)
    seq_map = None
    orf_nuc_fh = None
    orf_aa_fh = None
    if args.write_orf_fasta:
        logging.info("Loading FASTA into memory for ORF FASTA output (mem mode).")
        seq_map = {rec.id: rec for rec in SeqIO.parse(str(args.in_fasta), "fasta")}
        prefix = args.out_csv.stem
        orf_nuc_path = Path(f"{prefix}_orfs.nuc.fa")
        orf_aa_path = Path(f"{prefix}_orfs.aa.fa")
        orf_nuc_fh = orf_nuc_path.open("w")
        orf_aa_fh = orf_aa_path.open("w")
        logging.info(f"ORF nucleotide FASTA will be written to: {orf_nuc_path}")
        logging.info(f"ORF amino-acid FASTA will be written to: {orf_aa_path}")

    # Main streaming parse for analysis & CSV writing
    records = SeqIO.parse(str(args.in_fasta), "fasta")
    written = 0
    skipped_length = 0
    skipped_prot = 0
    total = 0

    with args.out_csv.open("w", newline="") as fh:
        fieldnames = ["id","description","length","gc_percent","rev_comp_first50",
                      "prot_len","prot_seq_first50","prot_mw","prot_pI","aa_counts",
                      "orf_frame","orf_strand","orf_start","orf_end","orf_partial"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for rec in records:
            total += 1
            res = analyze_record(rec, use_orf=args.orf, both_strands=args.both_strands,
                                 min_orf_aa=args.min_orf_aa, allow_partial=args.allow_partial,
                                 codon_table=args.codon_table)

            # nucleotide-length filter
            if res["length"] < args.min_length:
                skipped_length += 1
                logging.debug(f"Skipping {res['id']} (len {res['length']} < min-length {args.min_length})")
                continue

            # protein-length filter
            if args.orf:
                prot_len_field = res["prot_len"]
                if prot_len_field == 0 or prot_len_field < args.min_orf_aa:
                    skipped_prot += 1
                    logging.debug(f"Skipping {res['id']} (orf prot_len {prot_len_field} < min-orf-aa {args.min_orf_aa})")
                    continue
            else:
                if res["prot_len"] < args.min_prot_len:
                    skipped_prot += 1
                    logging.debug(f"Skipping {res['id']} (prot_len {res['prot_len']} < min-prot-len {args.min_prot_len})")
                    continue

            # write CSV row
            writer.writerow(res)
            written += 1

            # write ORF FASTAs if requested
            if args.write_orf_fasta and res.get("orf_start") and res.get("orf_end"):
                sid = res["id"]
                seq_obj = seq_map.get(sid)
                if seq_obj is None:
                    logging.warning(f"Sequence {sid} not found in in-memory FASTA (unexpected).")
                else:
                    start = int(res["orf_start"])
                    end = int(res["orf_end"])
                    strand = res["orf_strand"]
                    frame = res["orf_frame"]
                    partial = res["orf_partial"]
                    aa_len = res["prot_len"]
                    nuc = seq_obj.seq[start-1:end]
                    if strand == "-":
                        nuc = nuc.reverse_complement()
                    aa = res["prot_seq_first50"] if res["prot_seq_first50"] else str(nuc.translate(to_stop=True))
                    header_base = f"{sid}|{start}-{end}|strand={strand}|frame={frame}|aa_len={aa_len}"
                    if partial:
                        header_base += "|partial"
                    orf_nuc_fh.write(f">{header_base}\n{str(nuc)}\n")
                    orf_aa_fh.write(f">{header_base}\n{aa}\n")

    logging.info(f"Processed {total} records. Written: {written}. Skipped (length): {skipped_length}. Skipped (prot_len): {skipped_prot}. Output CSV: {args.out_csv}")

    # integrated plotting using visualizer
    if args.plot_stats:
        try:
            outdir = args.out_csv.with_suffix("").stem + "_plots"
            logging.info("Generating plots (in-process) using bioseq.visualizer")
            generated = generate_plots_from_csv(args.out_csv, out_dir=outdir)
            logging.info(f"Plots generated: {generated}")
        except Exception as e:
            logging.warning(f"Failed to generate plots: {e}")

    if orf_nuc_fh:
        orf_nuc_fh.close()
    if orf_aa_fh:
        orf_aa_fh.close()


if __name__ == "__main__":
    main()
