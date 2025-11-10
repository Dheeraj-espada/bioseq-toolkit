[![CI](https://github.com/Dheeraj-espada/bioseq-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/Dheeraj-espada/bioseq-toolkit/actions/workflows/ci.yml)

# bioseq-toolkit

**Enhanced Sequence Analyzer** — a Python/BioPython tool for FASTA parsing, ORF finding, translation, and protein analysis.  
Includes unit tests, reproducible conda environment, and CI integration.

---

## 🔍 Features
- Parse FASTA and compute GC%, length, and base composition  
- Find the **longest ORF** (ATG..stop) on both strands  
- Allow **partial ORFs** (ATG..end)  
- Compute protein molecular weight, isoelectric point, and composition  
- Export results to CSV  
- Optionally write ORF nucleotide and amino-acid FASTA files  
- Includes **unit tests** (pytest) and **GitHub Actions CI**

---

## 🚀 Quick Start

```bash
# Clone this repository
git clone https://github.com/Dheeraj-espada/bioseq-toolkit.git
cd bioseq-toolkit

# Create the environment
conda env create -f environment.yml
conda activate bioinfo

# Run analysis
python3 enhanced_seq_analyzer_cli.py input.fasta \
  --orf --both-strands --allow-partial --write-orf-fasta -o results.csv

# Run tests
pytest -q

📂 Repository structure
bioseq-toolkit/
│
├── enhanced_seq_analyzer_cli.py    # Main CLI tool
├── tests/                          # Unit tests
│   └── test_orf_finder.py
├── environment.yml                 # Conda environment
├── .github/workflows/ci.yml        # Continuous Integration workflow
├── .gitignore
└── README.md

📜 License

This project is licensed under the MIT License.

👤 Author

Dheeraj Babu — aspiring Bioinformatician

GitHub: https://github.com/Dheeraj-espada

LinkedIn: https://www.linkedin.com/in/dheeraj-babu-m-74b787219/

## Output directories

- Plots are written under: \$BIOSEQ_OUTPUT_DIR/example_summary_plots (defaults to ./results/example_summary_plots)
- You can override with BIOSEQ_OUTPUT_DIR or pass absolute paths to save_plot functions.
