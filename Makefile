# Makefile — simple developer conveniences
PY := python
PIP := pip

.PHONY: help install test run clean

help:
	@echo "Targets:"
	@echo "  make install   - pip install -e ."
	@echo "  make test      - run pytest"
	@echo "  make run       - run batch on sample FASTA(s)"
	@echo "  make clean     - remove results and caches"

install:
	$(PIP) install -e .

test:
	pytest -q

run:
	bash scripts/run_analysis.sh test_data/*.fasta

clean:
	rm -rf results .pytest_cache
