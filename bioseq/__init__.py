# bioseq/__init__.py
"""bioseq package root."""

__version__ = "0.1.0"

# re-export common API for convenience
from .analyzer import SequenceRecord, analyze_record  # noqa: F401
