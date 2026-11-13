"""Held-out splits of the validation work (IMPROVEMENTS section H), fixed on 2026-11-11 before any comparison or tuning.

- proteins: sha256(systematic ORF name) mod 10 in 0-4 -> "dev", 5-9 -> "test". Error analysis and rule changes look
  only at dev; test numbers are computed once per finished change.
- diffusion papers: sha256(lower-case DOI) mod 10 < 4 -> "test", else "train". If a compartment would be left without a
  training paper, its test paper with the largest hash stays in train (see diffusion split in diffusion.py).
Do not change these functions: that would leak the test set into tuning.
"""
import hashlib


def _h(s):
    return int(hashlib.sha256(s.encode()).hexdigest(), 16) % 10


def protein_split(orf):
    return "dev" if _h(orf.strip().upper()) < 5 else "test"


def paper_split(doi):
    return "test" if _h(doi.strip().lower()) < 4 else "train"
