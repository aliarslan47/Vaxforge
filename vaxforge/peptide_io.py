"""Aday peptitlerin tam dökümü (`peptides.json`) — rapor YENİDEN üretilebilsin.

Neden var: rapor (HTML/PDF) `Peptide` NESNELERİNDEN üretilir, `run.json`'dan
değil. `run.json` yalnız CSV sütunlarını (`candidates`) taşıdığı için bir
rapordaki düzen/etiket hatasını düzeltmek koşuyu BAŞTAN ÇALIŞTIRMAK gerektiriyordu
(BoHV-1 için ~1 saat + NetMHCpan). Bu dosya, koşunun peptit metriklerini
kayıpsız yazar; `scripts/rebuild_report.py` bununla raporu yeniden üretir.

Biçim: {"version": 1, "n": N, "peptides": [{seq, parent, kind, start, metrics,
methods, candidacy, passed, notes}, ...]} — `Peptide` alanlarının birebir
karşılığı; JSON'a girmeyen değerler `default=str` ile yazılır.
"""

from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path

from .models import Peptide

VERSION = 1
FILENAME = "peptides.json"


def to_dict(p: Peptide) -> dict:
    return {f.name: getattr(p, f.name) for f in fields(Peptide)}


def dump(peptides, path: Path) -> Path:
    path = Path(path)
    path.write_text(json.dumps(
        {"version": VERSION, "n": len(peptides),
         "peptides": [to_dict(p) for p in peptides]},
        ensure_ascii=False, default=str), encoding="utf-8")
    return path


def load(path: Path) -> list[Peptide]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    names = {f.name for f in fields(Peptide)}
    out = []
    for d in data.get("peptides", []):
        out.append(Peptide(**{k: v for k, v in d.items() if k in names}))
    return out
