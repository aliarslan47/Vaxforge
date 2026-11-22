#!/usr/bin/env python3
"""Var olan bir koşu klasöründen raporu YENİDEN üretir — koşuyu TEKRARLAMADAN.

    .venv/bin/python scripts/rebuild_report.py outputs/run_2026-11-22T20-34-48

Okur:  run.json (meta) + peptides.json (aday peptitler, bkz. vaxforge/peptide_io)
Yazar: report.html, report.pdf  (aynı klasöre; diğer dosyalara DOKUNMAZ)

Neden: rapordaki bir düzen/etiket hatasını düzeltmek için hattı baştan
çalıştırmak gerekmiyor. Hat hiçbir biçimde ÇAĞRILMAZ — ne NetMHCpan, ne ağ,
ne de yeni bir koşu. `peptides.json` yoksa (o dosya eklenmeden önce üretilmiş
eski koşular) açıkça söyler ve çıkar.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vaxforge import peptide_io, report  # noqa: E402


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    outdir = Path(argv[1])
    meta_p, pep_p = outdir / "run.json", outdir / peptide_io.FILENAME
    for p in (meta_p, pep_p):
        if not p.exists():
            print(f"HATA: {p} yok.", file=sys.stderr)
            if p is pep_p:
                print("Bu koşu peptides.json'dan ÖNCE üretilmiş; raporu yeniden "
                      "üretmek için peptit metrikleri gerekiyor.", file=sys.stderr)
            return 1

    meta = json.loads(meta_p.read_text(encoding="utf-8"))
    meta.pop("candidates", None)          # rapor df'yi peptitlerden kurar
    peptides = peptide_io.load(pep_p)
    print(f"{len(peptides)} peptit, meta {meta.get('timestamp')} okundu.")

    df = report._candidates_df(peptides, meta.get("lang", "tr"))
    (outdir / "report.html").write_text(report._html(df, meta, peptides), encoding="utf-8")
    print("report.html yazıldı.")
    from vaxforge import report_pdf
    print(f"{report_pdf.build(outdir, peptides, meta)} yazıldı.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
