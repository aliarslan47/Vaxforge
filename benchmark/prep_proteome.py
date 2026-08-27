"""Bakteri proteomunu çıplak-aksesyon başlıklı dosyaya hazırla.

Tam UniProt proteomu `>sp|ACC|NAME ...` başlıklarıyla gelir; VaxForge/NERVE/Vaxign-ML
üçüne de AYNI çıplak-aksesyon başlıklı FASTA verilirse (MenB testset kalıbı) aksesyon
eşleştirmesi (recall@antijen) sorunsuz olur. Çıktı: data/validation/<bug>_prepared.faa
(ilk token = çıplak aksesyon).

Kullanım:  python3 benchmark/prep_proteome.py <bug>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent


def _acc(rec_id: str) -> str:
    return rec_id.split("|")[1] if "|" in rec_id else rec_id


def prepare(bug: str) -> Path:
    cfg = json.loads((ROOT / "benchmark" / "bacteria" / f"{bug}.json").read_text())
    src = ROOT / cfg["proteome"]
    out = ROOT / "data" / "validation" / f"{bug}_prepared.faa"
    n = 0
    with out.open("w") as fh:
        for rec in SeqIO.parse(str(src), "fasta"):
            acc = _acc(rec.id)
            desc = rec.description.split(None, 1)
            tail = desc[1][:70] if len(desc) > 1 else ""
            fh.write(f">{acc} {tail}\n{str(rec.seq)}\n")
            n += 1
    # koruyucu antijenlerin hazırlanan dosyada bulunduğunu doğrula
    accs = {r.id for r in SeqIO.parse(str(out), "fasta")}
    missing = [a for a in cfg["protective"] if a not in accs]
    print(f"{bug}: {n} protein -> {out.name} | koruyucu {len(cfg['protective'])-len(missing)}/{len(cfg['protective'])} bulundu")
    if missing:
        print(f"  ⚠️ EKSİK koruyucu aksesyon: {missing}")
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("kullanım: python3 benchmark/prep_proteome.py <saureus|listeria|salmonella>")
        raise SystemExit(2)
    prepare(sys.argv[1])
