"""Vaxign-ML eğitim seti ile benchmark ground-truth'umuz arasındaki örtüşmeyi ölç.

Vaxign-ML (Ong ve ark. 2020, doi:10.1093/bioinformatics/btaa119) dağıtılan modeli
Protegen'in TAMAMI üzerinde eğitilmiş XGBoost sınıflandırıcısıdır (pozitifler:
Protegen, 31 Temmuz 2019; negatifler: aynı 50 bakterinin UniProt proteomlarından
BPAg'lere <%30 benzer proteinler). Bizim head-to-head karşılaştırmamızın
ground-truth antijenleri de literatürden gelen deneysel koruyucu antijenler.

Soru: hangi türlerde ground-truth antijenimiz Vaxign-ML'in eğitim etiketinin
kendisi? Orada "başarı" genelleme değil, hatırlamadır.

Yöntem: Protegen v4.0 bakteri FASTA'sı (660 kayıt) indirilir; benchmark'taki 48
antijenin dizileri hazırlanmış proteomlardan çekilir; diamond blastp ile
karşılıklı en iyi hit ve %kimlik ölçülür. Eşikler: ≥%90 kimlik = pratikte aynı
protein (eğitim etiketi), %30-90 = homolog (Vaxign-ML negatif seçimi <%30
benzerlik kuralı kullandığı için bu aralık hâlâ bilgi sızdırır), <%30 = bağımsız.

Koşum:
    .venv/bin/python3 benchmark/measure_vaxignml_circularity.py \\
        --protegen /yol/protegen-bacterium-4.0-2019-01-09.faa
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "benchmark" / "results" / "vaxignml_circularity.json"
BACT = ROOT / "benchmark" / "bacteria"
PREP = ROOT / "data" / "validation"


def read_faa(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protegen", required=True)
    ap.add_argument("--diamond", default=shutil.which("diamond") or "diamond")
    args = ap.parse_args()

    protegen = Path(args.protegen)
    if not protegen.exists():
        print(f"Protegen FASTA yok: {protegen}")
        return

    # 1) benchmark ground-truth antijenleri + dizileri
    queries: dict[str, dict] = {}
    for cfg_path in sorted(BACT.glob("*.json")):
        bug = cfg_path.stem
        cfg = json.loads(cfg_path.read_text())
        faa = PREP / f"{bug}_prepared.faa"
        if not faa.exists():
            continue
        proteome = read_faa(faa)
        # config şeması: {"protective": {"<accession>": "<antijen adı>"}}
        for acc, label in (cfg.get("protective") or {}).items():
            seq = proteome.get(acc)
            if not seq:
                continue
            queries[f"{bug}|{acc}"] = {"bug": bug, "acc": acc,
                                       "antigen": label, "seq": seq}
    print(f"ground-truth antijen (dizisi bulunan): {len(queries)}")
    if not queries:
        print("antijen dizisi bulunamadı — hazırlanmış proteomlar eksik olabilir")
        return

    # 2) diamond: antijenler -> Protegen
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        qf = td / "query.faa"
        qf.write_text("".join(f">{k}\n{v['seq']}\n" for k, v in queries.items()))
        db = td / "protegen"
        subprocess.run([args.diamond, "makedb", "--in", str(protegen),
                        "-d", str(db), "--quiet"], check=True)
        tsv = td / "hits.tsv"
        subprocess.run([args.diamond, "blastp", "-q", str(qf), "-d", str(db),
                        "-o", str(tsv), "--outfmt", "6", "qseqid", "sseqid",
                        "pident", "length", "evalue", "--max-target-seqs", "1",
                        "--very-sensitive", "--quiet"], check=True)
        best: dict[str, tuple] = {}
        for ln in tsv.read_text().splitlines():
            q, s, pid, alen, ev = ln.split("\t")[:5]
            pid = float(pid)
            if q not in best or pid > best[q][1]:
                best[q] = (s, pid, int(alen), float(ev))

    # 3) sınıflandır
    rows = []
    for k, v in queries.items():
        hit = best.get(k)
        pid = hit[1] if hit else 0.0
        verdict = ("eğitim etiketinin kendisi" if pid >= 90 else
                   "homolog — bilgi sızdırır" if pid >= 30 else
                   "bağımsız")
        rows.append({
            "bug": v["bug"], "acc": v["acc"], "antigen": v["antigen"],
            "protegen_hit": hit[0] if hit else None,
            "pident": round(pid, 1), "align_len": hit[2] if hit else None,
            "evalue": hit[3] if hit else None, "verdict": verdict,
        })

    per_bug: dict[str, dict] = {}
    for r in rows:
        d = per_bug.setdefault(r["bug"], {"n": 0, "same": 0, "homolog": 0, "clean": 0})
        d["n"] += 1
        d["same" if r["pident"] >= 90 else
          "homolog" if r["pident"] >= 30 else "clean"] += 1

    rep = {
        "question": ("benchmark ground-truth antijenlerimiz Vaxign-ML'in eğitim "
                     "setinde (Protegen) var mı?"),
        "vaxignml": {"citation": ("Ong E, Wang H, Wong MU, Seetharaman M, Valdez N, "
                                 "He Y. Bioinformatics 2020;36(10):3185-3191"),
                     "doi": "10.1093/bioinformatics/btaa119",
                     "training_positives": "397 BPAg (584'ten homolog temizliği sonrası), Protegen 31-07-2019",
                     "training_negatives": "aynı 50 bakterinin UniProt proteomlarından BPAg'lere <%30 benzer proteinler",
                     "deployed_model": "tüm Protegen verisiyle eğitilmiş XGBoost"},
        "protegen": {"file": protegen.name, "n_records": protegen.read_text(
            encoding="utf-8", errors="replace").count(">"),
            "citation": "Yang B, Sayers S, Xiang Z, He Y. Nucleic Acids Res 2011;39(Database issue):D1073-8",
            "doi": "10.1093/nar/gkq944"},
        "thresholds": {"same_protein": 90, "homolog_leak": 30},
        "per_bug": per_bug,
        "antigens": sorted(rows, key=lambda r: (-r["pident"], r["bug"])),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"\n{'tür':12} {'n':>3} {'eğitimde':>9} {'homolog':>8} {'bağımsız':>9}")
    for bug, d in sorted(per_bug.items()):
        print(f"{bug:12} {d['n']:>3} {d['same']:>9} {d['homolog']:>8} {d['clean']:>9}")
    print(f"\n{'tür':11} {'antijen':34} {'kimlik':>7}  karar")
    for r in rep["antigens"]:
        if r["pident"] >= 30:
            print(f"{r['bug']:11} {r['antigen'][:34]:34} {r['pident']:>6}%  {r['verdict']}")
    print(f"\nyazıldı: {OUT}")


if __name__ == "__main__":
    main()
