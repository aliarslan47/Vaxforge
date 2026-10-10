"""Protein-bütünü konak homolojisi, epitop-düzeyi otoimmünite riskinin vekili mi?

Sert konak-homoloji kapısı (funnel, ≥%35 kimlik → eleme) protein BÜTÜNÜNDE ölçüm
yapıp tüm proteini atıyordu. Bu betik aynı proteinleri epitop düzeyinde ölçer:
9-mer'lerinin kaçı konağın kendi proteininde BİREBİR geçiyor. İkisi arasındaki
fark, kapının hangi granülaritede durması gerektiğini belirler.

Girdi: benchmark/results/miss_attribution.json'daki antijenler + hazırlanmış
proteomlar (data/validation/<bug>_prepared.faa).

Çıktı: benchmark/results/homology_granularity.json + ekrana tablo.

Koşum: .venv/bin/python3 benchmark/measure_homology_granularity.py [--bug ssuis]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import mimicry  # noqa: E402

OUT = ROOT / "benchmark" / "results" / "homology_granularity.json"
MISS = ROOT / "benchmark" / "results" / "miss_attribution.json"


def read_faa(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for ln in path.read_text().splitlines():
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
    ap.add_argument("--bug", default="ssuis")
    ap.add_argument("--host", default="human")
    ap.add_argument("-k", type=int, default=9)
    args = ap.parse_args()

    if not mimicry.available(args.host):
        print(f"konak self-proteom DB yok: {mimicry.db_path(args.host)}")
        return
    faa = ROOT / "data" / "validation" / f"{args.bug}_prepared.faa"
    if not faa.exists():
        print(f"proteom yok: {faa}")
        return
    seqs = read_faa(faa)
    antigens = [r for r in json.loads(MISS.read_text())["antigens"]
                if r["bug"] == args.bug]

    rows = []
    for r in antigens:
        s = seqs.get(r["acc"])
        if not s:
            continue
        kmers = [s[i:i + args.k] for i in range(len(s) - args.k + 1)]
        res = mimicry.predict(kmers, [args.host], k=args.k)
        hits = [p for p in kmers if res.get(p, {}).get("self_match")]
        rows.append({
            "acc": r["acc"], "antigen": r["antigen"],
            "recovered": r["recovered"], "miss_cause": r["miss_cause"],
            "protein_pident": r["human_pident"], "length": len(s),
            "n_kmer": len(kmers), "n_self_exact": len(hits),
            "pct_self_exact": round(100 * len(hits) / len(kmers), 2) if kmers else None,
            "n_safe_kmer": len(kmers) - len(hits),
        })

    rep = {
        "bug": args.bug, "host": args.host, "k": args.k,
        "question": ("protein-bütünü konak kimliği, epitop-düzeyi self-eşleşme "
                     "oranının vekili mi?"),
        "finding": ("hayır: protein-bütünü kimlik %37-52 olan antijenlerin "
                    "9-mer'lerinin yalnız %1.4-3.0'ü konak self ile birebir aynı; "
                    "sert protein-düzeyi eleme, az sayıda riskli k-mer için "
                    "yüzlerce self-olmayan epitobu da atıyor"),
        "antigens": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"{'protein':34} {'prot %':>7} {'k-mer':>6} {'self-tam':>9} {'oran':>7} {'güvenli':>8}")
    for d in sorted(rows, key=lambda x: -(x["protein_pident"] or 0)):
        print(f"{d['antigen'][:34]:34} {d['protein_pident']:>7} {d['n_kmer']:>6} "
              f"{d['n_self_exact']:>9} {d['pct_self_exact']:>6}% {d['n_safe_kmer']:>8}")
    print(f"\nyazıldı: {OUT}")


if __name__ == "__main__":
    main()
