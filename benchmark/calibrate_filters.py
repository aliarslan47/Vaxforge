"""Filtre kalibrasyonu — her kapı, bağımsız kontrol setleri üzerinde ölçülür.

NEDEN: 2026-10-10'da toksisite kapısının (ToxinPred2 Model 2) 9-mer epitoplarda
insan self peptitlerinin %74'ünü "toksik" işaretlediği, yani epitop aşamasının
rastgele bir 4'te-3 kırımı olduğu bulundu (ESAT-6/L7-L12/ActA/LLO kayıpları).
Kök neden araç değil, KULLANIM: tam-protein üzerinde eğitilmiş bir sınıflandırıcı
9-mer'e uygulanıyor ve üreticinin varsayılan eşiği hiç kalibre edilmemişti.

Aynı hatanın başka kapılarda da olup olmadığını ölçmek için bu betik TÜM
peptit-düzeyi kapıları ve protein-düzeyi SERT filtreleri ortak kontrol setlerinde
tarar. Çıktı: results/filter_calibration.json + FILTER_CALIBRATION.md.

KONTROL SETLERİ
  self9 / self15   : insan Swiss-Prot'tan rastgele 9/15-mer. Bir aşı epitobu
                     filtresi bunları ELEMEMELİ (üst sınır yanlış-pozitif oranı;
                     küçük bir kesrin gerçekten sitotoksik olması mümkün).
  iedb9 / iedb15   : tools/db/iedb önbelleğindeki DENEYSEL DOĞRULANMIŞ epitoplar
                     (patojen kaynaklı, insanda bağışık yanıt verdiği gösterilmiş).
                     Boru hattının bulması GEREKEN peptitler; elenmeleri doğrudan
                     duyarlılık kaybıdır.
  toxprot / nontox : ToxinPred2 datasets/*_realistic_dataset — TAM UZUNLUK
                     toksin/non-toksin proteinler. Yalnız "araç kendi uzunluk
                     sınıfında çalışıyor mu" sağlamasıdır.
                     !! ToxinPred2/3 için DAİRESEL (eğitim verisi) — duyarlılık
                     sayıları bu yüzden iyimser; raporda böyle etiketlenir.

EPİTOP UZUNLUĞUNDA "TOKSİK PEPTİT" İÇİN GERÇEK BİR ALTIN STANDART YOKTUR: bir
toksinden kesilmiş 9-mer'in kendisi toksik değildir (toksisite katlanmış proteinin
işlevidir). Bu yüzden burada epitop düzeyinde duyarlılık RAPORLANMAZ; ölçülen şey
yanlış-pozitif oranıdır ve varılan sonuç mimaridir: epitop düzeyinde toksisite
YUMUŞAK olmalı, sert eleme protein ve nihai MEV konstruktu düzeyinde yapılmalıdır.

Koşum:  .venv/bin/python3 benchmark/calibrate_filters.py [--n 1000]
"""

from __future__ import annotations

import argparse
import glob
import gzip
import json
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import allergen, mimicry, toxinpred, toxinpred3  # noqa: E402

RESULTS = ROOT / "benchmark" / "results"
SEED = 0
AA = set("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------
# kontrol setleri
# --------------------------------------------------------------------------
def _read_fasta_gz(path: Path) -> list[str]:
    seqs, cur = [], []
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(ln.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def _read_fasta(path: Path) -> list[str]:
    seqs, cur = [], []
    for ln in path.read_text(errors="ignore").splitlines():
        if ln.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(ln.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def self_kmers(length: int, n: int, rng: random.Random) -> list[str]:
    """İnsan Swiss-Prot'tan rastgele, çakışmasız k-mer örneği."""
    prots = _read_fasta_gz(ROOT / "tools" / "db" / "human_sprot.fasta.gz")
    out: set[str] = set()
    guard = 0
    while len(out) < n and guard < n * 200:
        guard += 1
        s = prots[rng.randrange(len(prots))]
        if len(s) < length + 2:
            continue
        i = rng.randrange(0, len(s) - length)
        p = s[i:i + length]
        if set(p) <= AA:
            out.add(p)
    return sorted(out)


def iedb_epitopes(lo: int, hi: int, n: int, rng: random.Random) -> list[str]:
    """Yerel IEDB önbelleğinden doğrulanmış epitoplar (uzunluk [lo, hi])."""
    out: set[str] = set()
    for f in sorted(glob.glob(str(ROOT / "tools" / "db" / "iedb" / "*.json"))):
        try:
            recs = json.load(open(f))
        except Exception:
            continue
        if not isinstance(recs, list):
            continue
        for r in recs:
            s = (r.get("seq") or "").strip().upper()
            if lo <= len(s) <= hi and set(s) <= AA:
                out.add(s)
    pool = sorted(out)
    rng.shuffle(pool)
    return sorted(pool[:n])


def tox_proteins(kind: str, n: int, rng: random.Random) -> list[str]:
    """ToxinPred2 realistic dataset'ten tam uzunluk protein örneği."""
    name = "Positive_realistic_dataset" if kind == "pos" else "Negative_realistic_dataset"
    p = ROOT / "tools" / "toxinpred2" / "datasets" / name
    if not p.exists():
        return []
    seqs = [s for s in _read_fasta(p) if s and set(s) <= AA]
    rng.shuffle(seqs)
    return seqs[:n]


# --------------------------------------------------------------------------
# kapılar — her biri {peptit: bool(elenir mi)} döndürür
# --------------------------------------------------------------------------
def gate_toxinpred2_hybrid(peps, thr=0.6):
    r = toxinpred.predict(peps, threshold=thr)
    return {p: bool(v["toxic"]) for p, v in r.items()}, \
           {p: v["ml_score"] for p, v in r.items()}


def gate_toxinpred3_m1(peps, thr=0.38):
    r = toxinpred3.predict(peps, threshold=thr, model=1)
    return {p: bool(v["toxic"]) for p, v in r.items()}, \
           {p: v["ml_score"] for p, v in r.items()}


def gate_toxinpred3_m2(peps, thr=0.38):
    r = toxinpred3.predict(peps, threshold=thr, model=2)
    return {p: bool(v["toxic"]) for p, v in r.items()}, \
           {p: v["ml_score"] for p, v in r.items()}


def gate_allergen(peps, **_):
    r = allergen.predict(peps)
    return {p: bool(v["allergen"]) for p, v in r.items()}, {}


def gate_mimicry_human(peps, **_):
    """Self-mimikri kapısı (YUMUŞAK; elemez, skoru düşürür) — işaretleme oranı."""
    try:
        r = mimicry.predict(peps, ["human"])
    except Exception:
        return {}, {}
    # self_match = tam eşleşme (sert mimikri); near_match = tek konservatif ikame
    return {p: bool(v.get("self_match") or v.get("near_match"))
            for p, v in (r or {}).items() if isinstance(v, dict)}, {}


GATES = {
    "toxicity_ToxinPred2_hybrid_t0.6": gate_toxinpred2_hybrid,
    "toxicity_ToxinPred3_ET_t0.38": gate_toxinpred3_m1,
    "toxicity_ToxinPred3_hybrid_t0.38": gate_toxinpred3_m2,
    "allergenicity_FAO_WHO_6mer": gate_allergen,
    "mimicry_human_exact_kmer(soft)": gate_mimicry_human,
}

# Hangi kontrol setinde elenmek YANLIŞ-POZİTİF sayılır
FP_SETS = ("self9", "self15", "iedb9", "iedb15", "nontox_protein")
SENS_SETS = ("tox_protein",)


def sweep_thresholds(fn, sets: dict, grid) -> list[dict]:
    """Eşik taraması: her eşikte FP oranları + protein duyarlılığı."""
    rows = []
    cache: dict[str, dict] = {}
    for name, peps in sets.items():
        _, scores = fn(peps)
        cache[name] = scores
    for t in grid:
        row = {"threshold": round(t, 3)}
        for name, scores in cache.items():
            if not scores:
                continue
            flag = sum(1 for v in scores.values() if v >= t)
            row[name] = round(100 * flag / len(scores), 1)
        rows.append(row)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000, help="set başına peptit sayısı")
    ap.add_argument("--nprot", type=int, default=200, help="protein seti boyutu")
    args = ap.parse_args()
    rng = random.Random(SEED)

    print("kontrol setleri kuruluyor…", flush=True)
    sets = {
        "self9": self_kmers(9, args.n, rng),
        "self15": self_kmers(15, args.n, rng),
        "iedb9": iedb_epitopes(9, 9, args.n, rng),
        "iedb15": iedb_epitopes(13, 17, args.n, rng),
        "tox_protein": tox_proteins("pos", args.nprot, rng),
        "nontox_protein": tox_proteins("neg", args.nprot, rng),
    }
    for k, v in sets.items():
        print(f"  {k:16s} n={len(v)}", flush=True)

    report: dict = {"seed": SEED, "set_sizes": {k: len(v) for k, v in sets.items()},
                    "gates": {}, "sweeps": {}}

    for gname, fn in GATES.items():
        print(f"\n[{gname}]", flush=True)
        g: dict = {}
        for sname, peps in sets.items():
            if not peps:
                continue
            flags, scores = fn(peps)
            if not flags:
                g[sname] = {"n": len(peps), "note": "araç sonuç vermedi"}
                print(f"  {sname:16s} — sonuç yok", flush=True)
                continue
            n = len(flags)
            k = sum(1 for v in flags.values() if v)
            entry = {"n": n, "flagged": k, "flagged_pct": round(100 * k / n, 1)}
            if scores:
                sv = list(scores.values())
                entry["score_median"] = round(statistics.median(sv), 3)
                entry["score_mean"] = round(statistics.mean(sv), 3)
            role = "YANLIŞ-POZİTİF" if sname in FP_SETS else (
                "duyarlılık" if sname in SENS_SETS else "")
            g[sname] = entry
            print(f"  {sname:16s} n={n:5d} işaretlenen=%{entry['flagged_pct']:<5} {role}",
                  flush=True)
        report["gates"][gname] = g

    # toksisite araçları için eşik taraması (hangi eşik FP<=%5 verir?)
    grid = [i / 100 for i in range(5, 100, 5)]
    for gname, fn in (("ToxinPred2_hybrid", gate_toxinpred2_hybrid),
                      ("ToxinPred3_ET", gate_toxinpred3_m1)):
        print(f"\n[eşik taraması] {gname}", flush=True)
        report["sweeps"][gname] = sweep_thresholds(fn, sets, grid)

    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / "filter_calibration.json"
    out.write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"\nyazıldı: {out}")


if __name__ == "__main__":
    main()
