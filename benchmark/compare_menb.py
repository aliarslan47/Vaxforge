"""MenB 3'lü head-to-head: VaxForge vs NERVE 2.0 vs Vaxign-ML.

Aynı 64-protein test seti (4 koruyucu 4CMenB antijeni + 60 arka plan), HER ARAÇ DEFAULT
parametrelerle (adil koşul — eşik tuning yok). Üç aracın çıktısını okur, birleşik tablo +
istatistik (Fisher enrichment + ROC-AUC/bootstrap) üretir → results/menb_headtohead.json.

Çalıştırma:  . .venv/bin/activate && python3 benchmark/compare_menb.py
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "benchmark"))
import stats as bstats  # noqa: E402

BENCH = ROOT / "benchmark"
PROTECTIVE = {"Q9JXV4": "fHbp", "Q9JXK7": "NadA", "Q7DD37": "NHBA(VFDB-)", "P0DH58": "PorA"}
ACCS = set(PROTECTIVE)

# Env ile set seçimi (64 varsayılan / 204 flagship).
VXF_SUMMARY = Path(os.environ.get("VXF_SUMMARY", ROOT / "outputs" / "validation_menb" / "validation_summary.json"))
VXG_TSV = Path(os.environ.get("VXG_TSV", BENCH / "results" / "vaxignml_menb" / "menb_testset.result.tsv"))
NERVE_DIR = Path(os.environ.get("NERVE_DIR", BENCH / "tools" / "NERVE" / "out_menb"))
N_BG = int(os.environ.get("N_BG", "60"))
OUT_JSON = Path(os.environ.get("OUT_JSON", BENCH / "results" / "menb_headtohead.json"))
ALL_ACC_RE = re.compile(r"^([A-Z][0-9A-Z]{5}) ")


def _rank_map(scored: list[tuple[str, float]]) -> dict[str, int]:
    scored = sorted(scored, key=lambda t: -t[1])
    return {a: i + 1 for i, (a, _) in enumerate(scored)}


def vaxforge() -> dict:
    s = json.load(open(VXF_SUMMARY))
    rows = {r["accession"]: r for r in s["protective_detail"]}
    n_cand = s["n_candidate_proteins"]
    ranks = {a: rows[a]["protein_rank"] for a in ACCS}
    recovered = [a for a in ACCS if rows[a]["produced_candidate"]]
    bg_cand = s["background_producing_candidates"]
    fish = bstats.fisher_enrichment(len(recovered), 4, bg_cand, N_BG)
    # AUC: candidacy sıralamasından tam concordance. Adaylar tüm arka planın üstünde;
    # aday-olmayan arka plan candidacy=0 (tüm koruyucu adayların altında).
    auc = None
    if len(recovered) == 4:  # hepsi aday → temiz hesap
        prod_ranks = sorted(ranks[a] for a in recovered)  # aday-içi rank (1..n_cand)
        concordant = 0
        for r in prod_ranks:
            pos_above = sum(1 for rr in prod_ranks if rr < r)
            cand_neg_above = (r - 1) - pos_above         # üstteki arka-plan adayları
            concordant += (N_BG - cand_neg_above)         # kalan tüm negatifleri geçer
        auc = round(concordant / (4 * N_BG), 3)
    return {"tool": "VaxForge", "n_candidates": n_cand, "recall": f"{len(recovered)}/4",
            "ranks": {PROTECTIVE[a]: ranks[a] for a in ACCS},
            "recovered": sorted(PROTECTIVE[a] for a in recovered),
            "fold_enrichment": s["fold_enrichment"],
            "fisher_p": fish.p_value, "auc": auc,
            "note": "candidacy skoru; funnel soft-VFDB (NHBA korunur)"}


def vaxignml() -> dict:
    f = VXG_TSV
    scored, preds = [], {}
    with open(f) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            acc = x["sample"].split()[0]
            scored.append((acc, float(x["protegenicity"])))
            preds[acc] = float(x["prediction"])
    rank = _rank_map(scored)
    n_pred = sum(1 for a in preds if preds[a] >= 0.5)
    recovered = [a for a in ACCS if preds.get(a, 0) >= 0.5]
    # istatistik: protegenicity skoru sınıflandırıcı
    y = [1 if a in ACCS else 0 for a, _ in scored]
    sc = [s for _, s in scored]
    auc = bstats.roc_auc_ci(y, sc)
    prot_pred = sum(1 for a in ACCS if preds.get(a, 0) >= 0.5)
    bg_pred = n_pred - prot_pred
    fish = bstats.fisher_enrichment(prot_pred, 4, bg_pred, N_BG)
    return {"tool": "Vaxign-ML", "n_candidates": n_pred, "recall": f"{len(recovered)}/4",
            "ranks": {PROTECTIVE[a]: rank.get(a) for a in ACCS},
            "recovered": sorted(PROTECTIVE[a] for a in recovered),
            "fold_enrichment": fish.detail["fold_enrichment"],
            "fisher_p": fish.p_value, "auc": auc.statistic, "auc_ci": auc.ci,
            "note": "protegenicity (XGBoost); pred>=0.5 protective"}


def nerve() -> dict:
    cand_f = NERVE_DIR / "vaccine_candidates.csv"
    disc_f = NERVE_DIR / "discarded_proteins.csv"
    # aday accession + score (gömülü newline'lara karşı satır-başı accession ile)
    cand = {}
    for line in open(cand_f):
        m = ALL_ACC_RE.match(line)
        if m:
            parts = line.split(",")
            sc = 0.0
            for p in parts[2:5]:
                try:
                    sc = float(p); break
                except ValueError:
                    continue
            cand[m.group(1)] = sc
    n_cand = len(cand)
    rank = _rank_map(list(cand.items()))
    recovered = [a for a in ACCS if a in cand]
    # Fisher: koruyucu aday / arka plan aday
    prot_cand = len(ACCS & set(cand))
    bg_cand = n_cand - prot_cand
    fish = bstats.fisher_enrichment(prot_cand, 4, bg_cand, N_BG)
    # discard sebepleri (koruyucular için)
    disc_reason = {}
    for line in open(disc_f):
        m = ALL_ACC_RE.match(line)
        if m and m.group(1) in ACCS:
            cols = line.split(",")
            loc = cols[5] if len(cols) > 5 else "?"
            adh = cols[8] if len(cols) > 8 else "?"
            disc_reason[PROTECTIVE[m.group(1)]] = f"loc={loc} adhesin_p={adh}"
    return {"tool": "NERVE 2.0", "n_candidates": n_cand, "recall": f"{len(recovered)}/4",
            "ranks": {PROTECTIVE[a]: (rank.get(a) or "DISCARDED") for a in ACCS},
            "recovered": sorted(PROTECTIVE[a] for a in recovered),
            "fold_enrichment": fish.detail["fold_enrichment"],
            "fisher_p": fish.p_value, "auc": None,
            "discarded_reasons": disc_reason,
            "note": "select modülü sert filtreler (localization + adhesin pad>=0.5)"}


def main() -> int:
    tools = [vaxforge(), vaxignml(), nerve()]
    out = {"pathogen": "MenB (N. meningitidis MC58)", "input": "menb_testset.faa (4 koruyucu + 60 arka plan = 64)",
           "fairness": "her araç DEFAULT parametre; eşik tuning yok; özdeş girdi",
           "ground_truth": PROTECTIVE, "tools": tools}
    OUT_JSON.write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    print("=" * 78)
    print("MenB HEAD-TO-HEAD — VaxForge vs NERVE 2.0 vs Vaxign-ML (default, özdeş girdi)")
    print("=" * 78)
    hdr = f"{'Araç':12s} {'recall':7s} {'aday':5s} {'fHbp':>6s} {'NadA':>6s} {'NHBA':>6s} {'PorA':>6s} {'fold':>6s} {'Fisher p':>9s} {'AUC':>5s}"
    print(hdr); print("-" * len(hdr))
    for t in tools:
        r = t["ranks"]
        fp = f"{t['fisher_p']:.3g}" if t["fisher_p"] is not None else "-"
        au = f"{t['auc']:.2f}" if t["auc"] is not None else "-"
        fo = f"{t['fold_enrichment']:.2f}" if t["fold_enrichment"] is not None else "-"
        print(f"{t['tool']:12s} {t['recall']:7s} {t['n_candidates']:5d} "
              f"{str(r['fHbp']):>6s} {str(r['NadA']):>6s} {str(r['NHBA(VFDB-)']):>6s} {str(r['PorA']):>6s} "
              f"{fo:>6s} {fp:>9s} {au:>5s}")
    print("\nNERVE elenen koruyucu antijenler (sert filtre):")
    for k, v in tools[2].get("discarded_reasons", {}).items():
        print(f"   {k}: {v}")
    print(f"\n→ {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
