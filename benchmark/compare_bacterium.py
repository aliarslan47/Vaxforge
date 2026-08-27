"""Bakteri 3'lü head-to-head: VaxForge vs NERVE 2.0 vs Vaxign-ML (TAM PROTEOM, config-sürücülü).

compare_menb.py'nin genellemesi: 5 koruyucu antijen, config'ten. Her araç DEFAULT parametre
(adil), özdeş girdi (çıplak-aksesyon tam proteom). Çıktı: results/<bug>_headtohead.json.

Çalıştırma:  python3 benchmark/compare_bacterium.py <saureus|listeria|salmonella>
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "benchmark"))
import stats as bstats  # noqa: E402

BENCH = ROOT / "benchmark"
ALL_ACC_RE = re.compile(r"^([A-Z][0-9A-Z]{5}) ")


def _short(label: str) -> str:
    """Antijen etiketinden kısa ad (ilk token, parantez öncesi)."""
    return label.split(" (")[0].split()[0]


def _rank_map(scored):
    scored = sorted(scored, key=lambda t: -t[1])
    return {a: i + 1 for i, (a, _) in enumerate(scored)}


def load(bug: str):
    cfg = json.loads((BENCH / "bacteria" / f"{bug}.json").read_text())
    protective = cfg["protective"]           # acc -> uzun etiket
    accs = set(protective)
    short = {a: _short(protective[a]) for a in accs}
    n_prot = len(accs)
    n_bg = cfg["n_proteins"] - n_prot
    res = BENCH / "results" / bug
    return cfg, protective, accs, short, n_prot, n_bg, res


def vaxforge(res, accs, short, n_prot, n_bg):
    s = json.load(open(res / "vaxforge_summary.json"))
    rows = {r["accession"]: r for r in s["protective_detail"]}
    ranks = {a: rows[a]["protein_rank"] for a in accs}
    recovered = [a for a in accs if rows[a]["produced_candidate"]]
    bg_cand = s["background_producing_candidates"]
    fish = bstats.fisher_enrichment(len(recovered), n_prot, bg_cand, n_bg)
    auc = None
    if len(recovered) == n_prot:
        prod_ranks = sorted(ranks[a] for a in recovered)
        concordant = 0
        for r in prod_ranks:
            pos_above = sum(1 for rr in prod_ranks if rr < r)
            cand_neg_above = (r - 1) - pos_above
            concordant += (n_bg - cand_neg_above)
        auc = round(concordant / (n_prot * n_bg), 3)
    return {"tool": "VaxForge", "n_candidates": s["n_candidate_proteins"],
            "recall": f"{len(recovered)}/{n_prot}",
            "ranks": {short[a]: ranks[a] for a in accs},
            "recovered": sorted(short[a] for a in recovered),
            "fold_enrichment": s["fold_enrichment"],
            "fisher_p": fish.p_value, "auc": auc,
            "note": "candidacy skoru; funnel soft-VFDB"}


def vaxignml(res, bug, accs, short, n_prot, n_bg):
    # Vaxign-ML çıktısı girdi-adına göre isimlenir (bazı bakterilerde sanitize edilmiş
    # ayrı dosya kullanıldı) → tek *.result.tsv'yi glob ile bul (sağlam).
    tsvs = sorted((res / "vaxignml").glob("*.result.tsv"))
    if not tsvs:
        raise FileNotFoundError(f"{res}/vaxignml içinde *.result.tsv yok")
    f = tsvs[0]
    scored, preds = [], {}
    with open(f) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            acc = x["sample"].split()[0]
            scored.append((acc, float(x["protegenicity"])))
            preds[acc] = float(x["prediction"])
    rank = _rank_map(scored)
    n_pred = sum(1 for a in preds if preds[a] >= 0.5)
    recovered = [a for a in accs if preds.get(a, 0) >= 0.5]
    y = [1 if a in accs else 0 for a, _ in scored]
    sc = [v for _, v in scored]
    auc = bstats.roc_auc_ci(y, sc)
    prot_pred = sum(1 for a in accs if preds.get(a, 0) >= 0.5)
    bg_pred = n_pred - prot_pred
    fish = bstats.fisher_enrichment(prot_pred, n_prot, bg_pred, n_bg)
    return {"tool": "Vaxign-ML", "n_candidates": n_pred,
            "recall": f"{len(recovered)}/{n_prot}",
            "ranks": {short[a]: rank.get(a) for a in accs},
            "recovered": sorted(short[a] for a in recovered),
            "fold_enrichment": fish.detail["fold_enrichment"],
            "fisher_p": fish.p_value, "auc": auc.statistic, "auc_ci": auc.ci,
            "note": "protegenicity (XGBoost); pred>=0.5 protective"}


def nerve(res, accs, short, n_prot, n_bg):
    cand_f = res / "nerve" / "vaccine_candidates.csv"
    disc_f = res / "nerve" / "discarded_proteins.csv"
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
    recovered = [a for a in accs if a in cand]
    prot_cand = len(accs & set(cand))
    bg_cand = n_cand - prot_cand
    fish = bstats.fisher_enrichment(prot_cand, n_prot, bg_cand, n_bg)
    disc_reason = {}
    for line in open(disc_f):
        m = ALL_ACC_RE.match(line)
        if m and m.group(1) in accs:
            cols = line.split(",")
            loc = cols[5] if len(cols) > 5 else "?"
            adh = cols[8] if len(cols) > 8 else "?"
            disc_reason[short[m.group(1)]] = f"loc={loc} adhesin_p={adh}"
    return {"tool": "NERVE 2.0", "n_candidates": n_cand,
            "recall": f"{len(recovered)}/{n_prot}",
            "ranks": {short[a]: (rank.get(a) or "DISCARDED") for a in accs},
            "recovered": sorted(short[a] for a in recovered),
            "fold_enrichment": fish.detail["fold_enrichment"],
            "fisher_p": fish.p_value, "auc": None,
            "discarded_reasons": disc_reason,
            "note": "select modülü sert filtreler (localization + adhesin pad>=0.5)"}


def main() -> int:
    if len(sys.argv) != 2:
        print("kullanım: python3 benchmark/compare_bacterium.py <bug>")
        return 2
    bug = sys.argv[1]
    cfg, protective, accs, short, n_prot, n_bg, res = load(bug)
    order = list(accs)  # sabit sütun sırası
    tools = [vaxforge(res, accs, short, n_prot, n_bg),
             vaxignml(res, bug, accs, short, n_prot, n_bg),
             nerve(res, accs, short, n_prot, n_bg)]
    out = {"pathogen": f"{cfg['organism']} {cfg['strain']}",
           "input": f"{bug}_prepared.faa (TAM PROTEOM, {cfg['n_proteins']} protein)",
           "gram": cfg["gram"], "genome_mb": cfg["genome_mb"],
           "fairness": "her araç DEFAULT parametre; eşik tuning yok; özdeş tam-proteom girdi",
           "ground_truth": protective, "n_protective": n_prot, "n_background": n_bg,
           "tools": tools}
    out_json = res.parent / f"{bug}_headtohead.json"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    cols = [short[a] for a in order]
    print("=" * 92)
    print(f"{cfg['organism']} HEAD-TO-HEAD — VaxForge vs NERVE 2.0 vs Vaxign-ML (TAM PROTEOM, default)")
    print("=" * 92)
    hdr = f"{'Araç':11s} {'recall':7s} {'aday':6s} " + " ".join(f"{c:>7s}" for c in cols) + f" {'fold':>6s} {'Fisher p':>9s} {'AUC':>5s}"
    print(hdr); print("-" * len(hdr))
    for t in tools:
        r = t["ranks"]
        fp = f"{t['fisher_p']:.3g}" if t["fisher_p"] is not None else "-"
        au = f"{t['auc']:.2f}" if t["auc"] is not None else "-"
        fo = f"{t['fold_enrichment']:.2f}" if t["fold_enrichment"] is not None else "-"
        cells = " ".join(f"{str(r[c]):>7s}" for c in cols)
        print(f"{t['tool']:11s} {t['recall']:7s} {t['n_candidates']:6d} {cells} {fo:>6s} {fp:>9s} {au:>5s}")
    if tools[2].get("discarded_reasons"):
        print("\nNERVE elenen koruyucu antijenler (sert filtre):")
        for k, v in tools[2]["discarded_reasons"].items():
            print(f"   {k}: {v}")
    print(f"\n→ {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
