"""VaxForge virüs kolu — Rabies (NC_001542) × 4 konak (insan/sığır/fare/domuz).

VİRÜS KOLU'nun görevi recall değil, ÇOK-KONAK MHC yeteneğinin gösterimi: aynı virüs
proteomunda host'a göre epitop repertuarının nasıl değiştiğini ölçer. Ground-truth =
glycoprotein G (NP_056796.1). NERVE/Vaxign host seçimi yapamaz → bu kol yalnız VaxForge.

İstatistik ([[benchmark/stats.py]]): G rank (permütasyon, n=1 caveat) + host-çifti Jaccard
örtüşmesi + host başına bağlanan-epitop sayısı + IEDB popülasyon kapsamı.

Çalıştırma:  . .venv/bin/activate && python3 benchmark/run_rabies.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmark"))

import stats as bstats  # noqa: E402

from vaxforge import pipeline  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402

# Env-parametrize: aynı çok-konak virüs kolu, farklı virüsler için (VSV/EMCV vb.).
# Varsayılan = Rabies (geriye dönük uyumlu).
PROTEOME = Path(os.environ.get("VIRUS_PROTEOME", ROOT / "data" / "validation" / "rabies_proteome.faa"))
TAXON = os.environ.get("VIRUS_TAXON", "NCBITaxon:11292")          # Lyssavirus rabies
HOSTS = ["human", "bovine", "mouse", "pig"]
GROUND_TRUTH = os.environ.get("VIRUS_GT", "NP_056796")            # G glikoprotein (ID substring)
GT_NAME = os.environ.get("VIRUS_GT_NAME", "glycoprotein G")
GT_ACC = os.environ.get("VIRUS_GT_ACC", "NP_056796.1")
VIRUS_LABEL = os.environ.get("VIRUS_LABEL", "rabies_proteome.faa (NC_001542, 5 protein)")
VIRUS_NAME = os.environ.get("VIRUS_NAME", "Rabies")
OUTDIR = Path(os.environ.get("VIRUS_OUTDIR", ROOT / "benchmark" / "results" / "rabies"))


def main() -> int:
    if not PROTEOME.exists():
        print(f"HATA: proteom yok: {PROTEOME}")
        return 1
    cfg = ThresholdConfig.load()
    hosts = HostRegistry.load()
    det = detect(str(PROTEOME))
    det.filename = VIRUS_LABEL
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print(f"{VIRUS_NAME} × {len(HOSTS)} konak: {', '.join(HOSTS)} · profil=virus · taxon={TAXON}")
    t0 = time.time()
    result = None
    for ev in pipeline.run(str(PROTEOME), det, cfg, "virus",
                           host_names=HOSTS, has_gpu=False,
                           outdir=str(OUTDIR.parent), host_registry=hosts,
                           organism_taxon=TAXON):
        ph, stt, msg = ev["phase"], ev["status"], ev["msg"]
        if ph == "__result__":
            result = ev["data"]
        elif ph == "__error__":
            print(f"  ✗ HATA: {msg}")
            return 1
        else:
            print(f"  [{time.time()-t0:6.0f}s] {stt:9s} {ph:16s} {msg}")

    if not result:
        print("Sonuç yok.")
        return 1
    peptides, meta = result["peptides"], result["meta"]

    # --- protein-seviyesi aday sıralaması + G rank ---
    cand = {}
    for p in peptides:
        cand[p.parent] = max(cand.get(p.parent, 0.0), p.candidacy)
    ranked = sorted(cand, key=lambda a: -cand[a])
    rank_of = {a: i + 1 for i, a in enumerate(ranked)}
    g_id = next((a for a in cand if GROUND_TRUTH in a), None)
    g_rank = rank_of.get(g_id)
    g_produced = g_id is not None

    # G rank permütasyon (n=1 — betimsel; scores = tüm aday proteinler)
    perm = None
    if g_produced and len(cand) >= 2:
        scores = [cand[a] for a in ranked]
        perm = bstats.permutation_rank_test(scores, [ranked.index(g_id)],
                                            n_perm=10000, statistic="best_rank").as_dict()

    # --- ÇOK-KONAK epitop repertuarı ---
    # host -> set(peptit seq)  (MHC-I + MHC-II, o host'ta sunulan)
    host_eps: dict[str, set] = {h: set() for h in HOSTS}
    host_eps_G: dict[str, set] = {h: set() for h in HOSTS}
    for p in peptides:
        if p.kind not in ("MHC-I", "MHC-II"):
            continue
        presented = p.metrics.get("hosts_presented", []) or []
        on_G = g_id is not None and GROUND_TRUTH in p.parent
        for h in presented:
            if h in host_eps:
                host_eps[h].add(p.seq)
                if on_G:
                    host_eps_G[h].add(p.seq)

    host_binder_counts = {h: len(host_eps[h]) for h in HOSTS}
    host_binder_counts_G = {h: len(host_eps_G[h]) for h in HOSTS}

    # host-çifti Jaccard (tüm proteom + yalnız G)
    jacc = {}
    jacc_G = {}
    for a, b in combinations(HOSTS, 2):
        jacc[f"{a}~{b}"] = bstats.jaccard(host_eps[a], host_eps[b]).as_dict()
        jacc_G[f"{a}~{b}"] = bstats.jaccard(host_eps_G[a], host_eps_G[b]).as_dict()

    popcov = meta.get("population_coverage")

    summary = {
        "input": det.filename, "taxon": TAXON, "hosts": HOSTS,
        "runtime_s": round(time.time() - t0, 1),
        "n_input": meta.get("n_input"), "n_funnel": meta.get("n_funnel"),
        "n_candidate_proteins": len(cand),
        "ground_truth": {"antigen": GT_NAME, "accession": GT_ACC,
                         "matched_id": g_id, "produced_candidate": g_produced,
                         "protein_rank": g_rank, "best_candidacy": round(cand.get(g_id, 0.0), 4) if g_produced else None},
        "g_rank_permutation": perm,
        "multihost": {
            "binder_counts_all_proteins": host_binder_counts,
            "binder_counts_G_only": host_binder_counts_G,
            "jaccard_overlap_all": jacc,
            "jaccard_overlap_G": jacc_G,
        },
        "population_coverage": popcov,
        "report_paths": {k: str(v) for k, v in result["paths"].items()},
        "note": ("Virüs kolu = çok-konak MHC yetenek gösterimi (recall değil). "
                 "Tek antijen (G) → rank istatistiği betimsel. Host seçimi VaxForge'a özgü; "
                 "NERVE (bacterial-only) ve Vaxign bu kolda oynayamaz."),
    }
    (OUTDIR / f"{VIRUS_NAME.lower()}_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    import shutil
    for k, v in result["paths"].items():
        try:
            shutil.copy(v, OUTDIR / Path(v).name)
        except Exception:
            pass

    print("\n" + "=" * 70)
    print(f"VİRÜS KOLU ÖZETİ — {VIRUS_NAME} × 4 konak")
    print("=" * 70)
    print(f"Girdi 5 protein · funnel sonrası {summary['n_funnel']} · aday protein "
          f"{summary['n_candidate_proteins']} · süre {summary['runtime_s']}s")
    gt = summary["ground_truth"]
    mark = "✅" if gt["produced_candidate"] else "❌"
    print(f"\n{mark} Ground-truth {GT_NAME}: rank #{gt['protein_rank']}/"
          f"{summary['n_candidate_proteins']} · skor {gt['best_candidacy']}")
    if perm:
        print(f"   G rank permütasyon p={perm['p_value']:.4g} ({perm['caveat'] or 'n=1 betimsel'})")
    print(f"\nHost başına bağlanan epitop (tüm proteom): {host_binder_counts}")
    print(f"Host başına bağlanan epitop (yalnız G):     {host_binder_counts_G}")
    print(f"\nHost-çifti Jaccard örtüşmesi (G epitopları):")
    for k, v in jacc_G.items():
        print(f"   {k:18s} J={v['statistic']}  (ortak {v['detail']['intersection']}/{v['detail']['union']})")
    print(f"\nÇıktılar: {OUTDIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
