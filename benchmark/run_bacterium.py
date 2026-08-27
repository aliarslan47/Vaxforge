"""VaxForge bakteri kolu — TAM PROTEOM koşusu (config-sürücülü, MenB validate_menb genellemesi).

MenB'den farkı: downsample YOK — tüm proteom pipeline'dan geçer. Koruyucu antijenler
config'ten (5 adet). Arka plan = kalan tüm proteom. Ölçütler MenB ile aynı schema:
recall@antijen + antijen-başına protein_rank + fold-enrichment (koruyucu vs tüm proteom) +
Fisher + permütasyon-rank + AUC(caveat) + per-antijen VFDB durumu + IEDB.

Çalıştırma:  BUG=saureus . .venv/bin/activate && python3 benchmark/run_bacterium.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmark"))

import stats as bstats  # noqa: E402
from prep_proteome import prepare  # noqa: E402

from vaxforge import discovery, ingest, pipeline  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402

BUG = os.environ.get("BUG")
if not BUG:
    print("HATA: BUG env yok (saureus|listeria|salmonella)")
    raise SystemExit(2)
CFG = json.loads((ROOT / "benchmark" / "bacteria" / f"{BUG}.json").read_text())
PROTECTIVE = CFG["protective"]
GRAM = CFG["gram"]                       # "positive" | "negative"
TAXON = CFG["iedb_taxon"]                # NCBITaxon:<species>
OUTDIR = Path(os.environ.get("OUTDIR", ROOT / "outputs" / f"validation_{BUG}"))
RESULTS = ROOT / "benchmark" / "results" / BUG


def _smoke_subset(full: Path, n_bg: int) -> Path:
    """SMOKE_N testi: koruyucu + n_bg arka plan ile küçük alt-küme (kod yolunu hızlı doğrula)."""
    from Bio import SeqIO
    recs = list(SeqIO.parse(str(full), "fasta"))
    prot = [r for r in recs if r.id in PROTECTIVE]
    bg = [r for r in recs if r.id not in PROTECTIVE][:n_bg]
    out = full.with_name(f"{BUG}_smoke.faa")
    SeqIO.write(prot + bg, str(out), "fasta")
    print(f"[SMOKE] {len(prot)} koruyucu + {len(bg)} arka plan = {len(prot)+len(bg)} protein -> {out.name}")
    return out


def main() -> int:
    protective_accs = set(PROTECTIVE)
    prepared = ROOT / "data" / "validation" / f"{BUG}_prepared.faa"
    if os.environ.get("SKIP_PREP") and prepared.exists():
        proteome = prepared  # orkestratör önceden hazırladı (eşzamanlı koşuda yarışı önle)
        print(f"{BUG}: hazırlanmış proteom kullanılıyor (SKIP_PREP) -> {prepared.name}")
    else:
        proteome = prepare(BUG)  # çıplak-aksesyon başlıklı tam proteom
    smoke_n = os.environ.get("SMOKE_N")
    if smoke_n:
        proteome = _smoke_subset(proteome, int(smoke_n))

    cfg = ThresholdConfig.load()
    hosts = HostRegistry.load()
    det = detect(str(proteome))
    det.filename = f"{BUG}_prepared.faa ({CFG['organism']} {CFG['strain']}, {CFG['n_proteins']} protein TAM PROTEOM)"

    labels = {r.id: (PROTECTIVE.get(r.id, "arka plan")) for r in ingest.load_proteins(str(proteome), det)}
    n_total = len(labels)

    # --- koruyucu antijenlerin VFDB durumu (bağımsız discovery çağrısı) ---
    resolved = cfg.resolve("bacteria")
    prot_records = [p for p in ingest.load_proteins(str(proteome), det) if p.id in protective_accs]
    discovery.run(prot_records, resolved["discovery_vfdb"], profile="bacteria")
    vfdb_status = {p.id: bool(p.annotations.get("vf_hit")) for p in prot_records}

    print(f"{CFG['organism']} TAM PROTEOM: {n_total} protein · Gram={GRAM} · konak=human · taxon={TAXON}")
    t0 = time.time()
    result = None
    for ev in pipeline.run(str(proteome), det, cfg, "bacteria",
                           host_names=["human"], has_gpu=False,
                           outdir=str(OUTDIR.parent), host_registry=hosts,
                           organism_taxon=TAXON, gram=GRAM):
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

    peptides = result["peptides"]
    meta = result["meta"]

    cand_proteins: dict[str, float] = {}
    for p in peptides:
        cand_proteins[p.parent] = max(cand_proteins.get(p.parent, 0.0), p.candidacy)
    cand_accs = set(cand_proteins)

    prot_in_cand = protective_accs & cand_accs
    bg_accs = set(labels) - protective_accs
    bg_in_cand = bg_accs & cand_accs

    p_rate = len(prot_in_cand) / max(1, len(protective_accs))
    b_rate = len(bg_in_cand) / max(1, len(bg_accs))
    fold = round(p_rate / b_rate, 2) if b_rate > 0 else None

    ranked_accs = sorted(cand_proteins, key=lambda a: -cand_proteins[a])
    rank_of = {a: i + 1 for i, a in enumerate(ranked_accs)}

    prot_rows = []
    for acc, lbl in PROTECTIVE.items():
        produced = acc in cand_accs
        vf = vfdb_status.get(acc)
        if produced and vf is False:
            interp = "VFDB-negatif ama aday üretti → soft-filter"
        elif produced and vf is True:
            interp = "VFDB-pozitif, aday üretti"
        elif not produced:
            interp = "aday peptit üretmedi (funnel ya da funnel-sonrası epitop/alerjen/toksisite elemesi)"
        else:
            interp = None
        prot_rows.append({
            "accession": acc, "antigen": lbl,
            "in_testset": acc in labels,
            "vfdb_hit": vf,
            "produced_candidate": produced,
            "best_candidacy": round(cand_proteins.get(acc, 0.0), 4) if produced else None,
            "protein_rank": rank_of.get(acc),
            "n_candidate_proteins": len(cand_accs),
            "interpretation": interp,
        })

    # Fisher + permütasyon-rank + AUC(caveat)
    fish = bstats.fisher_enrichment(len(prot_in_cand), len(protective_accs),
                                    len(bg_in_cand), len(bg_accs))
    perm = None
    if prot_in_cand and len(cand_proteins) >= 2:
        scores = [cand_proteins[a] for a in ranked_accs]
        idxs = [ranked_accs.index(a) for a in prot_in_cand]
        perm = bstats.permutation_rank_test(scores, idxs, n_perm=10000,
                                            statistic="best_rank").as_dict()

    im = meta.get("iedb_match", {})
    summary = {
        "bug": BUG, "organism": CFG["organism"], "strain": CFG["strain"],
        "input": det.filename, "taxon": TAXON, "gram": GRAM,
        "full_proteome": True,
        "runtime_s": round(time.time() - t0, 1),
        "n_testset": n_total, "n_protective": len(PROTECTIVE),
        "n_background": len(bg_accs),
        "n_input_proteins": meta.get("n_input"),
        "n_discovery": meta.get("n_discovery"),
        "n_funnel": meta.get("n_funnel"),
        "n_candidate_proteins": len(cand_accs),
        "recall_at_antigen": f"{len(prot_in_cand)}/{len(protective_accs)}",
        "protective_producing_candidates": sorted(prot_in_cand),
        "background_producing_candidates": len(bg_in_cand),
        "fold_enrichment": fold,
        "fisher_p": fish.p_value,
        "permutation_rank": perm,
        "protective_detail": prot_rows,
        "iedb_source": im.get("source"),
        "iedb_matched_candidates": im.get("n_matched"),
        "report_paths": {k: str(v) for k, v in result["paths"].items()},
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (RESULTS / "vaxforge_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    import shutil
    for k, v in result["paths"].items():
        try:
            shutil.copy(v, OUTDIR / Path(v).name)
        except Exception:
            pass

    print("\n" + "=" * 70)
    print(f"VaxForge TAM PROTEOM ÖZETİ — {CFG['organism']} ({CFG['strain']})")
    print("=" * 70)
    print(f"Proteom: {n_total} · keşif: {summary['n_discovery']} · huni: {summary['n_funnel']} · "
          f"aday-protein: {summary['n_candidate_proteins']} · süre: {summary['runtime_s']}s")
    print(f"\nRecall@antijen: {summary['recall_at_antigen']}  ·  fold: {fold}  ·  Fisher p: {fish.p_value:.3g}")
    for r in prot_rows:
        mark = "✅" if r["produced_candidate"] else "❌"
        vf = {True: "VFDB+", False: "VFDB−", None: "VFDB?"}[r["vfdb_hit"]]
        rnk = (f"sıra #{r['protein_rank']}/{r['n_candidate_proteins']} skor {r['best_candidacy']}"
               if r["produced_candidate"] else "aday peptit üretmedi")
        print(f"  {mark} [{vf}] {r['antigen']:52s} {rnk}")
    print(f"\nÇıktılar: {RESULTS}/vaxforge_summary.json  +  {OUTDIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
