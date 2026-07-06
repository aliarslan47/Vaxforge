"""VaxForge validasyon koşusu — Neisseria meningitidis B (Bexsero/4CMenB).

ANTİJEN-SEVİYESİ POZİTİF KONTROL: reverse vaccinology'nin doğuş vakası. Bexsero
aşısının bilinen KORUYUCU antijenleri (fHbp, NadA, NHBA + porA) MC58 proteomundan
alınır, rastgele ARKA PLAN proteinleriyle karıştırılır ve bakteri pipeline'ından
geçirilir. Soru: pipeline (keşif→huni→skor) bilinen koruyucu antijenleri arka
plandan AYIRT edip yüksek sıralıyor mu? (Protegen fold-enrichment mantığının
tractable, downsample edilmiş versiyonu.)

Çalıştırma:
    . .venv/bin/activate
    python3 scripts/validate_menb.py
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO  # noqa: E402

from vaxforge import pipeline  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402

PROTEOME = ROOT / "data" / "validation" / "menb_mc58_proteome.faa"
TESTSET = ROOT / "data" / "validation" / "menb_testset.faa"
TAXON = "NCBITaxon:487"          # Neisseria meningitidis (tür seviyesi, IEDB)
OUTDIR = ROOT / "outputs" / "validation_menb"
N_BACKGROUND = 60
SEED = 42

# Bexsero'nun bilinen koruyucu antijenleri (UniProt aksesyon -> etiket).
PROTECTIVE = {
    "Q9JXV4": "fHbp (Factor H binding protein, NMB1870)",
    "Q9JXK7": "NadA (Neisseria adhesin A, NMB1994)",
    "Q7DD37": "NHBA (Neisserial Heparin Binding Antigen, NMB2132)",
    "P0DH58": "PorA (Major outer membrane protein, OMV bileşeni)",
}


def _acc(rec) -> str:
    return rec.id.split("|")[1] if "|" in rec.id else rec.id


def build_testset() -> dict[str, str]:
    """Test setini (koruyucu + rastgele arka plan) deterministik kur. -> {acc: label}."""
    recs = {_acc(r): r for r in SeqIO.parse(str(PROTEOME), "fasta")}
    prot = {a: r for a, r in recs.items() if a in PROTECTIVE}
    pool = [a for a in recs if a not in PROTECTIVE]
    random.Random(SEED).shuffle(pool)
    bg = pool[:N_BACKGROUND]
    chosen = list(prot) + bg
    with TESTSET.open("w") as fh:
        for a in chosen:
            r = recs[a]
            fh.write(f">{a} {r.description.split(None, 1)[-1][:60]}\n{r.seq}\n")
    labels = {a: PROTECTIVE.get(a, "arka plan") for a in chosen}
    return labels


def main() -> int:
    if not PROTEOME.exists():
        print(f"HATA: proteom yok: {PROTEOME}")
        return 1
    labels = build_testset()
    protective_accs = set(PROTECTIVE)
    n_total = len(labels)
    cfg = ThresholdConfig.load()
    hosts = HostRegistry.load()
    det = detect(str(TESTSET))
    det.filename = f"menb_testset.faa ({len(PROTECTIVE)} koruyucu + {N_BACKGROUND} arka plan)"

    print(f"Test seti: {n_total} protein ({len(PROTECTIVE)} koruyucu antijen + "
          f"{N_BACKGROUND} arka plan) · profil=bacteria · konak=human · taxon={TAXON}")
    t0 = time.time()
    result = None
    n_discovery = n_funnel = None
    survivors_after = {}  # phase -> set(acc) — hangi proteinler o adımı geçti
    for ev in pipeline.run(str(TESTSET), det, cfg, "bacteria",
                           host_names=["human"], has_gpu=False,
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

    peptides = result["peptides"]
    meta = result["meta"]

    # --- aday peptitleri kaynak proteine indir ---
    cand_proteins = {}  # acc -> en iyi aday skoru
    for p in peptides:
        acc = p.parent
        cand_proteins[acc] = max(cand_proteins.get(acc, 0.0), p.candidacy)
    cand_accs = set(cand_proteins)

    prot_in_cand = protective_accs & cand_accs
    bg_accs = set(labels) - protective_accs
    bg_in_cand = bg_accs & cand_accs

    # fold-enrichment: koruyucuların aday-üreten oranı / arka planın oranı
    p_rate = len(prot_in_cand) / max(1, len(protective_accs))
    b_rate = len(bg_in_cand) / max(1, len(bg_accs))
    fold = round(p_rate / b_rate, 2) if b_rate > 0 else None

    # koruyucu antijenlerin durum tablosu
    prot_rows = []
    # sıralama: her proteinin en iyi aday sırası
    ranked_accs = sorted(cand_proteins, key=lambda a: -cand_proteins[a])
    rank_of = {a: i + 1 for i, a in enumerate(ranked_accs)}
    for acc, lbl in PROTECTIVE.items():
        in_input = acc in labels
        produced = acc in cand_accs
        prot_rows.append({
            "accession": acc, "antigen": lbl,
            "in_testset": in_input,
            "produced_candidate": produced,
            "best_candidacy": round(cand_proteins.get(acc, 0.0), 4) if produced else None,
            "protein_rank": rank_of.get(acc),
            "n_candidate_proteins": len(cand_accs),
        })

    im = meta.get("iedb_match", {})
    summary = {
        "input": det.filename, "taxon": TAXON, "seed": SEED,
        "runtime_s": round(time.time() - t0, 1),
        "n_testset": n_total, "n_protective": len(PROTECTIVE),
        "n_background": N_BACKGROUND,
        "n_input_proteins": meta.get("n_input"),
        "n_discovery": meta.get("n_discovery"),
        "n_funnel": meta.get("n_funnel"),
        "n_candidate_proteins": len(cand_accs),
        "protective_producing_candidates": sorted(prot_in_cand),
        "background_producing_candidates": len(bg_in_cand),
        "fold_enrichment": fold,
        "protective_detail": prot_rows,
        "iedb_source": im.get("source"),
        "iedb_matched_candidates": im.get("n_matched"),
        "report_paths": {k: str(v) for k, v in result["paths"].items()},
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    # rapor paketini yayın klasörüne kopyala
    import shutil
    for k, v in result["paths"].items():
        try:
            shutil.copy(v, OUTDIR / Path(v).name)
        except Exception:
            pass

    print("\n" + "=" * 70)
    print("VALIDASYON ÖZETİ — Neisseria meningitidis B (Bexsero)")
    print("=" * 70)
    print(f"Test seti: {n_total} protein  ·  keşif sonrası: {summary['n_discovery']}  ·  "
          f"huni sonrası: {summary['n_funnel']}  ·  aday-üreten protein: "
          f"{summary['n_candidate_proteins']}  ·  süre: {summary['runtime_s']}s")
    print(f"\nBilinen koruyucu antijenlerin izi (pipeline onları buldu mu?):")
    for r in prot_rows:
        mark = "✅" if r["produced_candidate"] else "❌"
        rnk = (f"protein sırası #{r['protein_rank']}/{r['n_candidate_proteins']}, "
               f"en iyi skor {r['best_candidacy']}") if r["produced_candidate"] else "aday üretmedi (elendi)"
        print(f"  {mark} {r['antigen']:52s} {rnk}")
    print(f"\nFold-enrichment (koruyucu aday-oranı / arka plan aday-oranı): "
          f"{summary['fold_enrichment']}")
    print(f"  koruyucu aday üreten : {len(prot_in_cand)}/{len(protective_accs)}")
    print(f"  arka plan aday üreten: {len(bg_in_cand)}/{len(bg_accs)}")
    print(f"IEDB (taxon {TAXON}) eşleşen aday: {summary['iedb_matched_candidates']}")
    print(f"\nÇıktılar: {OUTDIR}/  (validation_summary.json + rapor paketi)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
