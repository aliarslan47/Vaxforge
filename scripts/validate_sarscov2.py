"""VaxForge validasyon koşusu — SARS-CoV-2 (yayın için, tekrarlanabilir).

Referans proteom (UniProt UP000464024, 17 protein) VaxForge pipeline'ından
(virus profili, insan konak) geçirilir; kalan aday peptitler IEDB'de deneysel
doğrulanmış SARS-CoV-2 epitoplarıyla eşleştirilir ve recall/precision benchmark'ı
hesaplanır. Tüm çıktılar outputs/validation_sarscov2/ altına kaydedilir.

Çalıştırma:
    . .venv/bin/activate
    python3 scripts/validate_sarscov2.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import pipeline  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402

PROTEOME = ROOT / "data" / "validation" / "sarscov2_proteome.faa"
TAXON = "NCBITaxon:2697049"
OUTDIR = ROOT / "outputs" / "validation_sarscov2"

# Bilinen/immünodominant SARS-CoV-2 epitopları (literatür landmark'ları) —
# pipeline bunları yakalıyor mu diye ayrıca kontrol edilir.
LANDMARKS = {
    "YLQPRTFLL": "Spike (HLA-A*02:01, immünodominant CD8)",
    "KCYGVSPTKL": "Spike (CD8)",
    "SPRWYFYYL": "Nucleoprotein (HLA-B*07:02, immünodominant CD8)",
    "LALLLLDRLNQL": "Nucleoprotein (CD8 bölge)",
    "NYNYLYRLF": "Spike (CD8)",
}


def main() -> int:
    if not PROTEOME.exists():
        print(f"HATA: proteom yok: {PROTEOME}")
        return 1
    cfg = ThresholdConfig.load()
    hosts = HostRegistry.load()
    det = detect(str(PROTEOME))
    det.filename = "sarscov2_proteome.faa (UniProt UP000464024)"

    print(f"Girdi: {det.num_seqs} dizi · profil=virus · konak=human · taxon={TAXON}")
    t0 = time.time()
    result = None
    for ev in pipeline.run(str(PROTEOME), det, cfg, "virus",
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
    im = meta.get("iedb_match", {})
    bm = (im or {}).get("benchmark") or {}

    # --- landmark kontrolü: bilinen epitoplar aday/örtüşme olarak yakalandı mı ---
    pred = [p.seq.upper() for p in peptides]
    landmark_rows = []
    for ep, desc in LANDMARKS.items():
        hit = any(ep in s or s in ep or _shares(ep, s, 8) for s in pred)
        landmark_rows.append({"epitope": ep, "desc": desc, "recovered": hit})

    matched = [{
        "rank": i + 1, "peptide": p.seq, "kind": p.kind,
        "candidacy": p.candidacy,
        "match_type": (p.metrics.get("iedb") or {}).get("match_type"),
        "iedb_epitope": (p.metrics.get("iedb") or {}).get("epitope_seq"),
        "antigen": ((p.metrics.get("iedb") or {}).get("antigens") or ["—"])[0],
        "pmids": (p.metrics.get("iedb") or {}).get("pmids") or [],
    } for i, p in enumerate(peptides) if (p.metrics.get("iedb") or {}).get("matched") is True]

    summary = {
        "input": det.filename, "taxon": TAXON,
        "runtime_s": round(time.time() - t0, 1),
        "n_input_proteins": meta.get("n_input"),
        "n_funnel": meta.get("n_funnel"),
        "n_candidates": len(peptides),
        "iedb_source": im.get("source"),
        "benchmark": bm,
        "n_candidates_matching_known": len(matched),
        "landmark_recovery": landmark_rows,
        "matched_candidates": matched,
        "report_paths": {k: str(v) for k, v in result["paths"].items()},
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n" + "=" * 70)
    print("VALIDASYON ÖZETİ — SARS-CoV-2")
    print("=" * 70)
    print(f"Girdi protein: {summary['n_input_proteins']}  ·  huni sonrası: "
          f"{summary['n_funnel']}  ·  aday: {summary['n_candidates']}  ·  "
          f"süre: {summary['runtime_s']}s")
    if bm.get("recall") is not None:
        print(f"\nBenchmark (IEDB ground truth):")
        print(f"  Bilinen benzersiz epitop : {bm['n_known']}")
        print(f"  Yakalanan (recall)       : {bm['n_known_hit']} ({round(bm['recall']*100,2)}%)")
        print(f"  Aday sayısı              : {bm['n_pred']}")
        print(f"  Bilinene eşleşen aday    : {bm['n_pred_matched']} "
              f"(precision-benzeri {round((bm['precision_like'] or 0)*100,1)}%)")
    print(f"\nLandmark (immünodominant) epitop kurtarma:")
    for r in landmark_rows:
        print(f"  {'✅' if r['recovered'] else '❌'} {r['epitope']:14s} {r['desc']}")
    print(f"\nBilinen epitopla eşleşen ilk adaylar:")
    for m in matched[:12]:
        pm = ", ".join(m["pmids"][:2]) or "—"
        print(f"  #{m['rank']:<3d} {m['peptide']:16s} {m['kind']:6s} "
              f"skor={m['candidacy']:.3f}  {m['match_type']:18s} PMID {pm}")
    print(f"\nÇıktılar: {OUTDIR}/  (validation_summary.json + rapor paketi)")
    print(f"Rapor PDF: {summary['report_paths'].get('pdf')}")
    return 0


def _shares(a: str, b: str, k: int) -> bool:
    return any(a[i:i + k] in b for i in range(len(a) - k + 1))


if __name__ == "__main__":
    raise SystemExit(main())
