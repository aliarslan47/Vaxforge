"""Örnek veri üzerinde ESKİ vs YENİ toksisite kapısını aynı peptit setinde karşılaştır.

Tek bir pipeline koşusu yapılır (örnek demo proteom); üretilen aday peptitler üzerinde
iki kapı yan yana değerlendirilir:

  ESKİ : ToxinPred2 Model 2 (hibrit), eşik 0.6, SERT filtre (eleme)
  YENİ : ToxinPred3 ExtraTrees, config eşikleri, YUMUŞAK (işaret + skor cezası),
         eleme yalnız hard_discard_score üstünde

Çıktı: results/demo_toxicity_gate.json + ekrana özet tablo.

Koşum: .venv/bin/python3 benchmark/demo_toxicity_gate.py [--input FASTA]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import pipeline, toxinpred, toxinpred3  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402

DEMO = ROOT / "data" / "samples" / "proteome_demo.faa"
OUT = ROOT / "benchmark" / "results" / "demo_toxicity_gate.json"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DEMO))
    ap.add_argument("--run-json", default=None,
                    help="var olan bir run.json üzerinden koş (pipeline tekrar koşmaz)")
    args = ap.parse_args()

    if args.run_json:
        rj = json.loads(Path(args.run_json).read_text(encoding="utf-8"))
        print(f"var olan koşu okundu: {args.run_json}", flush=True)
        compare(rj, args.run_json)
        return

    print(f"girdi: {args.input}\npipeline koşuyor…", flush=True)
    det = detect(args.input)
    det.filename = Path(args.input).name
    cfg0 = ThresholdConfig.load()
    outdir = ROOT / "outputs" / "demo_toxicity_gate"
    result = None
    for ev in pipeline.run(args.input, det, cfg0, "bacteria",
                           host_names=["human"], overrides={}, has_gpu=False,
                           outdir=str(outdir), host_registry=HostRegistry.load(),
                           organism_taxon=None, gram="negative", lang="tr"):
        if ev["phase"] == "__error__":
            print("pipeline hatası:", ev.get("msg")); return
        if ev["phase"] == "__result__":
            result = ev["data"]
        elif ev.get("status") == "done" and ev["phase"].startswith("survival"):
            print("  " + ev.get("msg", ""), flush=True)
    print("pipeline bitti.", flush=True)

    # run.json tek doğruluk kaynağı: __result__ olayının şeması rapor katmanı
    # için sadeleştirilmiş olabiliyor (candidates/survival taşımıyor).
    outs = sorted(outdir.glob("run_*"), key=lambda p: p.stat().st_mtime)
    rj = json.loads((outs[-1] / "run.json").read_text(encoding="utf-8")) if outs else (result or {})
    compare(rj, args.input)


def compare(rj: dict, source: str) -> None:
    peps = sorted({c["peptide"] for c in rj.get("candidates", []) if c.get("peptide")})
    surv = rj.get("survival") or {}
    print(f"aday peptit: {len(peps)} · sağ kalım özeti: {surv}", flush=True)
    if not peps:
        print("aday yok — karşılaştırma yapılamadı"); return

    cfg = ThresholdConfig.load()
    tox = cfg.resolve(rj.get("profile", "bacteria"))["toxicity"]
    soft = float(tox.params["max_toxicity_score"].value)
    hard = float(tox.params["hard_discard_score"].value)

    old = toxinpred.predict(peps, threshold=0.6)
    new = toxinpred3.predict(peps, threshold=soft, model=1)

    old_kill = sum(1 for p in peps if old.get(p, {}).get("toxic"))
    new_flag = sum(1 for p in peps if (new.get(p, {}).get("ml_score", 0) or 0) >= soft)
    new_kill = sum(1 for p in peps if (new.get(p, {}).get("ml_score", 0) or 0) >= hard)

    rep = {
        "input": source, "n_candidates": len(peps),
        "note": ("peptit seti = koşunun NİHAİ aday listesi; 'eski kapı' sayısı, "
                 "ToxinPred2'nin sert uygulanması halinde kaybedilecek nihai "
                 "aşı adaylarını gösterir"),
        "survival_summary": surv,
        "old_gate": {"engine": "ToxinPred2 hibrit", "threshold": 0.6, "mode": "sert",
                     "eliminated": old_kill,
                     "eliminated_pct": round(100 * old_kill / len(peps), 1)},
        "new_gate": {"engine": "ToxinPred3 ET", "soft_threshold": soft,
                     "hard_threshold": hard, "mode": "yumuşak + yüksek-güven eleme",
                     "flagged": new_flag,
                     "flagged_pct": round(100 * new_flag / len(peps), 1),
                     "eliminated": new_kill,
                     "eliminated_pct": round(100 * new_kill / len(peps), 1)},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print("\n" + "=" * 62)
    print(f"{len(peps)} aday peptit, aynı set üzerinde iki kapı")
    print("-" * 62)
    print(f"ESKİ  ToxinPred2 hibrit t=0.6 (sert)   → {old_kill:5d} elenir "
          f"(%{rep['old_gate']['eliminated_pct']})")
    print(f"YENİ  ToxinPred3 ET t={soft} (yumuşak)  → {new_flag:5d} işaretlenir "
          f"(%{rep['new_gate']['flagged_pct']}), skor ×ceza")
    print(f"      ToxinPred3 ET t={hard} (eleme)   → {new_kill:5d} elenir "
          f"(%{rep['new_gate']['eliminated_pct']})")
    print("=" * 62)
    print(f"yazıldı: {OUT}")


if __name__ == "__main__":
    main()
