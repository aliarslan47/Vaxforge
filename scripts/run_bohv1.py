"""BoHV-1 (bovine alphaherpesvirus 1, JX898220) — sığır konağı re-run.

16 Kasım 2026 koşusu (outputs/run_2026-11-16T21-59-16) B-hücre adaylarını
Kolaskar-Tongaonkar proxy'siyle üretmiş; IEDB / ImmForge / çözünürlük /
düzensizlik / ikincil-yapı adımları boş çıkmıştı. O tarihten beri:
  - araçlar gerçekten kuruldu ve çalıştırılarak doğrulandı (21/23),
  - ImmForge'un sessizce devre dışı kalmasına yol açan yol çözümü düzeltildi,
  - örtüşen aday kümeleme adımı eklendi.
Bu script o koşunun aynı girdiyle, aynı konakla tekrarıdır.

profil=virus, host=bovine, taxon=10320 (GenBank dosyasından doğrulandı).
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from vaxforge import pipeline
from vaxforge.config_loader import ThresholdConfig
from vaxforge.detect import detect
from vaxforge.hosts import HostRegistry

INPUT = Path(__file__).resolve().parent.parent / "data" / "uploads" / "bohv1.gb"
if not INPUT.is_file():
    raise SystemExit(f"girdi yok: {INPUT}")

cfg = ThresholdConfig.load()
problems = cfg.validate()
if problems:
    raise SystemExit("config geçersiz: " + "; ".join(problems))

reg = HostRegistry.load()
det = detect(str(INPUT))
print(f"[detect] molecule={det.molecule} num_seqs={det.num_seqs}", flush=True)

t0 = time.time()
result = None
for ev in pipeline.run(
    str(INPUT), det, cfg, "virus",
    host_names=["bovine"],
    organism_taxon="NCBITaxon:10320",
    outdir="outputs", host_registry=reg,
):
    if ev["phase"] == "__result__":
        result = ev["data"]
    elif ev["phase"] == "__error__":
        print(f"[HATA] {ev['msg']}", flush=True)
    else:
        print(f"[{time.time() - t0:7.1f}s] {ev['phase']:14} {ev['status']:8} {ev['msg']}",
              flush=True)

if not result:
    raise SystemExit("koşu sonuç üretmedi")

meta = result.get("meta", {})
print("\n=== BİTTİ (%.1f dk) ===" % ((time.time() - t0) / 60), flush=True)
print("aday peptit:", len(result.get("peptides", [])), flush=True)
for k in ("n_raw", "n_discovery", "n_funnel", "n_epitope", "n_scored", "n_survivors"):
    print(f"  {k} = {meta.get(k)}", flush=True)
print("  kümeleme =", meta.get("clustering"), flush=True)
print("  araç özeti =", meta.get("tool_status_summary"), flush=True)
print("çıktılar:", flush=True)
for k, v in result.get("paths", {}).items():
    print(f"  {k}: {v}", flush=True)
