"""ebola.gb (Zaire ebolavirus) — virüs profili, insan konağı. Küçük genom → hızlı.
profil=virus, host=human, gram=None, taxon=NCBITaxon:186538 (hızlı IEDB).
"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from vaxforge import pipeline
from vaxforge.config_loader import ThresholdConfig
from vaxforge.detect import detect
from vaxforge.hosts import HostRegistry

SRC = "/mnt/c/Users/msi-nb/Desktop/ebola.gb"
work = Path("/home/msi-nb/vaxforge/data/uploads/ebola.gb")
work.parent.mkdir(parents=True, exist_ok=True)
work.write_bytes(Path(SRC).read_bytes())

cfg = ThresholdConfig.load()
reg = HostRegistry.load()
det = detect(str(work))
print(f"[detect] molecule={det.molecule} num_seqs={det.num_seqs}", flush=True)

t0 = time.time()
result = None
for ev in pipeline.run(str(work), det, cfg, "virus",
                       host_names=["human"], gram=None,
                       organism_taxon="NCBITaxon:186538",
                       outdir="outputs", host_registry=reg):
    if ev["phase"] == "__result__":
        result = ev["data"]
    else:
        print(f"[{time.time()-t0:6.1f}s] {ev['phase']:12} {ev['status']:8} {ev['msg']}", flush=True)

if result:
    paths = result.get("paths", {})
    print("\n=== BİTTİ ===", flush=True)
    print("aday peptit sayısı:", len(result.get("peptides", [])), flush=True)
    for k, v in paths.items():
        print(f"  {k}: {v}", flush=True)
