"""br1.gb (Brucella abortus) — sığır konağı, arayüz ayarlarının birebir eşi.
profil=bacteria, host=bovine, gram=negative.
"""
import sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from vaxforge import pipeline
from vaxforge.config_loader import ThresholdConfig
from vaxforge.detect import detect
from vaxforge.hosts import HostRegistry

SRC = "/mnt/c/Users/msi-nb/Desktop/br1.gb"
work = Path("/home/msi-nb/vaxforge/data/uploads/br1.gb")
work.parent.mkdir(parents=True, exist_ok=True)
work.write_bytes(Path(SRC).read_bytes())

cfg = ThresholdConfig.load()
reg = HostRegistry.load()
det = detect(str(work))
print(f"[detect] molecule={det.molecule} num_seqs={det.num_seqs}", flush=True)

t0 = time.time()
result = None
for ev in pipeline.run(str(work), det, cfg, "bacteria",
                       host_names=["bovine"], gram="negative",
                       outdir="outputs", host_registry=reg):
    if ev["phase"] == "__result__":
        result = ev["data"]
    else:
        dt = time.time() - t0
        print(f"[{dt:6.1f}s] {ev['phase']:12} {ev['status']:8} {ev['msg']}", flush=True)

if result:
    meta = result.get("meta", {})
    paths = result.get("paths", {})
    print("\n=== BİTTİ ===", flush=True)
    print("aday peptit sayısı:", len(result.get("peptides", [])), flush=True)
    print("çıktı dosyaları:", flush=True)
    for k, v in paths.items():
        print(f"  {k}: {v}", flush=True)
