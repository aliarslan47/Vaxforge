"""DeepLoc-2.0 — ökaryot hücre-altı lokalizasyon (venv CLI, Fast/ESM1b).

PSORTb prokaryota-özeldir; parazit (ökaryot) ve virüs proteinleri için DeepLoc
kullanılır. DeepLoc 10 ökaryot bölmesini (sitoplazma / çekirdek / hücre-dışı /
hücre membranı / mitokondri / plastid / ER / lizozom-vakuol / Golgi / peroksizom)
ve sinyal tiplerini tahmin eder. RV açısından hücre-dışı (salgı) + hücre membranı
= erişilebilir yüzey antijeni adayları.

Kurulum: `deeploc2` venv CLI + ESM1b ağırlığı (~7GB, torch hub önbelleği). İkisi de
yoksa available()=False → funnel heuristiğe düşer (dürüst 'method' etiketi). Model
İNDİRME tetiklememek için available() önbelleği de kontrol eder (uydurma/askıda kalma yok).
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
_DEEPLOC = _REPO / ".venv" / "bin" / "deeploc2"
_ESM1B = Path.home() / ".cache" / "torch" / "hub" / "checkpoints" / "esm1b_t33_650M_UR50S.pt"

# DeepLoc bölmesi -> config'te kullanılan normalize lokalizasyon anahtarı.
# Parazit allowed = [extracellular, plasma_membrane, secreted, gpi_anchored]
# Virüs allowed   = [surface, envelope, membrane]  (DeepLoc virionlar için güvenilmez → yumuşak)
_LOC_MAP = {
    "extracellular": "extracellular",
    "cellmembrane": "plasma_membrane",
    "cytoplasm": "cytoplasm",
    "nucleus": "nucleus",
    "mitochondrion": "mitochondrion",
    "plastid": "plastid",
    "endoplasmicreticulum": "endoplasmic_reticulum",
    "lysosomevacuole": "lysosome_vacuole",
    "golgiapparatus": "golgi",
    "peroxisome": "peroxisome",
}


def available() -> bool:
    """deeploc2 CLI var VE ESM1b ağırlığı önbellekte mi? (indirme tetiklemeden)."""
    return _DEEPLOC.exists() and _ESM1B.exists()


def _norm(loc: str) -> str:
    key = loc.strip().replace(" ", "").replace("/", "").replace("-", "").lower()
    return _LOC_MAP.get(key, "unknown")


def _parse_csv(csv_path: str, id_map: dict[str, str]) -> dict[str, dict]:
    """DeepLoc results_*.csv → {gerçek_id: {localization, raw, signals, score}}.

    Sütunlar: Protein_ID,Localizations,Signals,Membrane types,<10 olasılık>...
    'Localizations' pipe-ayrık; primary = ilk. Olasılık sütunundan güven okunur.
    """
    import csv
    out: dict[str, dict] = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sid = (row.get("Protein_ID") or "").strip()
            if sid not in id_map:
                continue
            locs = (row.get("Localizations") or "").strip()
            primary = locs.split("|")[0].strip() if locs else ""
            mapped = _norm(primary)
            # primary bölmenin olasılığı (varsa) — güven skoru
            try:
                score = float(row.get(primary, "") or 0.0)
            except (ValueError, TypeError):
                score = 0.0
            out[id_map[sid]] = {
                "localization": mapped,
                "raw": locs or "—",
                "signals": (row.get("Signals") or "").strip(),
                "score": round(score, 2),
            }
    return out


def predict(seq_pairs: list[tuple[str, str]], model: str = "Fast",
            timeout: int = 3600) -> dict[str, dict]:
    """[(id, seq)] -> {id: {localization, raw, signals, score}}. Hata/eksik -> {}.

    Sentetik boşluksuz ID (s0, s1…) kullanılır (DeepLoc başlığı boşlukta kırpar),
    sonra gerçek ID'ye çevrilir. Hata halinde {} → çağıran heuristiğe düşer (dürüst).
    """
    if not available() or not seq_pairs:
        return {}
    workdir = Path(tempfile.mkdtemp(prefix="deeploc_"))
    infile = workdir / "input.faa"
    outdir = workdir / "out"
    outdir.mkdir(exist_ok=True)
    id_map = {f"s{i}": pid for i, (pid, _) in enumerate(seq_pairs)}
    infile.write_text(
        "".join(f">s{i}\n{seq}\n" for i, (_, seq) in enumerate(seq_pairs)),
        encoding="utf-8")
    try:
        cmd = [str(_DEEPLOC), "-f", str(infile), "-o", str(outdir),
               "-m", model, "-d", "cpu"]
        env = {**os.environ, "OMP_NUM_THREADS": "4"}
        subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                       cwd=str(workdir), env=env)
    except Exception:
        shutil.rmtree(workdir, ignore_errors=True)
        return {}
    csvs = glob.glob(str(outdir / "results_*.csv"))
    if not csvs:
        shutil.rmtree(workdir, ignore_errors=True)
        return {}
    try:
        result = _parse_csv(csvs[0], id_map)
    except Exception:
        result = {}
    shutil.rmtree(workdir, ignore_errors=True)
    return result
