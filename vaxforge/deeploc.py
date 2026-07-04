"""GERÇEK hücre-altı lokalizasyon — DeepLoc-2.1 (DTU, ESM tabanlı, yerel).

deeploc2 CLI kuruluysa gerçek lokalizasyon + sinyal + membran tipi tahmini yapar.
Yoksa available()=False -> funnel.py Kyte-Doolittle/sinyal heuristiğine düşer.

NOT: DeepLoc ökaryot-eğilimlidir; bakteride 'Extracellular'/'Cell membrane'
sinyalleri yine de yüzey-erişilebilirlik için anlamlıdır. İlk çalıştırmada ESM1b
modelini indirir (~2.5GB, önbelleğe). CPU'da protein başına birkaç saniye.
"""

from __future__ import annotations

import csv
import functools
import glob
import shutil
import subprocess
import tempfile
from pathlib import Path

# DeepLoc lokalizasyonu -> config allowed sözlüğü
_MAP = {
    "Extracellular": "extracellular",
    "Cell membrane": "outer_membrane",   # yüzey-erişilebilir
    "Cell wall": "cell_wall",
    "Periplasmic space": "periplasm",
}


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return shutil.which("deeploc2") is not None


def predict(proteins: list[tuple[str, str]], model: str = "Fast") -> dict[str, dict]:
    """[(id, seq)] -> {id: {localization(vocab), raw, signal, membrane}}.

    Girdi id'leri s0,s1.. temiz olarak verilir (özel karakter sorununu önler).
    """
    if not available() or not proteins:
        return {}
    idmap = {f"s{i}": pid for i, (pid, _) in enumerate(proteins)}
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        for i, (_, seq) in enumerate(proteins):
            fh.write(f">s{i}\n{seq}\n")
        fpath = fh.name
    outdir = tempfile.mkdtemp(prefix="deeploc_")
    try:
        subprocess.run(["deeploc2", "-f", fpath, "-o", outdir, "-m", model],
                       capture_output=True, timeout=1800, text=True)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        shutil.rmtree(outdir, ignore_errors=True)
        return {}
    Path(fpath).unlink(missing_ok=True)
    out: dict[str, dict] = {}
    csvs = sorted(glob.glob(f"{outdir}/*.csv"))
    if csvs:
        with open(csvs[-1]) as fh:
            for row in csv.DictReader(fh):
                sid = row.get("Protein_ID", "")
                pid = idmap.get(sid, sid)
                raw = row.get("Localizations", "")
                locs = [x.strip() for x in raw.split("|")]
                vocab = "cytoplasm"
                for loc in locs:
                    if loc in _MAP:
                        vocab = _MAP[loc]
                        break
                out[pid] = {"localization": vocab, "raw": raw,
                            "signal": row.get("Signals", ""),
                            "membrane": row.get("Membrane types", "")}
    shutil.rmtree(outdir, ignore_errors=True)
    return out
