"""GERÇEK lineer B-hücre epitop tahmini — BepiPred-3.0 (DTU, ESM-2, yerel).

BepiPred-3.0 (Clifford ve ark. 2022, Protein Science): ESM-2 (650M) protein dil
modeli kodlaması + 5-fold FFNN topluluğu ile protein başına per-residue B-hücre
epitop olasılığı. BepiPred-1.0'a (2006, HMM + Parker) göre belirgin daha yüksek
doğruluk (AUC).

Repo: tools/bepipred-3.0 (bp3 paketi + eğitilmiş Fold1-5 modelleri gömülü).
ESM-2 modeli ilk koşuda torch hub önbelleğine iner (~2.6GB). CPU'da yavaş ve
bellek O(L²) olduğundan:
  - TÜM (funnel'dan geçen) proteinler için TEK batch koşulur,
  - çok uzun proteinler (> MAX_LEN) batch'ten çıkarılır → çağıran BepiPred-1.0'a
    ya da klasik yönteme düşer (batch çökmesin diye).

Kuruluysa {id: [per-residue skor]} döndürür; yoksa available()=False →
epitope.py BepiPred-1.0'a, o da yoksa Kolaskar-Tongaonkar+Parker'a düşer.
"""

from __future__ import annotations

import csv
import functools
import subprocess
import sys
import tempfile
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent / "tools" / "bepipred-3.0"
_CLI = _ROOT / "bepipred3_CLI.py"

# BepiPred-3.0 web server varsayılan değişken eşiği (ortalama pozitif olasılık).
DEFAULT_THRESHOLD = 0.1512
# ESM-2 (650M) CPU'da O(L²) bellek; çok uzun proteinleri batch'ten çıkar.
MAX_LEN = 1024


@functools.lru_cache(maxsize=1)
def available() -> bool:
    if not _CLI.exists():
        return False
    import importlib.util
    # bp3 çekirdeğinin zincir bağımlılığı (plotly) kurulu mu?
    return importlib.util.find_spec("plotly") is not None


def predict_batch(seq_pairs: list[tuple[str, str]]) -> dict[str, list[float]]:
    """[(id, seq)] -> {id: [per-residue BepiPred-3.0 skoru]}.

    Yalnız uzunluğu MAX_LEN'i aşmayan proteinler koşulur; diğerleri sözlükte
    yer almaz (çağıran yedek yönteme düşürür). Herhangi bir hata → {} (tam yedek).
    """
    if not available() or not seq_pairs:
        return {}

    id_map: dict[str, str] = {}   # bp3protN -> gerçek id
    with tempfile.TemporaryDirectory() as tmpname:
        tmp = Path(tmpname)
        fasta = tmp / "in.fasta"
        n_written = 0
        with fasta.open("w") as fh:
            for n, (pid, seq) in enumerate(seq_pairs):
                if not seq or len(seq) > MAX_LEN:
                    continue
                tag = f"bp3prot{n}"
                id_map[tag] = pid
                fh.write(f">{tag}\n{seq}\n")
                n_written += 1
        if n_written == 0:
            return {}

        outdir = tmp / "out"
        outdir.mkdir()
        cmd = [sys.executable, str(_CLI), "-i", str(fasta), "-o", str(outdir),
               "-pred", "vt_pred", "-esm_dir", str(tmp / "esm")]
        try:
            subprocess.run(cmd, cwd=str(_ROOT), capture_output=True, text=True,
                           timeout=14400, check=True)
        except Exception:
            return {}

        csv_path = outdir / "raw_output.csv"
        if not csv_path.exists():
            return {}
        scores: dict[str, list[float]] = {}
        with csv_path.open() as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                acc = (row.get("Accession") or "").strip()
                pid = id_map.get(acc, acc)
                try:
                    scores.setdefault(pid, []).append(float(row["BepiPred-3.0 score"]))
                except (KeyError, ValueError):
                    continue
    return scores
