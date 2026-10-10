"""GERÇEK toksisite tahmini — ToxinPred3 (Raghava lab, pip paketi).

NEDEN TOXINPRED2'DEN AYRI BİR MODÜL: ToxinPred2'nin başlığı bile "toxicity of
*proteins*"tir; AAC+RF modeli tam uzunluk proteinler/uzun peptitler üzerinde
eğitilmiştir ve 9-mer epitoplarda dağılım dışı kalır. 2026-10-10 kalibrasyonunda
ToxinPred2 (Model 2 hibrit) insan self 9-mer'lerinin %74'ünü "toksik" işaretledi
ve influenza M1(58-66) GILGFVFTL gibi deneysel doğrulanmış zararsız epitopları da
eledi. ToxinPred3 açıkça PEPTİT düzeyi için yazılmıştır (varsayılan eşik 0.38).

Atıf: Rathore AS, Choudhury S, Arora A, Tijare P, Raghava GPS. ToxinPred 3.0: An
improved method for predicting the toxicity of peptides. Comput Biol Med.
2024;179:108926. doi:10.1016/j.compbiomed.2024.108926

Model dosyası eski sklearn ile pickle'landığı için bir kez dönüştürülür:
`scripts/convert_toxinpred3_model.py` (IFNepitope2 ile aynı köprü mantığı).

Model 1 = AAC+DPC ExtraTrees (yalnız ML). Model 2 = Hybrid (ML + MERCI motif).
ToxinPred2'nin hibritindeki BLAST toksin-homolojisi bileşeni ToxinPred3'te YOK;
bu bizim için önemli, çünkü ters aşıbilimde antijenlerin bir kısmı zaten toksindir
(listeriolizin O bir hemolizin; difteri/tetanoz aşıları toksoiddir) ve toksin
veritabanına homoloji araması en klasik koruyucu antijenleri eler.
"""

from __future__ import annotations

import csv
import functools
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Yazarların peptit düzeyi için önerdiği varsayılan (toxinpred3 -t default).
DEFAULT_THRESHOLD = 0.38


@functools.lru_cache(maxsize=1)
def _cli() -> str | None:
    """CLI yolu. PATH'te arar; bulamazsa çalışan yorumlayıcının yanına bakar.

    venv aktive edilmeden `.venv/bin/python3` ile çalıştırıldığında venv'in bin
    dizini PATH'te olmaz; o durumda shutil.which boş döner ve araç sessizce
    'kurulu değil' sanılırdı.
    """
    hit = shutil.which("toxinpred3")
    if hit:
        return hit
    # sys.executable bir sembolik bağ olabilir (resolve() sistem python'una
    # çıkar), bu yüzden venv kökünü sys.prefix üzerinden alıyoruz.
    for cand in (Path(sys.prefix) / "bin" / "toxinpred3",
                 Path(sys.executable).parent / "toxinpred3"):
        if cand.exists():
            return str(cand)
    return None


@functools.lru_cache(maxsize=1)
def available() -> bool:
    """CLI var VE model kurulu sklearn ile yüklenebiliyor mu?

    Dönüştürülmemiş model 'incompatible dtype' ile patlar; o durumda sessizce
    kullanılabilir saymak yerine False döndürüp çağıranın düşmesini sağlıyoruz.
    """
    if not _cli():
        return False
    return bool(predict(["GILGFVFTLA"], _probe=True))


def predict(peptides: list[str], threshold: float = DEFAULT_THRESHOLD,
            model: int = 1, _probe: bool = False) -> dict[str, dict]:
    """Peptit -> {ml_score, toxic, model}. Hata/eksik -> {}.

    ID eşlemesi SIRA ile yapılır: toxinpred3 girdi başlıklarını kendi
    'Seq_{i}' adlarıyla değiştiriyor. Model 1 çıktısında Sequence kolonu
    bulunduğu için eşleme orada ayrıca doğrulanır.
    """
    if not _cli() or not peptides:
        return {}
    peps = sorted({p for p in peptides if p})
    # Upstream kusuru: TEK dizilik girdide AAC/DPC matrisleri 1-B'ye çöküyor ve
    # np.concatenate(axis=1) AxisError veriyor. site-packages'ı yamalamak yerine
    # kukla bir dizi ekleyip sonucu ayıklıyoruz (pip reinstall'a dayanıklı).
    _pad = "AAAAAAAAAA"
    probe_pad = len(peps) == 1 and peps[0] != _pad
    run_peps = sorted(set(peps) | {_pad}) if probe_pad else peps
    workdir = Path(tempfile.mkdtemp(prefix="toxinpred3_"))
    fin, fout = workdir / "in.seq", workdir / "out.csv"
    fin.write_text("\n".join(run_peps) + "\n", encoding="utf-8")
    try:
        subprocess.run([_cli(), "-i", str(fin), "-o", str(fout),
                        "-m", str(model), "-d", "2", "-t", str(threshold)],
                       cwd=str(workdir), capture_output=True,
                       timeout=7200, text=True)
    except Exception:
        shutil.rmtree(workdir, ignore_errors=True)
        return {}
    out: dict[str, dict] = {}
    try:
        with open(fout, newline="") as fh:
            for i, row in enumerate(csv.DictReader(fh)):
                if i >= len(run_peps):
                    break
                pep = run_peps[i]
                seq_col = (row.get("Sequence") or "").strip()
                if seq_col and seq_col != pep:
                    # sıra varsayımı bozulmuş: sessiz yanlış eşleme yerine boş dön
                    return {}
                score_s = (row.get("Hybrid Score") or row.get("ML Score")
                           or row.get("ML_Score") or "0")
                try:
                    score = float(score_s)
                except ValueError:
                    score = 0.0
                pred = (row.get("Prediction") or "").strip().lower()
                toxic = (pred == "toxin") if pred else score >= threshold
                out[pep] = {"ml_score": round(score, 3), "toxic": toxic,
                            "model": "hybrid" if model == 2 else "et"}
    except FileNotFoundError:
        pass
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    if probe_pad:
        out.pop(_pad, None)
    return out
