"""Antijenite — z-descriptor Auto-Cross-Covariance (ACC) yöntemi.

VaxiJen otomatikleştirilemediği için (web 403, indirilebilir sürüm yok),
VaxiJen'in YAYINLANMIŞ ÖZELLİK YÖNTEMİNİ uygular: Sandberg (1998) z-skalaları +
auto/cross-covariance dönüşümü. Üstüne kompakt bir lojistik regresyon eğitir:
  pozitif = küratörlü bakteriyel virülans faktörleri (VFDB setA, yoksa UniProt
            KW-0843) — antijenik/immünojenik bakteri proteinleri
  negatif = insan Swiss-Prot (konağa yabancı olmayan -> antijenik değil)

Hangi pozitif setle eğitildiği modelin içine yazılır (`model_info()`); rapor
etiketi bunu okur, kaynak adı varsayılmaz.

DÜRÜST ETİKET: Bu, VaxiJen'in ÖZELLİK yöntemidir; VaxiJen'in resmi eğitilmiş SVM
modeli DEĞİLDİR. Yaklaşık bir antijenite indeksidir. Model bir kez eğitilip
önbelleğe alınır (tools/db/antigen_acc_model.pkl).
"""

from __future__ import annotations

import functools
import gzip
import pickle
from pathlib import Path

import numpy as np

# Sandberg (1998) z-skalaları (z1..z5): lipofiliklik, sterik, polarite, +2
Z5 = {
    "A": (0.24, -2.32, 0.60, -0.14, 1.30), "R": (3.52, 2.50, -3.50, 1.99, -0.17),
    "N": (3.05, 1.62, 1.04, -1.15, 1.61), "D": (3.98, 0.93, 1.93, -2.46, 0.75),
    "C": (0.84, -1.67, 3.71, 0.18, -2.65), "Q": (1.75, 0.50, -1.44, -1.34, 0.66),
    "E": (3.11, 0.26, -0.11, -3.04, -0.25), "G": (2.05, -4.06, 0.36, -0.82, -0.38),
    "H": (2.47, 1.95, 0.26, 3.90, 0.09), "I": (-3.89, -1.73, -1.71, -0.84, 0.26),
    "L": (-4.28, -1.30, -1.49, -0.72, 0.84), "K": (2.29, 0.89, -2.49, 1.49, 0.31),
    "M": (-2.85, -0.22, 0.47, 1.94, -0.98), "F": (-4.22, 1.94, 1.06, 0.54, -0.62),
    "P": (-1.66, 0.27, 1.84, 0.70, 2.00), "S": (2.39, -1.07, 1.15, -1.39, 0.67),
    "T": (0.75, -2.18, -1.12, -1.46, -0.40), "W": (-4.36, 3.94, 0.59, 3.44, -1.59),
    "Y": (-2.54, 2.44, 0.43, 0.04, -1.47), "V": (-2.59, -2.64, -1.54, -0.85, -0.02),
}
_LAG = 8
_DB = Path(__file__).resolve().parent.parent / "tools" / "db"
_MODEL = _DB / "antigen_acc_model.pkl"

# POZİTİF SET adayları, sırayla (küratörlü bakteriyel virülans proteinleri).
# discovery._VIRDB_SOURCES ile aynı mantık: VFDB varsa o, yoksa UniProt KW-0843.
_POS_SOURCES: list[tuple[Path, str]] = [
    (_DB / "VFDB_setA_pro.fas.gz", "VFDB setA"),
    (_DB / "uniprot_virulence_kw0843.fasta.gz", "UniProt KW-0843"),
]

# NEGATİF SET adayları (konak proteomu). İKİ BİÇİM: host_db otomatik indirmesi
# SIKIŞTIRILMAMIŞ `human_sprot_auto.fasta` yazıyor, eski sabit yol ise yalnız
# `human_sprot.fasta.gz`e bakıyordu. O dosya hiç var olmadığı için model, VFDB
# indirilmiş olsa bile ASLA eğitilemezdi — yani ACC/VaxiJen'in eksikliği tek
# başına VFDB'ye bağlı değildi (2026-11-17'de bulundu).
_NEG_SOURCES: list[tuple[Path, str]] = [
    (_DB / "human_sprot.fasta.gz", "insan Swiss-Prot"),
    (_DB / "human_sprot_auto.fasta", "insan Swiss-Prot (otomatik)"),
]


def _first_existing(cands: list[tuple[Path, str]]) -> tuple[Path, str] | None:
    for path, label in cands:
        if path.exists():
            return path, label
    return None


def training_sets() -> tuple[tuple[Path, str], tuple[Path, str]] | None:
    """Modelin GERÇEKTEN eğitilebileceği (pozitif, negatif) kaynaklar ya da None.

    Rapordaki yöntem etiketi bunu kullanır: hangi pozitif set koştuysa o yazılır
    (ör. 'UniProt KW-0843'), VFDB adı varsayılmaz.
    """
    pos, neg = _first_existing(_POS_SOURCES), _first_existing(_NEG_SOURCES)
    return (pos, neg) if pos and neg else None


def acc_features(seq: str, lag: int = _LAG) -> np.ndarray:
    """z-descriptor auto-cross-covariance -> 5*5*lag boyutlu vektör."""
    z = np.array([Z5[c] for c in seq.upper() if c in Z5])  # (N,5)
    n = len(z)
    feats = []
    for l in range(1, lag + 1):
        if n > l:
            a, b = z[:n - l], z[l:]
            cov = (a[:, :, None] * b[:, None, :]).sum(axis=0) / (n - l)  # (5,5)
        else:
            cov = np.zeros((5, 5))
        feats.append(cov.ravel())
    return np.concatenate(feats)


def _open_text(path: Path):
    """FASTA'yı aç — .gz olsun olmasın (eski hâli yalnız gzip açıyordu)."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="latin-1")
    return open(path, "rt", encoding="latin-1")


def _iter_seqs(path: Path, limit: int, min_len: int = 50):
    out = []
    with _open_text(path) as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                s = "".join(seq)
                if len(s) >= min_len:
                    out.append(s)
                    if len(out) >= limit:
                        return out
                seq = []
            else:
                seq.append(line.strip())
    return out


def _train():
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    sets = training_sets()
    if sets is None:
        return None
    (pos_path, pos_label), (neg_path, neg_label) = sets
    pos = _iter_seqs(pos_path, 500)
    neg = _iter_seqs(neg_path, 500)
    if not pos or not neg:
        return None
    X = np.array([acc_features(s) for s in pos + neg])
    y = np.array([1] * len(pos) + [0] * len(neg))
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=1000).fit(scaler.transform(X), y)
    # Eğitim kaynakları modelin İÇİNE yazılır: önbellekten yüklenen bir model de
    # hangi pozitif setle eğitildiğini söyleyebilsin (rapor etiketi buna bakar).
    model = {"scaler": scaler, "clf": clf, "n_pos": len(pos), "n_neg": len(neg),
             "pos_label": pos_label, "neg_label": neg_label}
    _MODEL.parent.mkdir(parents=True, exist_ok=True)
    _MODEL.write_bytes(pickle.dumps(model))
    return model


def model_info() -> dict | None:
    """Eğitilmiş modelin künyesi: pozitif/negatif set etiketi + örnek sayıları."""
    m = _model()
    if not m:
        return None
    return {k: m.get(k) for k in ("pos_label", "neg_label", "n_pos", "n_neg")}


@functools.lru_cache(maxsize=1)
def _model():
    if _MODEL.exists():
        try:
            return pickle.loads(_MODEL.read_bytes())
        except Exception:
            pass
    if training_sets() is not None:
        try:
            return _train()
        except Exception:
            return None
    return None


def available() -> bool:
    return _model() is not None


def predict(seq: str) -> float:
    """Protein -> antijenite olasılığı 0-1 (ACC + lojistik model)."""
    m = _model()
    if m is None:
        return 0.0
    x = m["scaler"].transform(acc_features(seq).reshape(1, -1))
    return round(float(m["clf"].predict_proba(x)[0, 1]), 3)
