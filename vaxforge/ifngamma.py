"""Adım 4 (MHC-II sonrası) — IFN-γ üreten epitop tahmini.

GERÇEK araç: IFNepitope2 (Dhall, Patyal, Raghava, Sci Rep 2024; DOI
10.1038/s41598-024-77957-8). Yardımcı T (HTL / MHC-II) epitoplarının hücresel
Th1 yanıtın itici gücü olan IFN-γ salgılatıp salgılatmadığını tahmin eder —
güncel çok-epitoplu aşı tasarımı literatüründe (Nipah/Salmonella/TBEV 2025)
neredeyse evrensel bir HTL seçim ölçütü.

Hibrit yöntem = ExtraTrees (dipeptit kompozisyonu, DPC) + BLAST düzeltmesi:
  * ml    — ExtraTrees P(IFN-γ indükleyici), 0-1.
  * blast — bilinen indükleyiciye (P_*) benzerse +0.5, indüklemeyene (N_*) -0.5,
            eşleşme yoksa 0 (paketteki blastp + human/mouse db).
  * total = ml + blast; total > eşik (varsayılan 0.49) → indükleyici.

Yalnız İNSAN ve FARE modeli vardır (araç sınırı). Diğer konaklar (sığır/domuz/
tavuk) için IFN-γ tahmini YAPILMAZ (dürüstçe atlanır). Bu yüzden IFN-γ skoru
peptit başına tek değerdir: seçili konaklar arasından insan → yoksa fare kullanılır.

MODEL SÜRÜM NOTU: IFNepitope2 modelleri scikit-learn 0.24 ile pickle'lanmış;
bu venv'de (IApred yüzünden) sklearn 1.5.2 var. sklearn 1.3+ ağaç düğümlerine
'missing_go_to_left' ekledi ve düğüm 'values'ının olasılık (normalize) olmasını
bekliyor (eski sürüm ham sınıf sayısı saklardı). PatchUnpickler ikisini de
çalışma anında düzeltir (idempotent: dönüştürülmüş model üzerinde no-op). Sonuç
1.2.1'in ürettiğiyle birebir aynıdır (doğrulandı; bkz scripts/convert_ifn_models.py).
"""

from __future__ import annotations

import functools
import importlib.util
import os
import pickle
import platform
import subprocess
import tempfile

import numpy as np

METHOD = "GERÇEK (IFNepitope2, ExtraTrees+BLAST, Dhall 2024)"
DEFAULT_THRESHOLD = 0.49
_STD = list("ACDEFGHIKLMNPQRSTVWY")
_ALLOWED = set(_STD)


# --------------------------------------------------------------------------- #
# Model sürüm-köprüsü (sklearn 0.24 pickle → 1.5.2)                             #
# --------------------------------------------------------------------------- #
class _PatchUnpickler(pickle._Unpickler):
    """Eski ExtraTrees ağaçlarını kurulu sklearn'e uyarlayarak yükler.

    pickle._Unpickler opcode'ları class-level `dispatch` dict'inden çağırır;
    metodu override etmek yetmez — dispatch girdisini yeniden bağlarız.
    """

    dispatch = dict(pickle._Unpickler.dispatch)

    def load_build(self):
        state = self.stack[-1]
        if isinstance(state, dict) and "nodes" in state:
            state = dict(state)
            nodes = state["nodes"]
            names = getattr(getattr(nodes, "dtype", None), "names", None)
            if names and "missing_go_to_left" not in names:
                new_dt = np.dtype(nodes.dtype.descr + [("missing_go_to_left", "u1")])
                new_nodes = np.zeros(nodes.shape, dtype=new_dt)
                for n in names:
                    new_nodes[n] = nodes[n]
                new_nodes["missing_go_to_left"] = 0
                state["nodes"] = new_nodes
            if "values" in state:  # ham sayı → olasılık (idempotent: zaten 1'e toplanıyorsa no-op)
                vals = np.asarray(state["values"], dtype=np.float64)
                ssum = vals.sum(axis=-1, keepdims=True)
                ssum[ssum == 0] = 1.0
                state["values"] = vals / ssum
            self.stack[-1] = state
        pickle._Unpickler.load_build(self)

    dispatch[pickle.BUILD[0]] = load_build


def load_extratrees(path: str):
    with open(path, "rb") as fh:
        return _PatchUnpickler(fh).load()


# --------------------------------------------------------------------------- #
# Paket konumu / kaynaklar                                                     #
# --------------------------------------------------------------------------- #
@functools.lru_cache(maxsize=1)
def _base() -> str | None:
    """Kurulu ifnepitope2 paketinin kök dizini (model/blast burada)."""
    spec = importlib.util.find_spec("ifnepitope2.python_scripts")
    if spec is None or not spec.origin:
        return None
    return os.path.dirname(os.path.dirname(spec.origin))


def _blastp() -> str | None:
    base = _base()
    if not base:
        return None
    sub = {"Linux": "linux", "Darwin": "mac", "Windows": "windows"}.get(platform.system(), "linux")
    exe = "blastp.exe" if sub == "windows" else "blastp"
    p = os.path.join(base, "blast_binaries", sub, exe)
    return p if os.path.exists(p) else None


@functools.lru_cache(maxsize=1)
def available() -> bool:
    base = _base()
    return bool(base) and os.path.exists(os.path.join(base, "model", "human_et.pkl"))


@functools.lru_cache(maxsize=2)
def _model(host: str):
    base = _base()
    name = "human_et.pkl" if host == "human" else "mouse_et.pkl"
    return load_extratrees(os.path.join(base, "model", name))


def model_host(host_names: list[str] | None) -> str | None:
    """Seçili konaklardan IFN-γ modeli olanı seç: insan > fare > yok."""
    names = set(host_names or [])
    if "human" in names:
        return "human"
    if "mouse" in names:
        return "mouse"
    return None


# --------------------------------------------------------------------------- #
# Özellikler + tahmin                                                          #
# --------------------------------------------------------------------------- #
def _sanitize(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in _ALLOWED)[:20]


def _dpc(seq: str) -> list[float]:
    """Dipeptit kompozisyonu (400 özellik, DPC1). IFNepitope2 feature_gen ile aynı."""
    n = len(seq) - 1
    if n <= 0:
        return [0.0] * 400
    counts = {a + b: 0 for a in _STD for b in _STD}
    for i in range(n):
        dp = seq[i:i + 2]
        if dp in counts:
            counts[dp] += 1
    return [counts[a + b] / n * 100 for a in _STD for b in _STD]


def _blast_scores(seqs: list[str], host: str) -> dict[str, float]:
    """seq -> ±0.5/0 (IFNepitope2 BLAST_processor mantığı). Hata/araç yoksa hepsi 0."""
    bp, base = _blastp(), _base()
    if not bp or not base:
        return {}
    db = os.path.join(base, "blast_db", "human_db" if host == "human" else "mouse_db")
    idx = {f"q{i}": s for i, s in enumerate(seqs)}
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for qid, s in idx.items():
            fh.write(f">{qid}\n{s}\n")
        qpath = fh.name
    opath = qpath + ".out"
    try:
        subprocess.run([bp, "-task", "blastp-short", "-db", db, "-query", qpath,
                        "-out", opath, "-outfmt", "6", "-evalue", "0.001"],
                       capture_output=True, text=True, timeout=600)
    except Exception:
        for p in (qpath, opath):
            try:
                os.unlink(p)
            except OSError:
                pass
        return {}
    best: dict[str, float] = {}
    try:
        with open(opath) as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) < 2 or parts[0] in best:
                    continue  # ilk (en iyi) hit tutulur
                prefix = parts[1].split("_")[0]
                best[parts[0]] = 0.5 if prefix == "P" else -0.5 if prefix == "N" else 0.0
    except OSError:
        pass
    for p in (qpath, opath):
        try:
            os.unlink(p)
        except OSError:
            pass
    return {idx[q]: v for q, v in best.items()}


def norm(total: float) -> float:
    """total (ml+blast, ~[-0.5, 1.5]) -> 0-1 skor bileşeni."""
    return round(max(0.0, min(1.0, total)), 4)


def predict(seqs: list[str], host: str,
            threshold: float = DEFAULT_THRESHOLD) -> dict[str, dict]:
    """[seq] + konak('human'/'mouse') -> {seq: {ml, blast, total, inducer, norm}}.

    Araç yoksa / konak desteklenmiyorsa boş döner (çağıran atlar).
    """
    if not available() or host not in ("human", "mouse") or not seqs:
        return {}
    uniq = list(dict.fromkeys(_sanitize(s) for s in seqs))
    uniq = [s for s in uniq if s]
    if not uniq:
        return {}
    clf = _model(host)
    X = np.array([_dpc(s) for s in uniq])
    proba = clf.predict_proba(X)[:, -1]
    blast = _blast_scores(uniq, host)
    out: dict[str, dict] = {}
    for s, ml in zip(uniq, proba):
        b = blast.get(s, 0.0)
        total = float(ml) + b
        out[s] = {
            "ml": round(float(ml), 3),
            "blast": b,
            "total": round(total, 3),
            "inducer": total > threshold,
            "norm": norm(total),
        }
    return out
