#!/usr/bin/env python3
"""ToxinPred2'nin RF_model'ini kurulu scikit-learn sürümüne dönüştürür.

SORUN
-----
ToxinPred2 (Sharma ve ark. 2022) `RF_model` dosyası eski bir scikit-learn ile
(joblib üzerinden) kaydedilmiş. sklearn 1.3+ karar-ağacı düğümlerine
`missing_go_to_left` (u1) alanını ekledi, bu yüzden kurulu 1.5.2 eski diziyi
reddediyor:

    ValueError: node array from the pickle has an incompatible dtype

Sonuç SESSİZ DEĞİL ama GÖRÜNMEZDİ: `toxinpred2.py` çöküyor, CSV üretilmiyor,
`toxinpred.predict()` BOŞ sözlük dönüyor ve pipeline bunu "hiçbir peptit toksik
değil" diye okuyordu. Yani ölçülmemiş bir güvenlik kriteri, OLUMLU sonuç gibi
görünüyordu — projenin dürüstlük ilkesine aykırı.

NEDEN DÖNÜŞÜM TAHMİNİ DEĞİŞTİRMEZ
---------------------------------
`missing_go_to_left` YALNIZ girdide NaN varken kullanılır. ToxinPred2'nin
özellikleri amino-asit kompozisyonu (AAC) — tanım gereği NaN içermez. Alanı 0
ile doldurmak orijinal modelin davranışını birebir korur.

Eski sklearn (0.24) ayrıca düğüm `values` dizisini HAM SINIF SAYISI olarak
saklıyordu ve `predict_proba` tahmin anında normalize ediyordu; 1.3+ ise
value'nun ZATEN normalize olmasını bekler. Bu yüzden düğüm başına sınıf
ekseninde normalize edilir (aynı gerekçe: `scripts/convert_ifn_models.py`).

KULLANIM
--------
    python3 scripts/convert_toxinpred_model.py
    python3 scripts/convert_toxinpred_model.py --check   # yalnız bildir

İlk dönüşümde `RF_model.orig` yedeği bırakılır ve dönüşüm HER ZAMAN o pristine
yedekten yapılır → idempotent, tekrar tekrar koşturulabilir.
"""
from __future__ import annotations

import argparse
import pickle
import shutil
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
TARGETS = [REPO / "tools" / "toxinpred2" / "RF_model"]


def _patch_state(state: dict) -> dict:
    """Eski ağaç durumunu kurulu sklearn'in beklediği biçime çevirir."""
    state = dict(state)
    nodes = state.get("nodes")
    names = getattr(getattr(nodes, "dtype", None), "names", None)
    if names and "missing_go_to_left" not in names:
        new_dt = np.dtype(nodes.dtype.descr + [("missing_go_to_left", "u1")])
        new_nodes = np.zeros(nodes.shape, dtype=new_dt)
        for n in names:
            new_nodes[n] = nodes[n]
        new_nodes["missing_go_to_left"] = 0       # NaN yok → davranış aynı
        state["nodes"] = new_nodes
    if "values" in state:
        vals = np.asarray(state["values"], dtype=np.float64)
        ssum = vals.sum(axis=-1, keepdims=True)
        ssum[ssum == 0] = 1.0
        state["values"] = vals / ssum             # ham sayı → olasılık
    return state


def _make_unpickler_class():
    """joblib'in NumpyUnpickler'ını yamalar (dosya joblib ile yazıldı).

    pickle._Unpickler opcode'ları CLASS düzeyindeki `dispatch` dict'inden
    çağırır; metodu override etmek yetmez, dispatch girdisini de yeniden
    bağlamak gerekir.
    """
    from joblib.numpy_pickle import NumpyUnpickler

    class _Patched(NumpyUnpickler):
        dispatch = dict(NumpyUnpickler.dispatch)

        def load_build(self):
            st = self.stack[-1]
            if isinstance(st, dict) and "nodes" in st:
                self.stack[-1] = _patch_state(st)
            NumpyUnpickler.load_build(self)

        dispatch[pickle.BUILD[0]] = load_build

    return _Patched


def _backfill_attrs(model) -> list[str]:
    """Yeni sklearn'in beklediği ama eski pickle'da OLMAYAN öznitelikleri ekler.

    Eski sürümlerde var olmayan her alan, o sürümün DAVRANIŞINA denk gelen
    varsayılanla doldurulur — yani model aynı tahmini üretmeye devam eder:

      * monotonic_cst (sklearn 1.4+) = None → monotonluk kısıtı yok (0.22'de
        böyle bir özellik hiç yoktu).
      * n_features_in_ (1.0+) = eski `n_features_` değeri.
      * feature_names_in_ bilerek EKLENMEZ: yokluğu "ad bilgisi yok" demektir,
        uydurma ad koymak yanlış olur.

    Dönüş: eklenen alanların listesi (şeffaflık için raporlanır).
    """
    added: list[str] = []
    trees = list(getattr(model, "estimators_", []) or [])
    for obj in [model, *trees]:
        if not hasattr(obj, "monotonic_cst"):
            obj.monotonic_cst = None
            added.append(f"{type(obj).__name__}.monotonic_cst")
        if not hasattr(obj, "n_features_in_"):
            nf = getattr(obj, "n_features_", None)
            if nf is not None:
                obj.n_features_in_ = nf
                added.append(f"{type(obj).__name__}.n_features_in_")
    # aynı ad birden çok ağaçta eklenir → tekilleştir, sayıyı koru
    uniq = sorted(set(added))
    return [f"{a} (×{added.count(a)})" if added.count(a) > 1 else a for a in uniq]


def _loads_cleanly(path: Path) -> bool:
    """Model yükleniyor VE gerçekten tahmin üretiyor mu?

    Yalnız `joblib.load` başarısına bakmak YETMEZ: dtype düzeltildikten sonra
    yükleme geçiyor ama `predict_proba` eksik öznitelik yüzünden
    (monotonic_cst) hâlâ düşüyordu. Bu yüzden minik bir gerçek tahmin denenir.
    """
    try:
        import joblib
        m = joblib.load(path)
        nf = getattr(m, "n_features_in_", None) or getattr(m, "n_features_", None)
        if nf:
            m.predict_proba(np.zeros((1, int(nf))))
        return True
    except Exception:
        return False


def _load_patched(path: Path):
    """Yamalı unpickler ile yükler (joblib'in imza farklarına dayanıklı)."""
    cls = _make_unpickler_class()
    with open(path, "rb") as fh:
        try:
            up = cls(filename=str(path), file_handle=fh,
                     ensure_native_byte_order="auto")
        except TypeError:                 # eski joblib: bu parametre yok
            up = cls(filename=str(path), file_handle=fh)
        return up.load()


def convert_one(path: Path, *, check: bool) -> str:
    if not path.exists():
        return f"YOK (araç kurulu değil): {path}"
    backup = path.with_name(path.name + ".orig")
    if _loads_cleanly(path):
        return f"OK (zaten uyumlu): {path.name}"
    if check:
        return f"DÖNÜŞTÜRÜLMELİ: {path.name} (kurulu sklearn ile yüklenemiyor)"
    if not backup.exists():
        shutil.copy2(path, backup)
    model = _load_patched(backup)
    added = _backfill_attrs(model)
    import joblib
    joblib.dump(model, path)
    if not _loads_cleanly(path):
        shutil.copy2(backup, path)        # geri al, yanlış model bırakma
        return f"BAŞARISIZ (geri alındı): {path.name}"
    detail = ("nodes+values" if not added
              else "nodes+values, eklenen alan: " + ", ".join(added))
    return f"DÖNÜŞTÜRÜLDÜ ({detail}): {path.name}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="değişiklik yapma, yalnız bildir (gerekiyorsa çıkış 1)")
    args = ap.parse_args()
    msgs = [convert_one(t, check=args.check) for t in TARGETS]
    for m in msgs:
        print(m)
    need = any(m.startswith(("DÖNÜŞTÜRÜLMELİ", "BAŞARISIZ")) for m in msgs)
    return 1 if (args.check and need) else (1 if any(m.startswith("BAŞARISIZ") for m in msgs) else 0)


if __name__ == "__main__":
    sys.exit(main())
