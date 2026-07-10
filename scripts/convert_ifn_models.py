"""IFNepitope2'nin ExtraTrees modellerini kurulu sklearn sürümüne dönüştürür.

IFNepitope2 (Dhall ve ark. 2024) modelleri scikit-learn 1.2.1 ile pickle'lanmış;
bu venv'de IApred nedeniyle sklearn 1.5.2 var. sklearn 1.3+ ağaç düğümlerine
'missing_go_to_left' (u1) alanını ekledi → eski pickle'lar 'incompatible dtype'
hatasıyla yüklenmiyor. Bu alan yalnız NaN girdide kullanılır; DPC özelliklerinde
NaN olmaz → 0 ile doldurmak tahmini DEĞİŞTİRMEZ (birebir aynı model).

Yamalı saf-Python unpickler ağaç 'nodes' dizisine alanı ekler, model kurulu
sklearn'de yeniden pickle'lanır (yanında .orig yedeği). idempotent: dönüşmüş
model zaten düz pickle.load ile açılıyorsa dokunmaz. pip reinstall sonrası
yeniden çalıştırılabilir.
"""

from __future__ import annotations

import os
import pickle
import shutil

import numpy as np


def _pkgdir() -> str:
    # NOT: ifnepitope2.py'yi import ETME — modül düzeyinde argparse çalıştırıyor.
    import ifnepitope2.python_scripts as ps
    return os.path.dirname(os.path.dirname(ps.__file__))


class _PatchUnpickler(pickle._Unpickler):
    """BUILD sırasında eski ağaç 'nodes' dizisine missing_go_to_left ekler.

    NOT: pickle._Unpickler opcode'ları class-level `dispatch` dict'inden çağırır;
    metodu override etmek yetmez — dispatch girdisini de yeniden bağlarız.
    """

    dispatch = dict(pickle._Unpickler.dispatch)

    def load_build(self):
        state = self.stack[-1]
        if isinstance(state, dict) and "nodes" in state:
            state = dict(state)
            nodes = state["nodes"]
            names = getattr(getattr(nodes, "dtype", None), "names", None)
            if names and "missing_go_to_left" not in names:
                # (1) sklearn 1.3+ ağaç düğümlerine missing_go_to_left (u1) ekledi.
                new_dt = np.dtype(nodes.dtype.descr + [("missing_go_to_left", "u1")])
                new_nodes = np.zeros(nodes.shape, dtype=new_dt)
                for n in names:
                    new_nodes[n] = nodes[n]
                new_nodes["missing_go_to_left"] = 0
                state["nodes"] = new_nodes
            # (2) Eski sklearn (0.24) düğüm 'values'ını HAM SINIF SAYISI olarak
            # saklıyordu; predict_proba tahmin anında normalize ederdi. sklearn 1.3+
            # ise value'nun ZATEN normalize olmasını bekler → ham sayılar 1'den
            # büyük "olasılık" verir. Düğüm başına normalize et (sınıf ekseni) →
            # 1.2.1'in ürettiğiyle birebir aynı, ağaçlar eşit ağırlıklı ortalanır.
            if "values" in state:
                vals = np.asarray(state["values"], dtype=np.float64)
                ssum = vals.sum(axis=-1, keepdims=True)
                ssum[ssum == 0] = 1.0
                state["values"] = vals / ssum
            self.stack[-1] = state
        pickle._Unpickler.load_build(self)

    dispatch[pickle.BUILD[0]] = load_build


def _loads_cleanly(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            pickle.load(fh)
        return True
    except Exception:
        return False


def convert_one(path: str) -> str:
    # Pristine kaynak: .orig yedeği (varsa) — yoksa path'i yedekle.
    backup = path + ".orig"
    if not os.path.exists(backup):
        if _loads_cleanly(path):
            return f"OK (zaten uyumlu, yedek yok): {path}"
        shutil.copy2(path, backup)
    with open(backup, "rb") as fh:
        model = _PatchUnpickler(fh).load()
    with open(path, "wb") as fh:
        pickle.dump(model, fh)
    assert _loads_cleanly(path), "dönüşüm sonrası yine yüklenemedi"
    return f"DÖNÜŞTÜRÜLDÜ (nodes+values normalize): {path}"


def main() -> None:
    base = _pkgdir()
    for name in ("human_et.pkl", "mouse_et.pkl"):
        print(convert_one(os.path.join(base, "model", name)))


if __name__ == "__main__":
    main()
