"""ToxinPred3'ün ExtraTrees modelini kurulu sklearn sürümüne dönüştürür.

ToxinPred3 (Rathore ve ark. 2024, Comput Biol Med 179:108926,
doi:10.1016/j.compbiomed.2024.108926) modeli eski sklearn ile pickle'lanmış;
bu venv'de sklearn 1.5.x var ve ağaç düğümlerinde 'missing_go_to_left' (u1)
alanını bekliyor → 'node array ... incompatible dtype' hatası.

Dönüşüm mantığı IFNepitope2 ile birebir aynı olduğundan yamalı unpickler
convert_ifn_models'tan yeniden kullanılır (tek kaynak). idempotent.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from convert_ifn_models import convert_one  # noqa: E402


def _model_path() -> str:
    import toxinpred3.python_scripts as ps
    base = os.path.dirname(os.path.dirname(ps.__file__))
    return os.path.join(base, "model", "toxinpred3.0_model.pkl")


def main() -> None:
    print(convert_one(_model_path()))


if __name__ == "__main__":
    main()
