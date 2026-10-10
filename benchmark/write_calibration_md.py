"""filter_calibration.json -> FILTER_CALIBRATION.md (elle tablo yazmayı bitirir).

Koşum: .venv/bin/python3 benchmark/write_calibration_md.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "benchmark" / "results" / "filter_calibration.json"
DST = ROOT / "benchmark" / "FILTER_CALIBRATION.md"

SET_LABEL = {
    "self9": "insan self 9-mer",
    "self15": "insan self 15-mer",
    "iedb9": "IEDB doğrulanmış epitop (9-mer)",
    "iedb15": "IEDB doğrulanmış epitop (13-17-mer)",
    "tox_protein": "toksin proteini (tam uzunluk) — **duyarlılık**",
    "nontox_protein": "non-toksin proteini (tam uzunluk)",
}
ORDER = list(SET_LABEL)


def main() -> None:
    d = json.loads(SRC.read_text(encoding="utf-8"))
    L: list[str] = []
    L.append("# Filtre kalibrasyonu — kapıların bağımsız kontrol setlerindeki davranışı\n")
    L.append("Üretim: `benchmark/calibrate_filters.py` · veri: "
             "`benchmark/results/filter_calibration.json` · seed: "
             f"{d.get('seed')}\n")
    L.append("Her hücre, o kontrol setinde **işaretlenen** (yani filtrenin eleyeceği) "
             "peptit/protein yüzdesidir. `tox_protein` dışındaki her satırda yüksek "
             "değer = **yanlış-pozitif**: filtrenin elememesi gereken diziyi eliyor.\n")
    L.append("Set boyutları: " + ", ".join(
        f"`{k}`={v}" for k, v in d.get("set_sizes", {}).items()) + "\n")

    L.append("## Kapılar\n")
    gates = d.get("gates", {})
    hdr = "| kapı | " + " | ".join(SET_LABEL[k] for k in ORDER) + " |"
    L.append(hdr)
    L.append("|" + "---|" * (len(ORDER) + 1))
    for gname, g in gates.items():
        cells = []
        for k in ORDER:
            e = g.get(k) or {}
            cells.append("—" if "flagged_pct" not in e else f"%{e['flagged_pct']}")
        L.append(f"| `{gname}` | " + " | ".join(cells) + " |")
    L.append("")

    L.append("## Eşik taraması\n")
    for gname, rows in d.get("sweeps", {}).items():
        L.append(f"### {gname}\n")
        keys = [k for k in ORDER if any(k in r for r in rows)]
        L.append("| eşik | " + " | ".join(SET_LABEL[k] for k in keys) + " |")
        L.append("|" + "---|" * (len(keys) + 1))
        for r in rows:
            L.append(f"| {r['threshold']} | "
                     + " | ".join(f"%{r.get(k, '—')}" for k in keys) + " |")
        L.append("")

    L.append("## Yorum ve karar\n")
    L.append("1. **Toksisite, epitop düzeyinde sert filtre olamaz.** ToxinPred2 "
             "(Model 2 hibrit, eşik 0.6) insan self 9-mer'lerinin ve deneysel "
             "doğrulanmış IEDB epitoplarının çoğunluğunu \"toksik\" işaretliyor; "
             "influenza M1(58-66) GILGFVFTL dahil. Araç tam uzunluk proteinde ise "
             "yüksek duyarlılık ve düşük yanlış-pozitifle çalışıyor. Sorun aracın "
             "kalitesi değil, **uzunluk sınıfı**: ToxinPred2 'toxicity of proteins' "
             "için eğitilmiştir.\n")
    L.append("2. **Epitop uzunluğunda altın standart yok.** Bir toksinden kesilen "
             "9-mer'in kendisi toksik değildir (toksisite katlanmış proteinin "
             "işlevidir), dolayısıyla bu uzunlukta duyarlılık ölçülemez. Ölçülebilen "
             "tek şey yanlış-pozitif oranıdır; bu da mimari kararı dayatıyor: "
             "epitopta **yumuşak işaretleme**, proteinde/MEV'de **sert eleme**.\n")
    L.append("3. **ToxinPred3 peptit düzeyi için doğru araç** (Rathore ve ark. 2024, "
             "Comput Biol Med 179:108926) ama protein düzeyinde ToxinPred2'nin yerini "
             "almaz — tablodaki `tox_protein` satırı bunu gösteriyor.\n")
    L.append("4. **Alerjenite ölçütü uzunluğa göre seçilmeli.** FAO/WHO 6-mer kuralı "
             "9-mer'de özgül; tam proteinde non-alerjen proteinleri de işaretliyor. "
             "MEV konstruktu bir protein olduğu için orada Codex/FAO-WHO birincil "
             "ölçütü (80 aa pencerede ≥%35 kimlik, diamond) kullanılır.\n")
    L.append("5. **Mimikri kapısı sağlam.** İnsan self k-mer'lerinde ~%100, IEDB "
             "epitoplarında ~%0 işaretleme — tasarlandığı gibi ayırıyor ve zaten "
             "yumuşak.\n")
    L.append("\n> Kural: bir sert kapı, bu setlerde kalibre edilmeden üretime girmez.\n")

    DST.write_text("\n".join(L), encoding="utf-8")
    print(f"yazıldı: {DST}")


if __name__ == "__main__":
    main()
