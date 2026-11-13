#!/usr/bin/env python3
"""ToxinPred2'yi bu makinede ÇALIŞIR hâle getirir (üç ayrı kusur).

ToxinPred2 kutudan çıktığı gibi bu ortamda hiç çalışmıyordu ve bu SESSİZ bir
başarısızlıktı: `toxinpred.predict()` boş sözlük dönüyor, huni de bunu
"hiçbir peptit toksik değil" diye okuyordu. Yani ölçülmemiş bir güvenlik
kriteri OLUMLU sonuç gibi görünüyordu. Üç bağımsız neden vardı:

1. `envfile` GELİŞTİRİCİNİN kendi makine yollarını içeriyor
   (`/home/neelam/ncbi-blast-2.2.29+/bin/blastp`,
   `/home/neelam/TOXINPRED2/progs/...`). Dosyanın kendisi "User can change the
   paths according to their machines" diyor — yani elle düzeltilmesi
   bekleniyor. Düzeltilmezse Model 2 (Hybrid) BLAST ve MERCI adımlarında düşer.

2. `toxinpred2.py` pandas 1.x API'si kullanıyor:
     * `DataFrame.append` pandas 2.0'da KALDIRILDI → `pd.concat`.
     * `df.sum(axis=1)` artık metin kolonlarını sessizce atlamıyor →
       `numeric_only=True` (eski pandas davranışının birebir karşılığı).

3. `RF_model` eski scikit-learn ile kaydedilmiş → ayrı betik:
   `scripts/convert_toxinpred_model.py`.

Bu betik 1 ve 2'yi yapar, sonra 3'ü çağırır. Her değiştirilen dosyanın yanına
ilk seferinde `.orig` yedeği bırakılır; IDEMPOTENT.

KULLANIM
--------
    python3 scripts/setup_toxinpred.py
    python3 scripts/setup_toxinpred.py --check    # yalnız bildir (çıkış 1)

NOT — UZUNLUK ALANI: araç çalıştıktan SONRA bile 8-15 aa epitoplarda sert
filtre olarak kullanılmamalı; altı referans setinin tamamında en kısa dizi
35 aa. Bunu pipeline tarafı hallediyor (`toxinpred.MIN_VALIDATED_LENGTH`,
`survival.py`). Ayrıntı: DURUM.md bulgu 22.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TP = REPO / "tools" / "toxinpred2"

# (eski metin, yeni metin, kısa açıklama) — toxinpred2.py için
PANDAS_PATCHES: list[tuple[str, str, str]] = [
    (
        "            df3 = df3.append(df2.loc[df2.Subject==i][0:5]).reset_index(drop=True)",
        "            # pandas 2.0 DataFrame.append'i kaldirdi -> pd.concat (ayni davranis;\n"
        "            # VaxForge kurulum yamasi, bkz. scripts/setup_toxinpred.py)\n"
        "            df3 = pd.concat([df3, df2.loc[df2.Subject==i][0:5]],\n"
        "                            ignore_index=True)",
        "DataFrame.append -> pd.concat",
    ),
    (
        "    df8['Hybrid Score'] = df8.sum(axis=1)",
        "    # pandas 2+ artik string kolonlari sessizce atlamiyor -> numeric_only\n"
        "    # ('Subject' metin kolonu toplama girmesin; eski pandas da onu atliyordu)\n"
        "    df8['Hybrid Score'] = df8.sum(axis=1, numeric_only=True)",
        "sum(axis=1) -> numeric_only=True",
    ),
    (
        "    data_test = np.loadtxt(file_name, delimiter=',')\n    X_test = data_test",
        "    data_test = np.loadtxt(file_name, delimiter=',')\n"
        "    # TEK diziyle cagrildiginda np.loadtxt 1B dizi dondurur; eski sklearn bunu\n"
        "    # tolere ediyordu, 1.x \"Expected 2D array\" ile duser -> tek-peptit girdide\n"
        "    # arac komple cokuyor ve toksisite SESSIZCE olculmemis kaliyordu.\n"
        "    # (VaxForge kurulum yamasi, bkz. scripts/setup_toxinpred.py)\n"
        "    X_test = np.atleast_2d(data_test)",
        "loadtxt 1B -> np.atleast_2d (tek-peptit girdi)",
    ),
]


def _backup(path: Path) -> None:
    b = path.with_name(path.name + ".orig")
    if not b.exists():
        shutil.copy2(path, b)


def fix_envfile(*, check: bool) -> str:
    env = TP / "envfile"
    if not env.exists():
        return "YOK: envfile (ToxinPred2 kurulu değil)"
    blast = shutil.which("blastp")
    if not blast:
        return "ATLANDI: blastp PATH'te yok → Model 2 (Hybrid) kullanılamaz"
    want = (
        "#Path information for BLAST and MERCI commands\n"
        "#Path information of the data required to run the BLAST and MERCI\n"
        "#User can change the paths according to their machines\n"
        f"BLAST:{blast}\n"
        f"BLAST database:{TP}/Database/data\n"
        f"MERCI:{TP}/progs/MERCI_motif_locator.pl\n"
        f"MERCI motif file:{TP}/Database/pos_motif.txt\n"
    )
    if env.read_text() == want:
        return "OK: envfile yolları doğru"
    if check:
        return "DÜZELTİLMELİ: envfile geliştirici makine yollarını içeriyor"
    _backup(env)
    env.write_text(want)
    return "DÜZELTİLDİ: envfile → yerel blastp/MERCI/Database yolları"


def fix_pandas(*, check: bool) -> list[str]:
    script = TP / "toxinpred2.py"
    if not script.exists():
        return ["YOK: toxinpred2.py"]
    txt = script.read_text()
    out, changed = [], False
    for old, new, label in PANDAS_PATCHES:
        if new.splitlines()[-1].strip() in txt:
            out.append(f"OK: {label} (zaten yamalı)")
            continue
        if old not in txt:
            out.append(f"ATLANDI: {label} (beklenen satır bulunamadı)")
            continue
        if check:
            out.append(f"DÜZELTİLMELİ: {label}")
            continue
        if not changed:
            _backup(script)
            changed = True
        txt = txt.replace(old, new, 1)
        out.append(f"DÜZELTİLDİ: {label}")
    if changed:
        script.write_text(txt)
    return out


def fix_model(*, check: bool) -> str:
    conv = REPO / "scripts" / "convert_toxinpred_model.py"
    if not conv.exists():
        return "YOK: convert_toxinpred_model.py"
    cmd = [sys.executable, str(conv)] + (["--check"] if check else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    last = [ln for ln in r.stdout.splitlines() if ln.strip()]
    return "RF_model → " + (last[-1] if last else f"çıkış {r.returncode}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="değişiklik yapma, yalnız bildir (gerekiyorsa çıkış 1)")
    args = ap.parse_args()
    if not TP.is_dir():
        print(f"ToxinPred2 kurulu değil: {TP}")
        return 0
    msgs = [fix_envfile(check=args.check), *fix_pandas(check=args.check),
            fix_model(check=args.check)]
    for m in msgs:
        print("  " + m)
    need = any(m.startswith(("DÜZELTİLMELİ", "DÖNÜŞTÜRÜLMELİ")) or "DÖNÜŞTÜRÜLMELİ" in m
               for m in msgs)
    return 1 if (args.check and need) else 0


if __name__ == "__main__":
    sys.exit(main())
