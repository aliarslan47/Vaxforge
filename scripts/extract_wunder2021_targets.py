#!/usr/bin/env python3
"""Wunder 2021 (eLife 64166) Ek Tablo 4 → izlenen etiket dosyası.

    .venv/bin/python scripts/extract_wunder2021_targets.py

Okur : data/labels/leptospira/wunder_2021_immunity_correlates/
       elife-64166-supp4-v1.xlsx   (İZLENMEZ — yayıncı ek dosyası, indirilir)
Yazar: vaxforge/refdata/leptospira/wunder2021_targets.tsv  (izlenir)

NEDEN AYRI BİR SCRIPT
---------------------
Yayıncı ek dosyası repoda tutulmaz (bkz. data/labels/.../manifest.json). Ama
etiket seti koşudan koşuya DEĞİŞMEMELİ ve hakem/okur onu görebilmeli. Çözüm:
türetilmiş, atıflı, küçük bir TSV izlenir; onu üreten bu script de izlenir.
Böylece etiketin nereden geldiği ve nasıl türetildiği denetlenebilir.

LİSANS: eLife makaleleri CC-BY 4.0'dır; türetilmiş tablo atıfla yeniden
dağıtılabilir. Atıf TSV başlığına yazılır.

DÜRÜSTLÜK NOTLARI (bu repodaki desen)
-------------------------------------
  * Hiçbir değer hesaplanıp "iyileştirilmez"; sütunlar olduğu gibi taşınır.
  * Etiket KORUMA DEĞİL: "zayıflatılmış aşı modelinde bağışıklığın KORELATI
    olarak proteome array ile saptanmış hedef". Korumayla ilişkilidir, koruma
    ölçümü değildir.
  * Anlamlılık KONAĞA GÖRE ayrıdır (hamster / fare / 'fare tümü'); tek bir
    payda kullanmak yanlış olur. Fare sütunlarında 18 satır BOŞ — o satırlarda
    ölçüm YOK, 'anlamsız' değil; TSV'de boş kalır, 0 yazılmaz.
  * 154 satır ama 152 BENZERSİZ ORF: LIC10464 (ligB) dizide üç ayrı parça
    olarak basılmış (ligB 1-6, 7-12, C-terminal). Satırlar korunur; protein
    düzeyi payda 152'dir.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = (ROOT / "data" / "labels" / "leptospira" /
       "wunder_2021_immunity_correlates" / "elife-64166-supp4-v1.xlsx")
OUT = ROOT / "vaxforge" / "refdata" / "leptospira" / "wunder2021_targets.tsv"

CITATION = ("Wunder EA Jr, Adhikarla H, Hamond C, et al. A live attenuated-vaccine "
            "model confers cross-protective immunity against different species of "
            "the Leptospira genus. eLife. 2021;10:e64166. doi:10.7554/eLife.64166 "
            "(Supplementary file 4; CC-BY 4.0)")

COLS = [("ORF ID", "orf"), ("Product name", "product"), ("Gene", "gene"),
        ("PsortB", "psortb"), ("t-rank Hamster", "t_rank_hamster"),
        ("BHp Hamster", "bhp_hamster"), ("BHp_sig Hamster", "sig_hamster"),
        ("t-rank Mouse", "t_rank_mouse"), ("BHp Mouse", "bhp_mouse"),
        ("BHp_sig mouse", "sig_mouse"), ("t-rank Mouse All", "t_rank_mouse_all"),
        ("BHp Mouse All", "bhp_mouse_all"), ("BHp_sig mouse all", "sig_mouse_all")]


def _cell(v) -> str:
    """Hücreyi metne çevirir. BOŞ hücre BOŞ kalır — 0/False'a çevrilmez."""
    import pandas as pd
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    s = str(v).strip()
    return "" if s in {"-", "_", "nan", "NaN"} else s


def main() -> int:
    if not SRC.is_file():
        print(f"Kaynak dosya yok: {SRC}\n"
              "Önce indir: .venv/bin/python scripts/fetch_lepto_labels.py",
              file=sys.stderr)
        return 1
    import pandas as pd

    d = pd.read_excel(SRC, sheet_name="154 proteins", header=1)
    d = d.dropna(subset=["ORF ID"])
    missing = [src for src, _ in COLS if src not in d.columns]
    if missing:
        print(f"Beklenen sütun(lar) yok, tablo düzeni değişmiş: {missing}",
              file=sys.stderr)
        return 1

    rows = [[_cell(r[src]) for src, _ in COLS] for _, r in d.iterrows()]
    n_orf = len({r[0] for r in rows})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as fh:
        fh.write(f"# Kaynak: {CITATION}\n")
        fh.write("# Etiket: zayıflatılmış aşı modelinde BAĞIŞIKLIĞIN KORELATI olarak "
                 "proteome array ile saptanan hedef. KORUMA ÖLÇÜMÜ DEĞİLDİR.\n")
        fh.write("# sig_* sütunları konağa göre ayrıdır; BOŞ hücre = o konakta "
                 "ÖLÇÜM YOK ('anlamsız' demek değil).\n")
        fh.write(f"# {len(rows)} satır (array özelliği), {n_orf} benzersiz ORF. "
                 "LIC10464 (ligB) üç ayrı parça olarak basılmıştır.\n")
        fh.write(f"# Üretildi: scripts/extract_wunder2021_targets.py "
                 f"(kaynak dosya repoda İZLENMEZ: {SRC.name})\n")
        fh.write("\t".join(dst for _, dst in COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"{OUT.relative_to(ROOT)} yazıldı: {len(rows)} satır, {n_orf} benzersiz ORF.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
