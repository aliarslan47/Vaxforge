#!/usr/bin/env python3
"""DTU araç sarmalayıcılarındaki satıcı YER TUTUCU yollarını gerçek kuruluma çevirir.

NEDEN BU BETİK VAR
------------------
netMHCpan / netMHCIIpan / netCTL / BepiPred paketleri, kurulum sırasında ELLE
düzenlenmek üzere satıcı yer tutucu yoluyla gelir, örneğin:

    setenv NMHOME /tools/src/netMHCpan-4.2
    setenv BP     /usr/cbs/packages/bepipred/1.0/bepipred-1.0
    setenv CTLHOME /var/www/packages/netCTL/1.2a/netCTL-1.2a
    setenv TMPDIR  /scratch

Bu satır düzeltilmezse araç HİÇ çalışmaz ama çoğu zaman SESSİZCE başarısız olur:
sarmalayıcı dosyası vardır, `available()` True döner, tahmin fonksiyonu boş
sonuç döndürür ve pipeline bunu "bulgu yok" sanır. VaxForge artık bu durumu
yakalıyor (netmhc_local._wrapper_home_ok, bepipred._configured,
netctl.available duman testi, vaxforge/tool_status.py) — ama en iyisi yolu
baştan doğru kurmaktır.

`tools/` dizini .gitignore'da (lisanslı araçlar repoya girmez), dolayısıyla bu
düzeltmeler sürüm kontrolünde SAKLANMAZ. Araçları yeniden çıkarınca bu betiği
koşturmak gerekir.

KULLANIM
--------
    python3 scripts/configure_dtu_wrappers.py            # düzelt
    python3 scripts/configure_dtu_wrappers.py --check     # yalnız bildir (çıkış 1)

Her dosyanın ilk değişiminde yanına `<ad>.orig` yedeği bırakılır. Betik
IDEMPOTENT: zaten doğru olan yolu tekrar yazmaz.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TOOLS = REPO / "tools"

# (sarmalayıcı glob, ortam değişkeni, hedefin sarmalayıcıya göre konumu)
# "." = sarmalayıcıyı içeren dizin;  "tmp" = o dizinin altındaki tmp/
WRAPPERS: list[tuple[str, str, str]] = [
    ("netMHCpan-*/netMHCpan",     "NMHOME",  "."),
    ("netMHCIIpan-*/netMHCIIpan", "NMHOME",  "."),
    ("netCTL-*/netCTL",           "CTLHOME", "."),
    ("netCTL-*/netCTL",           "TMPDIR",  "tmp"),
    ("bepipred-*/bepipred",       "BP",      "."),
    ("bepipred-*/bepipred",       "TMP",     "tmp"),
]


def _target(wrapper: Path, where: str) -> Path:
    base = wrapper.parent
    return base if where == "." else base / where


def fix_one(wrapper: Path, var: str, where: str, *, check: bool) -> str | None:
    """Tek bir `setenv <var> <yol>` satırını düzeltir. Dönüş: açıklama ya da None."""
    try:
        txt = wrapper.read_text(errors="ignore")
    except Exception as e:
        return f"OKUNAMADI {wrapper}: {type(e).__name__}"
    pat = re.compile(rf"^(\s*setenv\s+{var}\s+)(\S+)", re.M)
    m = pat.search(txt)
    if not m:
        return None                      # bu araçta bu değişken yok — sorun değil
    current = Path(m.group(2))
    want = _target(wrapper, where)
    if current == want and current.is_dir():
        return None                      # zaten doğru
    if check:
        state = "dizin değil" if not current.is_dir() else "yanlış yol"
        return f"DÜZELTİLMELİ {wrapper.name}: {var}={current} ({state}) → {want}"
    want.mkdir(parents=True, exist_ok=True)
    backup = wrapper.with_suffix(wrapper.suffix + ".orig")
    if not backup.exists():
        shutil.copy2(wrapper, backup)
    wrapper.write_text(pat.sub(lambda mm: mm.group(1) + str(want), txt, count=1))
    return f"DÜZELTİLDİ {wrapper.name}: {var} → {want}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="değişiklik yapma, yalnız bildir (sorun varsa çıkış 1)")
    args = ap.parse_args()

    if not TOOLS.is_dir():
        print(f"tools/ yok: {TOOLS} — lisanslı araçlar kurulmamış.")
        return 0

    msgs: list[str] = []
    seen: set[Path] = set()
    for glob, var, where in WRAPPERS:
        for wrapper in sorted(TOOLS.glob(glob)):
            if wrapper.name.endswith((".orig", ".orig2")) or not wrapper.is_file():
                continue
            seen.add(wrapper)
            msg = fix_one(wrapper, var, where, check=args.check)
            if msg:
                msgs.append(msg)

    if not seen:
        print("Hiç DTU sarmalayıcısı bulunamadı (tools/ altında kurulu paket yok).")
        return 0
    print(f"İncelenen sarmalayıcı: {len(seen)}")
    for m in msgs:
        print("  " + m)
    if not msgs:
        print("  Tüm yollar doğru (değişiklik gerekmedi).")
    # --check modunda sorun bulunduysa CI/kurulum fark etsin diye 1 dön
    return 1 if (args.check and msgs) else 0


if __name__ == "__main__":
    sys.exit(main())
