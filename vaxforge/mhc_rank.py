"""MHC %rank -> 0-1 bağlanma skoru: TEK kaynak (tüm omurgalar için).

NEDEN BU MODÜL VAR (2026-11-22'de ölçülerek bulundu)
-----------------------------------------------------
İki ayrı ölçek hatası vardı:

1. `netmhc_local` skoru `1 - rank/100` ile üretiyordu. %rank bir YÜZDELİK
   dilimdir ve biyolojik olarak anlamlı aralığı 0.01-2 arasıdır; bu doğrusal
   dönüşüm o aralığın TAMAMINI 0.98-1.00'e sıkıştırır. Sığır BoLA paneli +
   1500 doğal BoHV-1 9-mer'i ile ölçüldü: eşiği geçen peptitlerde eski skorun
   toplam yayılımı 0.0196 — yani `mhc_score` pratikte SABİT. Adaylık puanı
   (scoring.py) bu bileşeni ağırlıklı ortalamaya katıyordu, dolayısıyla
   %rank 0.003'lük gerçek bir bağlayıcı ile %rank 1.13'lük sınırda bir peptit
   skorda ayırt edilemiyordu (0.99997 vs 0.9887).

2. `iedb` omurgası ise rank'ı hiç kullanmayıp API'nin kendi `score` alanını
   döndürüyordu — netMHCpan'ın yerel skoruyla KARŞILAŞTIRILAMAZ bir ölçek.
   Aynı koşuda iki omurga karışırsa (yerel düşer, IEDB devralır) adaylık
   puanları sessizce kıyaslanamaz hale gelir.

ÇÖZÜM: %rank zaten omurgalar ve aleller ARASI karşılaştırılabilir olması için
tanımlanmış büyüklüktür. O yüzden tek normalize edici burada: log10 ölçekte,
sınıfın KENDİ eşiklerine (güçlü/zayıf) demirlenmiş.

  rank <= strong/10  -> 1.0   (doyma; MHC-I'de %0.05, MHC-II'de %0.2)
  rank >= weak       -> 0.0   (zayıf eşiğin ötesi bağlanma kanıtı değil)
  arası              -> log10 doğrusal

Ölçülen etki (MHC-I, strong=0.5 weak=2.0): eşik geçen doğal peptitlerde
yayılım 0.0196 -> 0.9948. `AAAAARAA` (rank 1.131) 0.9887 -> 0.1545; en iyi
doğal bağlayıcı (rank 0.003) 1.0'da kalır.

Referans: NetMHCpan-4.1 / NetMHCIIpan-4.0 %rank tanımı ve önerilen eşikleri
(Reynisson et al., Nucleic Acids Res 2020, doi:10.1093/nar/gkaa379).
"""

from __future__ import annotations

import math

#: Sınıf başına standart "güçlü bağlayıcı" eşiği (%rank). NetMHCpan ailesinin
#: kendi önerisi; config bunları ezebilir (bkz. config/thresholds.yaml).
DEFAULT_STRONG = {"mhc_i": 0.5, "mhc_ii": 2.0}
DEFAULT_WEAK = {"mhc_i": 2.0, "mhc_ii": 10.0}

#: Doyma tavanı = güçlü eşiğin bu kesri. 0.1 -> MHC-I'de %rank 0.05.
SATURATION_FRACTION = 0.1

METHOD = ("%rank -> skor: log10 ölçek, sınıf eşiklerine demirli "
          "(NetMHCpan %rank tanımı, doi:10.1093/nar/gkaa379)")


def rank_to_score(rank: float, mhc_class: str,
                  strong: float | None = None,
                  weak: float | None = None) -> float:
    """%rank'ı 0-1 bağlanma skoruna çevir (log10, eşiklere demirli).

    `strong`/`weak` verilmezse sınıfın standart eşikleri kullanılır. Eşikler
    bozuksa (strong >= weak, ya da <= 0) skor 0.0 döner — sessiz yanlış skor
    üretmekten iyidir.
    """
    s = DEFAULT_STRONG.get(mhc_class, 0.5) if strong is None else float(strong)
    w = DEFAULT_WEAK.get(mhc_class, 2.0) if weak is None else float(weak)
    try:
        r = float(rank)
    except (TypeError, ValueError):
        return 0.0
    if not (0 < s < w) or r <= 0:
        return 0.0
    hi = s * SATURATION_FRACTION
    if r <= hi:
        return 1.0
    if r >= w:
        return 0.0
    val = (math.log10(w) - math.log10(r)) / (math.log10(w) - math.log10(hi))
    return round(max(0.0, min(1.0, val)), 4)
