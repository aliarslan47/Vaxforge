"""Epitop dizi karmaşıklığı — düşük karmaşıklıklı peptitleri AŞI ADAYI OLMAKTAN ÇIKAR.

NEDEN (biyolojik gerekçe)
--------------------------
Homopolimerik / düşük karmaşıklıklı peptitler (ör. `AAAAARAA`) bağlanma
tahmininde eşiği geçebilir ama aşı adayı olarak kullanılamaz:

* **TCR tanıması:** TCR, MHC yarığındaki peptidin yukarı bakan yan zincir
  ÇEŞİTLİLİĞİNİ okur. Tek aminoasidin tekrarı ayırt edici bir yüzey sunmaz;
  düşük karmaşıklıklı peptitlere karşı ölçülen T-hücre yanıtları zayıf ve
  tekrar üretilemez niteliktedir.
* **Konak self-benzerliği:** Ala/Gly/Pro/Glu bakımından zengin düşük
  karmaşıklıklı bölgeler konak proteomunda yaygındır (kollajen, prolin-zengin
  bölgeler, poli-Ala uzantıları). Bu da ya merkezi/periferik toleransla yanıt
  yokluğuna ya da otoimmünite riskine işaret eder.
* **Yapısal bağlam:** düşük karmaşıklıklı diziler çoğunlukla içsel düzensiz
  bölgelerden gelir; konformasyonel epitop bütünlüğü ve antikor erişimi
  belirsizdir.
* **Üretilebilirlik:** poli-Ala / poli-Gln uzantıları peptit sentezinde
  agregasyona ve düşük verime yol açar; MEV konstrüktünde tekrar bölgeleri
  ayrıca kodon tekrarı ve ekspresyon sorunu çıkarır.

Bu yüzden düşük karmaşıklık bir SKOR CEZASI değil, **ayrı bir dışlama ölçütü**
olarak uygulanır (bağlanma tahmini doğru olsa bile aday olamaz) ve raporda
GÖRÜNÜR kalır — sessizce atılmaz (bkz. araç dürüstlüğü kuralı: ölçülmemiş ya
da elenmiş olan, olumsuz bulgu gibi gösterilemez).

EŞİKLER NEREDEN GELDİ (2026-11-22, ölçüldü)
--------------------------------------------
BoHV-1 proteomundan 1500 rastgele 9-mer gerçek netMHCpan (sığır BoLA, 12
allel) ile skorlandı; 145'i güçlü bağlayıcı (%rank <= 0.5). Kural adayları bu
küme üzerinde denendi — amaç: homopolimerleri ele, GERÇEK güçlü bağlayıcıyı
kaybetme.

    kural                               elenen güçlü      elenen sonda
    maxrun>=4                           1/145  (0.7%)     12/14
    uniq<=3                             2/145  (1.4%)     14/14
    H<1.5                               3/145  (2.1%)     14/14
    uniq<=4                             7/145  (4.8%)     14/14   (fazla agresif)
    H<2.0                               8/145  (5.5%)     14/14   (fazla agresif)
    SEÇİLEN (maxrun>=4 | uniq<=3 | H<1.5)  3/145 (2.1%)   14/14

Seçilen kural 14 düşük karmaşıklıklı sondanın 14'ünü eler, gerçek güçlü
bağlayıcıların yalnız %2.1'ini kaybeder. Kaybedilen 3 peptit de ölçütlere göre
gerçekten düşük karmaşıklıklıdır (yani yanlış pozitif değil, kuralın kapsamı).

DİKKAT: Ala'dan ZENGİN olmak tek başına eleme sebebi DEĞİL. Ölçümde `AGAPAYVAA`
(%rank 0.008, H=1.88, 5 farklı aa) ölçülmüş EN GÜÇLÜ bağlayıcılardan biriydi ve
bu kural onu GEÇİRİR. Elenen şey tekrar/homopolimer, kompozisyon değil.

Referans yaklaşım: düşük karmaşıklık filtreleme (Wootton & Federhen, Comput
Chem 1993, doi:10.1016/0097-8485(93)85006-X) — SEG/DUST'ın dizi analizinde
standart gerekçesi; burada peptit ölçeğine uyarlandı.
"""

from __future__ import annotations

import math

from .sequtils import sanitize

#: Aynı aminoasidin ardışık tekrar sınırı (>= eleme).
MAX_RUN = 4
#: Peptitte bulunması gereken en az farklı aminoasit sayısı (<= eleme).
MIN_DISTINCT = 4
#: En az Shannon entropisi, bit (< eleme).
MIN_ENTROPY = 1.5

METHOD = ("Dizi karmaşıklığı: Shannon entropisi + homopolimer tekrar uzunluğu + "
          "farklı aminoasit sayısı (BoHV-1 proteomunda kalibre edildi; "
          "yaklaşım: Wootton & Federhen 1993, doi:10.1016/0097-8485(93)85006-X)")


def shannon_entropy(pep: str) -> float:
    """Peptitteki aminoasit dağılımının Shannon entropisi (bit)."""
    s = sanitize(pep)
    n = len(s)
    if n == 0:
        return 0.0
    h = 0.0
    for c in set(s):
        p = s.count(c) / n
        h -= p * math.log2(p)
    return round(h, 3)


def max_run(pep: str) -> int:
    """En uzun ardışık aynı-aminoasit tekrarı."""
    s = sanitize(pep)
    if not s:
        return 0
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def distinct(pep: str) -> int:
    """Peptitteki farklı aminoasit sayısı."""
    return len(set(sanitize(pep)))


def metrics(pep: str) -> dict:
    """Karmaşıklık ölçütleri + karar + (varsa) TÜRKÇE gerekçe.

    Döner: {entropy, max_run, distinct, low_complexity: bool, reason: str|None}
    """
    s = sanitize(pep)
    h, mr, nd = shannon_entropy(s), max_run(s), distinct(s)
    reasons = []
    if mr >= MAX_RUN:
        reasons.append(f"{mr} aa homopolimer tekrar (>= {MAX_RUN})")
    if nd <= MIN_DISTINCT - 1:
        reasons.append(f"yalnız {nd} farklı aminoasit (<= {MIN_DISTINCT - 1})")
    if h < MIN_ENTROPY:
        reasons.append(f"entropi {h} bit (< {MIN_ENTROPY})")
    return {"entropy": h, "max_run": mr, "distinct": nd,
            "low_complexity": bool(reasons),
            "reason": "; ".join(reasons) if reasons else None}


def is_low_complexity(pep: str) -> bool:
    return metrics(pep)["low_complexity"]


def annotate(p) -> None:
    """Peptide nesnesine karmaşıklık metriklerini + dışlama notunu yaz.

    `passed` alanına DOKUNMAZ: `passed` bağlanma kararıdır, bu ayrı bir
    ölçüttür. MEV kaseti `metrics['low_complexity']` olanı dışlar, rapor da
    gerekçesiyle gösterir.
    """
    m = metrics(p.seq)
    p.metrics["seq_entropy"] = m["entropy"]
    p.metrics["seq_max_run"] = m["max_run"]
    p.metrics["seq_distinct"] = m["distinct"]
    p.metrics["low_complexity"] = m["low_complexity"]
    if m["low_complexity"]:
        p.metrics["low_complexity_reason"] = m["reason"]
        p.notes.append(f"düşük karmaşıklık → aşı adaylığından çıkarıldı ({m['reason']})")
    p.methods["complexity"] = METHOD
