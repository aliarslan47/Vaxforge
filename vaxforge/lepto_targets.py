"""Leptospira protein düzeyi etiket seti — Wunder 2021 bağışıklık korelatları.

NEDEN BU MODÜL VAR
------------------
`lepto_validation.py` PEPTİT düzeyinde, 19 kayıtlı küçük bir set. Ama bu
hattın ürettiği asıl şey PROTEİN düzeyi aday sıralamasıdır (hangi antijen
seçilsin?) ve o düzeyde Leptospira'da ölçülecek bir payda yoktu: IEDB bu
organizmada pratikte boş (2026-11-22'de canlı ölçüldü: cins 2, L. interrogans
14 kayıt), yani recall ölçümü anlamsız.

Wunder 2021 (eLife 64166) zayıflatılmış canlı aşı modelinde ÇAPRAZ KORUMA ile
ilişkilendirilmiş 154 proteome-array hedefi verir. Elimizdeki, korumaya EN
YAKIN protein düzeyi etikettir.

ETİKETİN NE OLDUĞU — VE NE OLMADIĞI
-----------------------------------
  * Etiket = "korunan hayvanların serumunda array üzerinde anlamlı reaktivite".
    Yani KORELAT. Protein başına KORUMA ölçümü DEĞİLDİR: tek tek bu 154
    proteinin aşı olarak koruduğu gösterilmemiştir.
  * Bu ayrım bu organizmada hayati: Barbosa 2024'te ters vaksinolojiyle bulunan
    bir antijen (MCP) güçlü IgG uyardı, KORUMADI, aşılanan grup hastalığı
    kontrolden DAHA HIZLI geliştirdi. Yani yüksek reaktivite iyi aday
    demek değildir; bkz. `correlate.py`.
  * Anlamlılık KONAĞA GÖRE ayrıdır. Benzersiz PROTEİN olarak payda:
    hamster 131 (133 satır), fare 54 (56 satır), 'fare tümü' 13; fare
    sütunlarında 18 satırda ÖLÇÜM YOK. Tek paydayla recall hesaplamak bu
    yüzden yasaktır — `denominator(host)` konağı zorunlu alır ve satır değil
    BENZERSİZ PROTEİN sayar (satır saymak paydayı şişirip recall'ı yapay
    olarak düşürür).
    (Aynı ders peptit düzeyinde de çıktı: Mariano 2023 Pep2 hamsterda özgül
    pozitif, insanda özgül değil.)
  * Payda, proteomun TAMAMI taranarak kurulmuştur (array = 3359 ORF, Lessa-
    Aquino 2015), ama bu modül yalnız POZİTİF hedefleri içerir. Dolayısıyla
    buradan recall (duyarlılık) ölçülür, ÖZGÜLLÜK ÖLÇÜLEMEZ. Özgüllük için
    array'in negatif kısmı gerekir — henüz eklenmedi.

VERİ KAYNAĞI
------------
`refdata/leptospira/wunder2021_targets.tsv` — türetilmiş, atıflı, izlenen dosya.
Yayıncı ek dosyası repoda tutulmaz; TSV'yi `scripts/extract_wunder2021_targets.py`
üretir. 154 satır = array özelliği, 152 benzersiz ORF (LIC10464/ligB üç parça).

ORF kimlikleri L. interrogans sv Copenhageni Fiocruz L1-130 locus_tag'leridir
(LIC…). Bir koşunun proteinleriyle eşleştirme locus_tag üzerindendir; eşleşme
bulunamazsa bu DÜRÜSTÇE söylenir, ad benzerliğiyle tahmin YAPILMAZ.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

METHOD = "GERÇEK (hakemli literatür, proteome array bağışıklık korelatı)"
VERSION = 1
COMPILED = "2026-11-23"

_TSV = Path(__file__).resolve().parent / "refdata" / "leptospira" / "wunder2021_targets.tsv"

CITATION = ("Wunder EA Jr, Adhikarla H, Hamond C, et al. eLife. 2021;10:e64166. "
            "doi:10.7554/eLife.64166 (Supplementary file 4)")

#: Etiketin okuması — raporda bu cümle AYNEN geçmelidir.
READOUT = ("zayıflatılmış canlı aşı modelinde çapraz korumayla İLİŞKİLENDİRİLMİŞ "
           "proteome array reaktivitesi (korelat); protein başına koruma ölçümü değil")

#: Hangi konaklarda anlamlılık sütunu var. 'mouse_all' makalenin kendi
#: ayrımıdır (tüm fare grupları birlikte); yeniden adlandırılmadı.
HOSTS = ("hamster", "mouse", "mouse_all")

_LOCUS_RE = re.compile(r"\bLIC[_]?(\d{5})\b", re.I)


@dataclass(frozen=True)
class Target:
    """Tek bir array hedefi (satır). Ölçülmemiş alan None kalır."""

    orf: str
    product: str
    gene: str
    psortb: str
    sig: dict[str, bool | None]      # konak -> anlamlı mı (None = ölçüm yok)
    bhp: dict[str, float | None]     # konak -> BH düzeltilmiş p
    t_rank: dict[str, int | None]    # konak -> t-istatistiği sırası

    def significant_in(self, host: str) -> bool | None:
        """O konakta anlamlı mı? Ölçülmediyse None (False DEĞİL)."""
        return self.sig.get(host)


def _num(s: str):
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _flag(s: str) -> bool | None:
    """'1'/'0' → bool; BOŞ → None. Ölçüm yokluğu negatife çevrilmez."""
    if s == "":
        return None
    v = _num(s)
    return None if v is None else bool(v)


@lru_cache(maxsize=1)
def load() -> tuple[Target, ...]:
    """TSV'yi okur. Dosya yoksa BOŞ döner (uydurma veri yok)."""
    if not _TSV.is_file():
        return ()
    out: list[Target] = []
    header: list[str] = []
    for line in _TSV.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if not header:
            header = parts
            continue
        row = dict(zip(header, parts + [""] * (len(header) - len(parts))))
        out.append(Target(
            orf=row["orf"].strip(),
            product=row.get("product", ""),
            gene=row.get("gene", ""),
            psortb=row.get("psortb", ""),
            sig={"hamster": _flag(row.get("sig_hamster", "")),
                 "mouse": _flag(row.get("sig_mouse", "")),
                 "mouse_all": _flag(row.get("sig_mouse_all", ""))},
            bhp={"hamster": _num(row.get("bhp_hamster", "")),
                 "mouse": _num(row.get("bhp_mouse", "")),
                 "mouse_all": _num(row.get("bhp_mouse_all", ""))},
            t_rank={h: (int(v) if (v := _num(row.get(f"t_rank_{h}", ""))) is not None
                        else None) for h in HOSTS},
        ))
    return tuple(out)


def available() -> bool:
    """Etiket seti bu kurulumda var mı?"""
    return bool(load())


def n_features() -> int:
    """Array özelliği (satır) sayısı — 154."""
    return len(load())


def n_proteins() -> int:
    """BENZERSİZ protein sayısı — 152 (LIC10464 üç parça olarak basılmış)."""
    return len({t.orf for t in load()})


def denominator(host: str) -> tuple[str, ...]:
    """`host` konağında ANLAMLI olan benzersiz ORF kümesi (recall paydası).

    Konak zorunludur ve 'hepsi' diye bir seçenek yoktur: anlamlılık konağa
    göre değişir (hamster 133, fare 56, 'fare tümü' 13) ve fare sütunlarında
    18 satırda ölçüm yoktur. Tek payda kullanmak ölçümü bozar.
    """
    if host not in HOSTS:
        raise ValueError(f"bilinmeyen konak {host!r}; seçenekler: {HOSTS}")
    return tuple(sorted({t.orf for t in load() if t.significant_in(host) is True}))


def unmeasured(host: str) -> tuple[str, ...]:
    """O konakta ÖLÇÜLMEMİŞ ORF'ler — 'anlamsız' sayılanlardan ayrı tutulur."""
    if host not in HOSTS:
        raise ValueError(f"bilinmeyen konak {host!r}; seçenekler: {HOSTS}")
    return tuple(sorted({t.orf for t in load() if t.significant_in(host) is None}))


def locus_tags(text: str | None) -> set[str]:
    """Bir protein başlığı/tanımından LIC locus_tag'lerini çıkarır.

    Yalnız açık LIC##### kalıbı kabul edilir. Ürün adı benzerliğiyle eşleştirme
    YAPILMAZ: 'Hypothetical protein' 85 satırda geçiyor, ada göre eşleştirme
    uydurma olurdu.
    """
    if not text:
        return set()
    return {f"LIC{m.group(1)}" for m in _LOCUS_RE.finditer(text)}


def tags_of(proteins) -> set[str]:
    """Protein nesnelerinden LIC locus_tag kümesi (id / description / annotations)."""
    tags: set[str] = set()
    for pr in proteins or []:
        for src in (getattr(pr, "id", None), getattr(pr, "description", None),
                    (getattr(pr, "annotations", None) or {}).get("locus_tag")):
            tags |= locus_tags(str(src) if src is not None else None)
    return tags


def match(proteins, input_tags: set[str] | None = None) -> dict:
    """Koşunun proteinlerini etiket setiyle eşleştirir (protein düzeyi).

    `proteins`: `.id` / `.description` / `.annotations` taşıyan nesneler
    (vaxforge.models.Protein ile uyumlu; ördek tiplemesi yeterli). Genellikle
    hunİden GEÇEN proteinler verilir.

    `input_tags`: girdide (eleme ÖNCESİNDE) bulunan locus_tag kümesi. Verilirse
    ikinci bir oran daha hesaplanır: "girdide bulunan etiketli proteinlerin
    kaçı elemeden geçti". Bu ayrım şart, çünkü aksi hâlde "girdide hiç yoktu"
    ile "elemede düştü" aynı düşük recall'a karışır — ikisi bambaşka kusurdur.

    Dönen sözlük her konak için AYRI recall verir ve eşleşme kurulamadıysa
    bunu açıkça söyler — sıfır eşleşme "recall %0" olarak raporlanmaz, çünkü
    çoğu zaman kimlik biçimi tutmamıştır, aday kötü değildir.
    """
    tags = tags_of(proteins)

    out: dict = {"method": METHOD, "citation": CITATION, "readout": READOUT,
                 "n_features": n_features(), "n_proteins": n_proteins(),
                 "n_run_locus_tags": len(tags), "hosts": {}}
    if not available():
        out["note"] = ("etiket seti kurulumda yok (refdata TSV eksik) — "
                       "ölçüm yapılmadı")
        return out
    if not tags:
        out["note"] = ("koşunun proteinlerinde LIC locus_tag bulunamadı — "
                       "eşleştirme YAPILAMADI (recall %0 DEĞİL: kimlik biçimi "
                       "tutmuyor; girdi Copenhageni L1-130 annotasyonu mu?)")
        return out
    for h in HOSTS:
        den = set(denominator(h))
        hit = sorted(den & tags)
        rec: dict = {
            "denominator": len(den),
            "recovered": len(hit),
            "recall_percent": (round(100 * len(hit) / len(den), 1) if den else None),
            "recovered_orfs": hit,
            "unmeasured_in_host": len(unmeasured(h)),
        }
        if input_tags is not None:
            # "Girdide hiç yoktu" ile "elemede düştü" ayrı raporlanır.
            in_input = den & set(input_tags)
            rec["in_input"] = len(in_input)
            rec["kept_of_in_input"] = len(in_input & tags)
            rec["kept_percent_of_in_input"] = (
                round(100 * len(in_input & tags) / len(in_input), 1)
                if in_input else None)
            rec["lost_in_filtering"] = sorted(in_input - tags)
        out["hosts"][h] = rec
    out["specificity_note"] = ("Bu set yalnız POZİTİF hedefleri içerir; buradan "
                              "özgüllük ÖLÇÜLEMEZ, yalnız duyarlılık (recall).")
    return out
