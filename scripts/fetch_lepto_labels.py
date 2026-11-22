#!/usr/bin/env python3
"""Leptospira etiket kaynaklarını indirir — DÜRÜST durum raporuyla.

NE YAPAR
--------
Leptospira için *in silico* immün öngörüsünü eğitmek/değerlendirmek üzere
kullanılabilecek NİCEL etiket kaynaklarını (tüm-proteom protein mikroarray,
ORFeome fay gösterimi, immünoproteomik, sığır IFN-γ paneli) kaynağından çekmeye
çalışır ve her kaynak için SONUCU AÇIKÇA bildirir:

    indirildi      -> dosya diskte, boyutu yazılır
    erisilemez     -> yayıncı açık erişim vermiyor / 403 / ödeme duvarı
    ag_yok         -> internet yok
    bulunamadi     -> URL 404

NE YAPMAZ
---------
Veri UYDURMAZ. Bir ek dosya çekilemediyse o kaynak 'erisilemez' olarak raporlanır
ve bu repoda hiçbir yere sahte sayı yazılmaz. Bu, repodaki araç-dürüstlüğü
desenin gereğidir: bir yeteneğin "var" görünüp iş anında sessizce düşmesi
yasaktır (bkz. [[sessiz-arac-dususu-deseni]]).

Ayrıca: indirilen ek dosyalar GÜVENİLMEYEN VERİDİR. Her biri kendi boş
dizinine yazılır; bu betik onların içinden hiçbir şey ÇALIŞTIRMAZ ve içeriğini
yorumlamaz, yalnız indirir ve özetler.

ETİKETLERİN ANLAMI (önemli)
---------------------------
Bu kaynakların hepsi TANIMA (antikor bağlanması ya da IFN-γ geri çağrımı)
ölçer, KORUMA ölçmez. Leptospira'da doğrulanmış koruma korelatı anti-LPS
aglütine eden IgG'dir ve protein tabanlı bir pipeline onu modellemez
(bkz. vaxforge/correlate.py). Bu yüzden indirilen her kaynağın 'readout' alanı
raporda taşınır.

KULLANIM
--------
    python3 scripts/fetch_lepto_labels.py                 # indirmeyi dener
    python3 scripts/fetch_lepto_labels.py --dry-run       # yalnız planı yazar
    python3 scripts/fetch_lepto_labels.py --out DIZIN
"""

from __future__ import annotations

import argparse
import json
import re
import socket
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUT = REPO / "data" / "labels" / "leptospira"
UA = "Vaxforge-label-fetcher/1.0 (bioinformatics pipeline; contact via repository)"
TIMEOUT = 60


@dataclass
class Wanted:
    """İndirilmek istenen tek bir ek dosya."""

    name: str                 # yayıncıdaki dosya adı (eutils'ten çözülür)
    why: str                  # bu dosya NEDEN etiket kaynağı
    status: str = "denenmedi"
    detail: str = ""
    path: str | None = None
    size: int | None = None
    caption: str = ""
    tried: list[str] = field(default_factory=list)


@dataclass
class Source:
    """Tek bir etiket kaynağı (bir makale ve onun ek dosyaları)."""

    key: str
    citation: str
    pmid: str
    pmcid: str                # ek dosya ADLARINI çözmek için (eutils)
    doi: str
    publisher: str            # plos | elife | acs | elsevier
    readout: str              # ETİKET NE ÖLÇÜYOR
    scale: str                # kapsam (kaç protein/peptit)
    want: list[Wanted]        # istenen ek dosyalar (ad kalıbıyla eşleştirilir)
    note: str = ""
    status: str = "denenmedi"
    detail: str = ""
    available: list[str] = field(default_factory=list)  # eutils'in bildirdiği


# --------------------------------------------------------------------------- #
# Kaynaklar. URL'ler PMC açık erişim ek dosya yollarıdır; yayıncıda kalanlar
# için erişilemezlik BEKLENEN sonuçtur ve öyle raporlanır.
# --------------------------------------------------------------------------- #
def sources() -> list[Source]:
    """Etiket kaynakları. Dosya ADLARI eutils'ten CANLI çözülür, elle yazılmaz."""
    return [
        Source(
            key="wunder_2021_immunity_correlates",
            citation=("Wunder EA Jr, ve ark. A live attenuated-vaccine model "
                      "confers cross-protective immunity against different "
                      "species of the Leptospira genus. eLife. 2021;10:e64166."),
            pmid="33522488", pmcid="PMC7837694", doi="10.7554/eLife.64166",
            publisher="elife",
            readout=("ÇAPRAZ-KORUYUCU bağışıklıkla ilişkili proteome array "
                     "hedefleri (hamster + fare); korumaya en yakın etiket"),
            scale="154 protein hedefi; ayrıca 41 hedefin tür-arası kimliği",
            want=[
                Wanted("supp4", "bağışıklık korelatı olarak belirlenen 154 "
                                "protein hedefinin TAM listesi — ana etiket"),
                Wanted("fig4-data1", "proteome array HAM verisi (nicel sinyal)"),
                Wanted("fig5-data1", "41 seçili hedefin Leptospira türleri "
                                     "arasındaki amino asit kimliği (korunmuşluk)"),
                Wanted("fig3-data1", "MAT ve pasif transfer ham verisi "
                                     "(titre-koruma ilişkisini sınamak için)"),
            ],
            note=("EN DEĞERLİ KAYNAK: etiket 'tanıma' değil, çapraz-KORUMA ile "
                  "ilişkilendirilmiş repertuvar. Yine de koruma korelatı "
                  "DEĞİL, korelasyondur."),
        ),
        Source(
            key="lessa_aquino_2017_severity",
            citation=("Lessa-Aquino C, ve ark. Distinct antibody responses of "
                      "patients with mild and severe leptospirosis determined "
                      "by whole proteome microarray analysis. PLoS Negl Trop "
                      "Dis. 2017;11(1):e0005349."),
            pmid="28141801", pmcid="PMC5302828", doi="10.1371/journal.pntd.0005349",
            publisher="plos",
            readout="insan hasta serumu IgG/IgM reaktivitesi; HAFİF vs ŞİDDETLİ",
            scale="proteom ölçekli mikroarray; akut + konvalesan faz",
            want=[
                Wanted("s005", "akut faz serotanısal antijenler (hafif/şiddetli)"),
                Wanted("s006", "konvalesan faz serotanısal antijenler"),
                Wanted("s004", "mikroarray'de TEMSİL EDİLMEYEN proteinler — "
                               "paydayı bilmek için zorunlu"),
            ],
            note=("Şiddet katmanı, 'yüksek antikor = iyi' varsayımını sınamak "
                  "için kritik: şiddetli hastalarda reaktivite daha yüksek "
                  "olabilir."),
        ),
        Source(
            key="lessa_aquino_2013_seroreactive",
            citation=("Lessa-Aquino C, ve ark. Identification of seroreactive "
                      "proteins of Leptospira interrogans serovar Copenhageni "
                      "using a high-density protein microarray approach. "
                      "PLoS Negl Trop Dis. 2013;7(10):e2499."),
            pmid="24132137", pmcid="PMC3798601", doi="10.1371/journal.pntd.0002499",
            publisher="plos",
            readout="insan hasta serumu IgG/IgM reaktivitesi (protein mikroarray)",
            scale="~2700 protein ifade edildi; reaktif alt küme",
            want=[
                Wanted("s005", "reaktif antijen listesi / ek tablo"),
                Wanted("s006", "ek tablo"),
                Wanted("s007", "ek tablo"),
                Wanted("s008", "ek tablo"),
                Wanted("s009", "ek tablo"),
            ],
            note="Tam açık erişim (PLoS).",
        ),
        Source(
            key="lessa_aquino_2015_predictive",
            citation=("Lessa-Aquino C, ve ark. Proteomic features predict "
                      "seroreactivity against leptospiral antigens in "
                      "leptospirosis patients. J Proteome Res. 2015;14(2):549-56."),
            pmid="25358092", pmcid="PMC4286151", doi="10.1021/pr500718t",
            publisher="acs",
            readout="serore­aktiviteyi ÖNGÖREN proteomik özellikler",
            scale="3359 ORF = proteomun ~%91'i; 191 reaktif antijen (~%5)",
            want=[Wanted("si_001", "destekleyici bilgi (PDF)")],
            note=("ACS yayıncısı otomatik indirmeyi engelliyor (403 bekleniyor) "
                  "-> 'erisilemez' raporlanır, veri UYDURULMAZ. Sayılar "
                  "makalenin kendisinden okunmalıdır."),
        ),
        Source(
            key="deveson_lucas_2014_bovine_ifng",
            citation=("Deveson Lucas DS, ve ark. Recombinant Leptospira "
                      "borgpetersenii serovar Hardjo proteins elicit "
                      "interferon-gamma from cattle vaccinated against "
                      "leptospirosis. Vet Microbiol. 2014;169(3-4):163-70."),
            pmid="24418243", pmcid="", doi="10.1016/j.vetmic.2013.12.016",
            publisher="elsevier",
            readout="aşılanmış SIĞIR hücrelerinden IFN-γ geri çağrımı (Th1)",
            scale="238 rekombinant protein; PROTEİN düzeyi (peptit değil)",
            want=[Wanted("—", "ek tablo (açık erişim değil)")],
            note=("TEK protein-ölçekli T-hücre etiketi ve alanın en iyi "
                  "doğrulanmış protein korelatına (sığır Hardjo Th1) bağlı. "
                  "PMC'de YOK, Elsevier kapalı -> 'erisilemez'. Bu kaynak "
                  "gerekiyorsa kurumsal erişimle elle indirilmeli."),
        ),
    ]


def publisher_urls(src: Source, fname: str) -> list[str]:
    """Yayıncıya göre denenecek indirme adresleri (sırayla)."""
    if src.publisher == "plos":
        journal = "plosntds" if "pntd" in fname else "plosone"
        stem = fname.rsplit(".", 1)[0]          # pntd.0005349.s005
        sid = stem.replace("pntd.", "journal.pntd.").replace("pone.", "journal.pone.")
        return [f"https://journals.plos.org/{journal}/article/file"
                f"?id=10.1371/{sid}&type=supplementary"]
    if src.publisher == "elife":
        base, ext = fname.rsplit(".", 1)
        art = src.doi.rsplit(".", 1)[-1]        # 64166
        # eLife CDN sürüm eki ister; v1..v3 denenir.
        return [f"https://cdn.elifesciences.org/articles/{art}/{base}-v{v}.{ext}"
                for v in (1, 2, 3)]
    if src.publisher == "acs":
        return [f"https://pubs.acs.org/doi/suppl/{src.doi}/suppl_file/{fname}"]
    return [f"https://doi.org/{src.doi}"]


def resolve_filenames(src: Source) -> tuple[list[tuple[str, str]], str]:
    """Ek dosya adlarını NCBI eutils'ten çözer (elle yazılmaz).

    Dönen: ([(dosya_adi, caption), ...], hata_mesaji)
    """
    if not src.pmcid:
        return [], "PMC kaydı yok (açık erişim değil)"
    pid = src.pmcid.replace("PMC", "")
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=pmc&id={pid}")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            xml = r.read().decode("utf-8", errors="replace")
    except Exception as e:                                    # noqa: BLE001
        return [], f"eutils erişilemedi: {e}"
    out: list[tuple[str, str]] = []
    for m in re.finditer(
            r"<supplementary-material[^>]*>(.{0,1500}?)</supplementary-material>",
            xml, re.S):
        blk = m.group(1)
        h = re.search(r'xlink:href="([^"]+)"', blk)
        if not h:
            continue
        cap = ""
        c = re.search(r"<caption>(.*?)</caption>", blk, re.S)
        if c:
            cap = " ".join(re.sub(r"<[^>]+>", " ", c.group(1)).split())
        out.append((h.group(1), cap))
    if not out:
        return [], "eutils kaydında ek dosya bulunamadı"
    return out, ""


# --------------------------------------------------------------------------- #
def online(host: str = "pmc.ncbi.nlm.nih.gov", port: int = 443) -> bool:
    try:
        socket.create_connection((host, port), timeout=10).close()
        return True
    except OSError:
        return False


def _download(urls: list[str], dest_dir: Path, w: Wanted) -> None:
    """Adresleri sırayla dener; sonucu w üzerine DÜRÜSTÇE yazar."""
    for url in urls:
        w.tried.append(url)
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                data = r.read()
                ctype = r.headers.get("Content-Type", "")
        except urllib.error.HTTPError as e:
            w.status = "bulunamadi" if e.code == 404 else "erisilemez"
            w.detail = f"HTTP {e.code}"
            continue
        except urllib.error.URLError as e:
            w.status, w.detail = "ag_yok", str(e.reason)
            continue
        except (TimeoutError, socket.timeout):
            w.status, w.detail = "erisilemez", "zaman aşımı"
            continue

        if not data:
            w.status, w.detail = "erisilemez", "boş yanıt"
            continue
        # Veri yerine HTML = ödeme duvarı / yönlendirme / hata sayfası.
        head = data[:400].lstrip().lower()
        if head.startswith((b"<!doctype html", b"<html")):
            w.status = "erisilemez"
            w.detail = (f"veri yerine HTML döndü ({len(data)} bayt) — ödeme "
                        f"duvarı/bot engeli olasılığı")
            continue
        if head.startswith(b"<?xml") and len(data) < 2000:
            w.status, w.detail = "bulunamadi", f"XML hata yanıtı ({len(data)} bayt)"
            continue

        path = dest_dir / url.rsplit("/", 1)[-1].split("?")[0]
        if not path.suffix:
            path = dest_dir / w.name
        path.write_bytes(data)
        w.status, w.path, w.size = "indirildi", str(path), len(data)
        w.detail = f"{len(data)} bayt, {ctype or '?'}"
        return

    if w.status == "denenmedi":
        w.status, w.detail = "erisilemez", "hiçbir adres sonuç vermedi"


def fetch(src: Source, out_dir: Path) -> None:
    """Kaynağın istenen ek dosyalarını indirir; her sonucu kaydeder."""
    # Her kaynak KENDİ dizinine: indirilen dosyalar güvenilmeyen veridir.
    dest_dir = out_dir / src.key
    dest_dir.mkdir(parents=True, exist_ok=True)

    found, err = resolve_filenames(src)
    src.available = [n for n, _ in found]
    if err:
        src.status, src.detail = "erisilemez", err
        for w in src.want:
            w.status, w.detail = "erisilemez", err
        return

    for w in src.want:
        # eutils'in bildirdiği gerçek dosya adını kalıpla eşleştir.
        hit = next(((n, c) for n, c in found if w.name in n), None)
        if hit is None:
            w.status = "bulunamadi"
            w.detail = (f"'{w.name}' kalıbı makalenin ek dosyaları arasında yok "
                        f"(mevcut: {', '.join(src.available[:6])})")
            continue
        fname, w.caption = hit
        _download(publisher_urls(src, fname), dest_dir, w)

    got = [w for w in src.want if w.status == "indirildi"]
    src.status = ("indirildi" if len(got) == len(src.want)
                  else "kismi" if got else "erisilemez")
    src.detail = f"{len(got)}/{len(src.want)} dosya"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help=f"çıktı dizini (varsayılan {DEFAULT_OUT})")
    ap.add_argument("--dry-run", action="store_true",
                    help="indirme yapmaz, yalnız planı yazar")
    args = ap.parse_args()

    srcs = sources()

    print("=" * 78)
    print("Leptospira ETİKET KAYNAKLARI")
    print("=" * 78)
    print("UYARI: Bu kaynakların çoğu TANIMA ölçer, KORUMA ölçmez.")
    print("Leptospira'nın doğrulanmış korelatı anti-LPS aglütine eden IgG'dir ve")
    print("protein tabanlı bir pipeline onu modellemez (vaxforge/correlate.py).")
    print("Tek istisna Wunder 2021: etiketi çapraz-korumayla İLİŞKİLENDİRİLMİŞ")
    print("repertuvardır — ama bu da korelasyondur, korelat değildir.")
    print()

    if args.dry_run:
        for s in srcs:
            print(f"[plan] {s.key}  ({s.publisher}, {s.pmcid or 'PMC yok'})")
            print(f"       okuma   : {s.readout}")
            print(f"       kapsam  : {s.scale}")
            print(f"       atıf    : {s.citation}")
            print(f"                 doi:{s.doi} (PMID {s.pmid})")
            for w in s.want:
                print(f"       istenen : {w.name} — {w.why}")
            if s.note:
                print(f"       not     : {s.note}")
            print()
        print(f"{len(srcs)} kaynak planlandı; --dry-run olduğu için İNDİRME YAPILMADI.")
        return 0

    if not online():
        for s in srcs:
            s.status, s.detail = "ag_yok", "ağ erişimi yok (bağlantı testi başarısız)"
            for w in s.want:
                w.status, w.detail = "ag_yok", "ağ erişimi yok"
        print("AĞ YOK: hiçbir kaynak indirilemedi. Hiçbir veri uydurulmadı.\n")
    else:
        args.out.mkdir(parents=True, exist_ok=True)
        for s in srcs:
            print(f"[{s.key}]", flush=True)
            fetch(s, args.out)
            for w in s.want:
                mark = {"indirildi": "+", "erisilemez": "-",
                        "bulunamadi": "?", "ag_yok": "x"}.get(w.status, "·")
                print(f"   {mark} {w.name:12s} {w.status:12s} {w.detail}")
            print(f"   => {s.status}: {s.detail}")

    # Özet
    print()
    print("=" * 78)
    print("ÖZET")
    print("=" * 78)
    files_ok = 0
    for s in srcs:
        got = [w for w in s.want if w.status == "indirildi"]
        files_ok += len(got)
        mark = {"indirildi": "+", "kismi": "~", "erisilemez": "-",
                "ag_yok": "x"}.get(s.status, "·")
        nbytes = sum(w.size or 0 for w in got)
        extra = f"  ({nbytes} bayt)" if nbytes else ""
        print(f" {mark} {s.key:34s} {len(got)}/{len(s.want)} dosya{extra}")

    print()
    if not files_ok:
        print("  HİÇBİR dosya indirilemedi. Etiket katmanı KURULAMADI.")
        print("  Bu bir başarısızlıktır ve öyle raporlanıyor; sahte etiket üretilmedi.")
    else:
        print(f"  {files_ok} dosya diskte: {args.out}")
        print("  İndirilen dosyalar GÜVENİLMEYEN VERİDİR: her biri kendi")
        print("  dizininde, çalıştırılmadı, içeriği bu betik tarafından")
        print("  yorumlanmadı — yalnız indirildi ve boyutu doğrulandı.")
    blocked = [s.key for s in srcs if s.status == "erisilemez"]
    if blocked:
        print(f"\n  ERİŞİLEMEYEN kaynaklar ({len(blocked)}): {', '.join(blocked)}")
        print("  Bunlar için veri ÜRETİLMEDİ. Gerekiyorsa kurumsal erişimle")
        print("  elle indirilip aynı dizine konmalıdır.")

    # Makine-okunur durum kaydı
    manifest = {
        "generated": date.today().isoformat(),
        "out_dir": str(args.out),
        "readout_warning": (
            "Kaynakların çoğu TANIMA ölçer, KORUMA ölçmez; Leptospira'nın "
            "doğrulanmış korelatı anti-LPS aglütine eden IgG'dir."),
        "n_files_downloaded": files_ok,
        "sources": [
            {"key": s.key, "status": s.status, "detail": s.detail,
             "pmid": s.pmid, "pmcid": s.pmcid, "doi": s.doi,
             "publisher": s.publisher, "readout": s.readout, "scale": s.scale,
             "citation": s.citation, "note": s.note,
             "available_in_pmc": s.available,
             "files": [{"want": w.name, "why": w.why, "status": w.status,
                        "detail": w.detail, "path": w.path, "size": w.size,
                        "caption": w.caption, "urls_tried": w.tried}
                       for w in s.want]}
            for s in srcs
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    mf = args.out / "manifest.json"
    mf.write_text(json.dumps(manifest, ensure_ascii=False, indent=2),
                  encoding="utf-8")
    print(f"\n  durum kaydı: {mf}")

    return 0 if files_ok else 1


if __name__ == "__main__":
    sys.exit(main())
