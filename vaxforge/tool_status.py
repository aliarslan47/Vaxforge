"""Koşu başında her dış aracın GERÇEK kullanılabilirliğini kaydeder.

NEDEN BU MODÜL VAR: rapordaki Yöntemler/Referanslar tablosu `citations.py`'den
KOŞULSUZ besleniyordu, yani kurulu olmayan ya da çöken bir araç da atıfla
birlikte "kullanıldı" gibi listeleniyordu. Kuduz koşusunda işleme aracı tam olarak
bunu yaptı: binary segfault verdi, 85 adayın tamamında `processing_norm` None
kaldı, ama rapor o aracı kendi atfıyla yöntem olarak gösterdi.

Burada toplanan durum `meta["tool_status"]` olarak run.json'a yazılır ve
`citations.for_report()` her atıf satırını `used` alanıyla işaretler. Böylece
"araç yok" bir negatif BULGU gibi görünmez.

Parametre isteyen `available()` fonksiyonları (conservation, mimicry, structure,
host_db, strain_fetch) koşuya özgüdür; burada SORGULANMAZ, çünkü yanıtları
girdiye bağlıdır ve pipeline onları kendi yerinde dürüstçe etiketler.
"""
from __future__ import annotations

# citation'daki `tool` metnine karşılık gelen modül adı -> görünen etiket.
# Yalnız parametresiz available() sunan, dışa bağımlı araçlar.
_PROBES: list[tuple[str, str, str]] = [
    # Etiket kaynak adı İÇERMEZ: kaynak VFDB ya da UniProt KW-0843 olabilir
    # (bkz. discovery.virulence_source) ve hangisi koştuysa rapordaki yöntem
    # satırı + atıf tablosu onu yazar. Sabit "VFDB" demek yanlış iddia olurdu.
    ("discovery", "Virülans DB eşleşmesi (DIAMOND)", "Virulence DB match (DIAMOND)"),
    ("psortb", "PSORTb 3.0", "PSORTb 3.0"),
    ("tmhmm_local", "TMHMM-2.0", "TMHMM-2.0"),
    ("signalp", "SignalP-5.0", "SignalP-5.0"),
    ("iapred", "IApred", "IApred"),
    ("netmhc_local", "NetMHCpan / NetMHCIIpan", "NetMHCpan / NetMHCIIpan"),
    ("iedb_processing", "IEDB processing API (kesim+TAP)",
     "IEDB processing API (cleavage+TAP)"),
    ("bepipred", "BepiPred-1.0", "BepiPred-1.0"),
    ("bepipred3", "BepiPred-3.0", "BepiPred-3.0"),
    ("allergen", "Alerjen referans seti", "Allergen reference set"),
    ("toxinpred", "ToxinPred2", "ToxinPred2"),
    ("ifngamma", "IFNepitope2", "IFNepitope2"),
    ("vaxijen", "VaxiJen ACC", "VaxiJen ACC"),
    ("antigen_acc", "Antijenite ACC modeli", "Antigenicity ACC model"),
    ("population", "IEDB Population Coverage", "IEDB Population Coverage"),
    ("rsa", "RSA / Shrake-Rupley", "RSA / Shrake-Rupley"),
    ("secstruct", "İkincil yapı", "Secondary structure"),
    ("solubility", "Çözünürlük", "Solubility"),
    ("disorder", "metapredict (düzensizlik)", "metapredict (disorder)"),
    ("hla_disease", "HLA-hastalık anotasyonu", "HLA-disease annotation"),
    ("immforge_score", "ImmForge immünizasyon", "ImmForge immunization"),
]


def collect() -> dict[str, dict]:
    """{modül: {label, available, note}} — her araç için tek gerçek kontrol.

    Bir araç kendi `available()` çağrısında istisna atarsa bu YUTULMAZ; durum
    False + istisna tipi olarak kaydedilir (sessiz True'dan iyidir).
    """
    import importlib

    out: dict[str, dict] = {}
    for mod_name, label_tr, label_en in _PROBES:
        # Etiket İKİ DİLLİ tutulur: EN raporda Türkçe araç adı görünmesin
        # (proje kuralı: arayüz/rapor metinleri daima TR + EN).
        entry = {"label": label_tr, "label_en": label_en,
                 "available": False, "check": "shallow", "note": ""}
        try:
            mod = importlib.import_module(f"vaxforge.{mod_name}")
            # DERİN kontrol: modül `verify()` sunuyorsa aracı GERÇEKTEN
            # çalıştırır. Gerekçe: `available()` çoğu araçta yalnız dosya
            # varlığına bakıyor ve bu YALAN söyleyebiliyor — TMHMM (bozuk
            # shebang), IApred (eksik bağımlılık) ve ToxinPred2 (eski pickle
            # + pandas API) hepsi "hazır" diyip hiç çalışmıyordu. `check`
            # alanı raporun hangi kontrolün yapıldığını söylemesini sağlar;
            # "shallow" bir GARANTİ DEĞİLDİR, sadece dosya var demektir.
            if hasattr(mod, "verify"):
                entry["available"] = bool(mod.verify())
                entry["check"] = "verified"
            else:
                entry["available"] = bool(mod.available())
        except Exception as e:                       # noqa: BLE001
            entry["note"] = f"kontrol hatası: {type(e).__name__}"
        out[mod_name] = entry
    return out


def _lbl(entry: dict, lang: str) -> str:
    return entry.get("label_en", entry["label"]) if lang == "en" else entry["label"]


def available_labels(status: dict[str, dict] | None, lang: str = "tr") -> set[str]:
    """Kullanılabilir araçların görünen etiketleri (citations eşlemesi için)."""
    if not status:
        return set()
    return {_lbl(v, lang) for v in status.values() if v.get("available")}


def missing(status: dict[str, dict] | None, lang: str = "tr") -> list[str]:
    """Kullanılamayan araçların etiketleri (rapor uyarısı için, sıralı)."""
    if not status:
        return []
    return sorted(_lbl(v, lang) for v in status.values() if not v.get("available"))


def summarize(status: dict[str, dict] | None) -> dict:
    """Özet; eksik araç etiketleri İKİ DİLDE ayrı listelenir."""
    if not status:
        return {"total": 0, "available": 0, "missing": 0,
                "missing_labels": [], "missing_labels_en": []}
    miss = missing(status, "tr")
    verified = sum(1 for v in status.values() if v.get("check") == "verified")
    return {"total": len(status), "available": len(status) - len(miss),
            "missing": len(miss), "missing_labels": miss,
            "missing_labels_en": missing(status, "en"),
            "verified": verified, "shallow": len(status) - verified}
