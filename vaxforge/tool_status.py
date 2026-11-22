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

import time

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


#: `verify()` düşerse bu kadar bekleyip BİR kez daha denenir (geçici yük hatası).
_RETRY_DELAY = 2.0


def _verify_with_retry(mod) -> tuple[bool, str]:
    """`mod.verify()` — düşerse bir kez daha dener. Dönüş: (ok, son hata adı).

    Geçici kaynak darlığı (paralel araç koşuları, bellek baskısı) verify'ı
    düşürebiliyor; tek deneme bunu kalıcı "ARAÇ YOK" diye rapora yazıyordu.
    """
    last = ""
    for attempt in (1, 2):
        try:
            if bool(mod.verify()):
                return True, ""
            last = "verify() False döndü"
        except Exception as e:                       # noqa: BLE001
            last = type(e).__name__
        if attempt == 1:
            time.sleep(_RETRY_DELAY)
    return False, last


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
                # YÜK ALTINDA YANLIŞ "ARAÇ YOK" (2026-11-22): `verify()` aracı
                # gerçekten koşturduğu için GEÇİCİ bir hatayla da düşebilir.
                # Ölçüldü: 9 paralel netMHCIIpan + pytest koşarken ToxinPred2'nin
                # verify'ı düştü ve araç rapora kalıcı "eksik" yazıldı; tek başına
                # `verify() -> True`. Dürüstlük tek yönlü değil — "var" demek
                # kadar "yok" demek de yanlış olabilir, ve eksik sanılan araç
                # ÖLÇÜLMEMİŞ sonucu olumsuz bulgu gibi gösterir.
                # Bu yüzden: bir kez RETRY, ve "dosya yok" ile "verify düştü"
                # AYRI kaydedilir (ikincisi geçici olabilir, raporda öyle görünür).
                ok, err = _verify_with_retry(mod)
                entry["available"] = ok
                entry["check"] = "verified"
                if not ok:
                    # Dosya hiç yok mu, yoksa dosya var ama koşmadı mı?
                    shallow = None
                    if hasattr(mod, "available"):
                        try:
                            shallow = bool(mod.available())
                        except Exception:            # noqa: BLE001
                            shallow = None
                    if shallow is False:
                        entry["note"] = "dosya yok (kurulu değil)"
                        entry["fail_kind"] = "missing"
                    else:
                        entry["note"] = (
                            f"dosya var ama doğrulama koşusu düştü"
                            + (f" ({err})" if err else "")
                            + " — geçici yük hatası olabilir")
                        entry["fail_kind"] = "verify_failed"
                        entry["retried"] = True
            else:
                entry["available"] = bool(mod.available())
                if not entry["available"]:
                    entry["fail_kind"] = "missing"
        except Exception as e:                       # noqa: BLE001
            entry["note"] = f"kontrol hatası: {type(e).__name__}"
            entry["fail_kind"] = "probe_error"
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
    # EKSİKLİĞİN TÜRÜ AYRI RAPORLANIR: "kurulu değil" kalıcı bir olgu,
    # "doğrulama koşusu düştü" ise geçici yük hatası olabilir (bir kez retry
    # edildi ve yine düştü). İkisini aynı listede göstermek, çalışan bir aracı
    # kalıcı yokmuş gibi raporluyordu.
    def _k(kind, lang):
        return sorted(_lbl(v, lang) for v in status.values()
                      if not v.get("available") and v.get("fail_kind") == kind)
    return {"total": len(status), "available": len(status) - len(miss),
            "missing": len(miss), "missing_labels": miss,
            "missing_labels_en": missing(status, "en"),
            "verified": verified, "shallow": len(status) - verified,
            "not_installed": _k("missing", "tr"),
            "not_installed_en": _k("missing", "en"),
            "verify_failed": _k("verify_failed", "tr"),
            "verify_failed_en": _k("verify_failed", "en"),
            "probe_error": _k("probe_error", "tr"),
            "probe_error_en": _k("probe_error", "en")}
