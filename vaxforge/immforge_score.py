"""VaxForge → ImmForge glue: MEV/aday konstruktu → in silico immünizasyon → OLUR/OLMAZ.

VaxForge pipeline'ının SON adımı bunu çağırır: üretilen aşı konstruktunu (MEVConstruct.seq)
seçili konak(lar)ın MHC allelleri için GERÇEK NetMHCpan/NetMHCIIpan'den geçirip ImmForge
motorunda immünize eder → immünojenisite skoru + karar. ImmForge = adayları eleyen immün test.

ORTAM-BAĞIMSIZ: NetMHCpan/BepiPred nerede kuruluysa (yerel VEYA servis sunucusu) orada çalışır;
predict_fn enjekte edilir. immforge paketi IMMFORGE_PATH ile, o yoksa bilinen yerel
konumlar (~/immforge, ~/repos/ImmForge, repo komşusu ImmForge) taranarak bulunur
(serviste pip-kurulu olur — netmhc_local'ın path deseniyle aynı).

DÜRÜSTLÜK: NetMHCpan yoksa uydurma yapmaz → {"verdict": "ARAÇ_YOK"}. Tür-özgü ImmForge
kalibrasyonu olmayan konak (ör. chicken) → human motoru proxy + açık not.
"""
from __future__ import annotations

import functools
import importlib
import os
import sys
from pathlib import Path

from vaxforge import bepipred, netmhc_local

def _resolve_immforge_path() -> str:
    """ImmForge kurulumunun yolunu bulur.

    IMMFORGE_PATH verilmişse SORGUSUZ o kullanılır (servis/konteyner dağıtımı).
    Verilmemişse bilinen yerel konumlar sırayla denenir; hiçbiri yoksa ilk aday
    döner ve `available()` dürüstçe False olur.

    Tarihsel not: eski varsayılan yalnız ~/immforge idi; repo ~/repos/ImmForge
    altında klonlanmışsa kalibrasyon katmanı SESSİZCE devre dışı kalıyordu
    (tüm konaklar human-proxy). Aday listesi bu sessiz düşüşü kapatır.

    İkinci sessiz düşüş (2026-11-22): aday testi `config/` ya da `immforge/`
    dizinine bakıyordu, ama entegrasyon giriş noktası
    (`immforge/integrations/vaxforge.py`) YALNIZ phase-4 dalında var. ~/repos/ImmForge
    `main`'de duruyorsa test geçiyor, sonra import başarısız oluyor ve ImmForge adımı
    IMMFORGE_PATH elle verilmedikçe sessizce boş kalıyordu — `web/guardian/supervise.sh`
    onu kendisi set ettiği için web arayüzünde çalışıyor, CLI/script koşusunda
    çalışmıyordu. Artık aday, giriş noktasını GERÇEKTEN taşıyorsa seçilir; worktree
    kardeşleri (ImmForge-phase4) de taranır.
    """
    env = os.environ.get("IMMFORGE_PATH")
    if env:
        return env
    home = Path.home()
    repo_parent = Path(__file__).resolve().parents[2]
    candidates = [
        home / "immforge",
        home / "repos" / "ImmForge-phase4",
        home / "repos" / "ImmForge",
        home / "ImmForge-phase4",
        home / "ImmForge",
        repo_parent / "ImmForge-phase4",
        repo_parent / "ImmForge",
        repo_parent / "immforge",
    ]
    entry_rel = Path("immforge") / "integrations" / "vaxforge.py"
    # 1. tur: giriş noktasını gerçekten taşıyan kurulum (iş yapabilen).
    for c in candidates:
        if (c / entry_rel).is_file():
            return str(c)
    # 2. tur: ImmForge gibi duran ama entegrasyonu olmayan kurulum. Seçilir ki
    # available() "import edilemedi" diye DÜRÜST bir yokluk raporlasın.
    for c in candidates:
        if (c / "config").is_dir() or (c / "immforge").is_dir():
            return str(c)
    return str(candidates[0])


IMMFORGE_PATH = _resolve_immforge_path()

def _calib_file(species: str) -> Path:
    return Path(IMMFORGE_PATH) / "config" / f"calibration.{species}.yaml"


def calibrated_species() -> set[str]:
    """GERÇEKTE kalibrasyon YAML'ı olan türler (dosya-tabanlı — hardcoded liste DEĞİL).

    Dürüstlük: bir tür yalnız calibration.{tür}.yaml gerçekten varsa 'kalibre' sayılır;
    dosya yoksa o konak human motoruyla proxy'lenir (uydurma yok).
    """
    d = Path(IMMFORGE_PATH) / "config"
    if not d.is_dir():
        return set()
    return {p.stem.split(".", 1)[1] for p in d.glob("calibration.*.yaml")
            if "." in p.stem}


def _ensure_immforge() -> None:
    if IMMFORGE_PATH not in sys.path and Path(IMMFORGE_PATH).is_dir():
        sys.path.insert(0, IMMFORGE_PATH)


ENTRY_MODULE = "immforge.integrations.vaxforge"
ENTRY_FUNC = "evaluate_candidate"


@functools.lru_cache(maxsize=1)
def entrypoint_available() -> tuple[bool, str]:
    """Entegrasyon giriş noktası GERÇEKTEN import edilebiliyor mu? (neden, ile birlikte)

    NEDEN GEREKLİ: `available()` eskiden yalnız ImmForge DİZİNİNİN varlığına ve
    NetMHCpan'a bakıyordu. Dizin varsa "araç var" diyor, pipeline da adımı
    çalıştırmayı deniyor ve `immforge.integrations.vaxforge` paketi mevcut
    olmadığı için ModuleNotFoundError ile düşüyordu. Sonuç: dürüst bir
    ARAÇ_YOK yerine BEKLENMEYEN bir istisna — projenin "araç yoksa dürüstçe
    söyle" ilkesine aykırı.

    DÜZELTME (2026-11-13): Bu alt paket ImmForge'da YAZILMIŞ ve eksiksiz —
    `immforge/integrations/vaxforge.py`, `origin/phase-4.13b-anchors` dalında
    (commit d916af8 + 3 takip commit'i); `evaluate_candidate` ve
    `bepipred_residues_to_bepitopes` imzaları buradaki çağrılarla birebir uyuyor.
    Bu makinede görünmemesinin sebebi kod eksiği DEĞİL: PC formatı sonrası
    ImmForge shallow (tek commit) ve yalnız `main` klonlandı, o dal hiç inmedi.
    Çözüm kod yazmak değil, dalı getirmek:
        git -C <ImmForge> fetch --unshallow origin '+refs/heads/*:refs/remotes/origin/*'
        git -C <ImmForge> worktree add ../ImmForge-phase4 origin/phase-4.13b-anchors
        IMMFORGE_PATH=<...>/ImmForge-phase4
    Bu kontrol yine de gerekli: dal yoksa çökme yerine etiketli bir yokluk verir.
    """
    if not Path(IMMFORGE_PATH).is_dir():
        return False, f"ImmForge dizini yok: {IMMFORGE_PATH}"
    _ensure_immforge()
    try:
        mod = importlib.import_module(ENTRY_MODULE)
    except Exception as e:
        return False, f"{ENTRY_MODULE} import edilemedi ({type(e).__name__})"
    if not hasattr(mod, ENTRY_FUNC):
        return False, f"{ENTRY_MODULE}.{ENTRY_FUNC} yok"
    return True, "hazır"


def available() -> bool:
    """ImmForge skorlaması bu ortamda çalışabilir mi?

    Üç koşul: ImmForge dizini + entegrasyon giriş noktası + NetMHCpan.
    """
    ok, _ = entrypoint_available()
    return ok and netmhc_local.available()


def _immforge_species(host_name: str) -> str:
    """Konağın kendi kalibrasyonu varsa onu, yoksa 'human' (proxy) kullan."""
    return host_name if _calib_file(host_name).exists() else "human"


def _provenance(species: str) -> str:
    """Kalibrasyon YAML meta.provenance'ından kısa köken bilgisi (varsa)."""
    cfg = _calib_file(species)
    if not cfg.exists():
        return ""
    try:
        import yaml
        meta = (yaml.safe_load(cfg.read_text()) or {}).get("meta", {})
        prov = str(meta.get("provenance", "")).replace("\n", " ").strip()
        return (prov[:180] + "…") if len(prov) > 180 else prov
    except Exception:
        return ""


def _load_params(species: str):
    _ensure_immforge()
    from immforge.params import load_params
    cfg = _calib_file(species)
    return load_params(calibration_path=str(cfg)) if cfg.exists() else load_params()


def _bepitopes(seq: str) -> list[tuple[str, float]]:
    """BepiPred artık-skorları → lineer B-epitop peptitleri (varsa)."""
    if not bepipred.available():
        return []
    _ensure_immforge()
    from immforge.integrations.vaxforge import bepipred_residues_to_bepitopes
    res = bepipred.predict_residues(seq)
    return bepipred_residues_to_bepitopes(seq, res) if res else []


def score_construct(seq: str, hosts, *, name: str = "MEV") -> dict:
    """Aday dizisi + konaklar → her konak için OLUR/OLMAZ kararı + genel özet.

    hosts: vaxforge.hosts.Host listesi (.name, .label, .alleles("mhc_i"/"mhc_ii")).
    Dönüş (JSON): {tool, bepipred, headline, per_host:[...], note}.
    """
    if not seq or not hosts:
        return {"tool": "none", "verdict": "ARAÇ_YOK",
                "note": "Boş konstrukt ya da konak yok."}
    if not available():
        ok, why = entrypoint_available()
        reason = why if not ok else "NetMHCpan bu ortamda çalıştırılamıyor"
        return {"tool": "none", "verdict": "ARAÇ_YOK",
                "note": f"İmmünizasyon testi yapılamadı: {reason}.",
                "reason": reason}

    _ensure_immforge()
    from immforge.integrations.vaxforge import evaluate_candidate

    b_eps = _bepitopes(seq)
    per_host: list[dict] = []
    for h in hosts:
        species = _immforge_species(h.name)
        params = _load_params(species)
        v = evaluate_candidate(
            name, seq,
            mhci_alleles=list(h.alleles("mhc_i")),
            mhcii_alleles=list(h.alleles("mhc_ii")),
            predict_fn=netmhc_local.predict,
            species=species, params=params,
            b_epitopes=b_eps or None,
        )
        d = v.to_dict()
        d["host"] = h.name
        d["host_label"] = getattr(h, "label", h.name)
        d["species_calibrated"] = (species == h.name)
        d["calibration_species"] = species   # gerçekte kullanılan kalibrasyon türü
        d["calibration_provenance"] = _provenance(species)
        if not d["species_calibrated"]:
            d["reasons"] = list(d.get("reasons", [])) + [
                f"'{h.name}' için ImmForge tür-özgü kalibrasyonu YOK → '{species}' motoru "
                "proxy (dürüst; türe-özgü immün farkları modellenmedi)."]
        per_host.append(d)

    headline = per_host[0]
    return {
        "tool": "netmhcpan",
        "bepipred": bool(b_eps),
        "n_bepitopes": len(b_eps),
        "headline": {"verdict": headline["verdict"], "score": headline["score"],
                     "host": headline["host"]},
        "per_host": per_host,
        "note": "In silico immünojenisite ÖNGÖRÜSÜ (garanti değil); eşikler operasyonel "
                "(OLUR≥60/ZAYIF≥35). Motor kinetiği peptid-aşısı held-out'unda doğrulandı.",
    }
