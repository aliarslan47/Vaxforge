"""VaxForge web backend — mevcut vaxforge pipeline'ını saran ince katman.

Bilim yeniden yazılmaz: burada yalnız `vaxforge.pipeline.run()` generator'ı
çağrılır, config/host kayıtları yüklenir ve `outputs/` altındaki run'lar okunur.
Streamlit `app.py` ile birebir aynı çağrı sözleşmesi kullanılır.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterator

# Repo kökü: web/backend/runner.py -> ../../  (vaxforge paketi ve outputs/ burada)
REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = REPO_ROOT / "outputs"

# vaxforge paketi repo kökünde; import edilebilmesi için sys.path'e eklenir.
import sys
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from vaxforge import pipeline  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402
from vaxforge.plan import build_plan, plan_table  # noqa: E402

# Config ve host kayıtları bir kez yüklenir (app.py:28-29 ile aynı).
CFG = ThresholdConfig.load()
HOSTS = HostRegistry.load()

# ── Pipeline adımı → GERÇEK aktif yöntem kataloğu ──────────────────────────
# Grid'in kodun gerçeğini yansıtması için: her satır bir pipeline adımı; o adımda
# GERÇEKTE çalışan aracı gösterir. Araç kuruluysa `real`, değilse o adımın gerçek
# `fallback` yöntemi (heuristik de olsa gerçek bir yöntemdir). "kurulu değil" değil.
# alan: (id, module, step_tr, step_en, real, fallback)   fallback=None → hep gerçek
_PIPELINE_STEPS = [
    ("discovery", "discovery", "Virülans keşfi", "Virulence discovery",
     "DIAMOND 2.1.9 + VFDB", "VFDB anahtar-kelime + yerel hizalama"),
    ("localization", "psortb", "Hücre-altı lokalizasyon", "Subcellular localization",
     "PSORTb 3.0 (Gram)", "Kyte-Doolittle + sinyal peptid"),
    ("tm", "tmhmm_local", "Transmembran topoloji", "Transmembrane topology",
     "TMHMM-2.0", "Kyte-Doolittle hidropati"),
    ("signal", "signalp", "Sinyal peptidi", "Signal peptide",
     "SignalP-5.0", "Heuristik sinyal skoru"),
    ("antigenicity", "iapred", "Antijenite", "Antigenicity",
     "IApred (Miles 2025)", "antigen_acc (VaxiJen-tarzı proxy)"),
    ("bcell", "bepipred", "B-hücre epitop", "B-cell epitope",
     "BepiPred-1.0", "Parker hidrofilisite"),
    ("mhc", "netmhc_local", "MHC-I/II bağlanma", "MHC-I/II binding",
     "NetMHCpan / NetMHCIIpan", "SMM/proxy bağlanma"),
    ("allergen", "allergen", "Alerjenite", "Allergenicity",
     "FAO/WHO Codex 6-mer", None),
    ("toxicity", "toxinpred", "Toksisite", "Toxicity",
     "ToxinPred2 (Hybrid)", "ToxinPred2 (Model 1, AAC+RF)"),
    ("solubility", "solubility", "Çözünürlük (MEV)", "Solubility (MEV)",
     "Protein-Sol", None),
    ("secstruct", "secstruct", "İkincil yapı (MEV)", "Secondary structure (MEV)",
     "S4PRED", None),
    ("disorder", "disorder", "Düzensizlik (MEV)", "Disorder (MEV)",
     "metapredict V3", None),
    ("population", "population", "Popülasyon kapsamı", "Population coverage",
     "IEDB Population Coverage", "hesaplanamadı"),
]


def _module_available(mod: str) -> bool:
    try:
        m = __import__(f"vaxforge.{mod}", fromlist=["available"])
        return bool(m.available())
    except Exception:
        return False


def pipeline_steps() -> list[dict]:
    """Her pipeline adımı için GERÇEKTE çalışan yöntemi döndürür (grid için)."""
    out = []
    for sid, mod, step_tr, step_en, real, fallback in _PIPELINE_STEPS:
        installed = _module_available(mod)
        # fallback None ise yöntem zaten yereldir → her zaman gerçek/kurulu say.
        active_installed = installed or fallback is None
        method = real if active_installed else fallback
        out.append({
            "id": sid, "module": mod,
            "step_tr": step_tr, "step_en": step_en,
            "method": method,
            "real": real, "fallback": fallback,
            "installed": active_installed,
        })
    return out


def tool_status() -> list[dict]:
    """Geriye-uyum: eski 'tools' alanı (module/label/available)."""
    return [{"module": s["module"], "label": s["method"], "available": s["installed"]}
            for s in pipeline_steps()]


def get_config() -> dict:
    """Frontend'in ihtiyaç duyduğu tüm seçenekler: profiller, konaklar, araçlar."""
    hosts = []
    for name in HOSTS.names():
        h = HOSTS.get(name)
        hosts.append({
            "name": name,
            "label": h.label,
            "source": getattr(h, "source", ""),
            "n_mhc_i": len(h.mhc_i),
            "n_mhc_ii": len(h.mhc_ii),
            "default": name in HOSTS.default_hosts,
        })
    try:
        from vaxforge import mev as _mev
        adjuvants = _mev.list_adjuvants()
    except Exception:
        adjuvants = []
    return {
        "profiles": CFG.profiles,
        "default_profile": CFG.default_profile,
        "hosts": hosts,
        "default_hosts": list(HOSTS.default_hosts),
        "tools": tool_status(),
        "steps": pipeline_steps(),
        "adjuvants": adjuvants,
        "default_adjuvant": "beta_defensin",
        "config_file": CFG.path.name,
        "hosts_file": HOSTS.path.name,
    }


def _run_summary(run_dir: Path, meta: dict) -> dict:
    """outputs/<run> için kart-özeti (liste görünümü)."""
    cands = meta.get("candidates") or []
    hosts = [h.get("label") or h.get("name") for h in meta.get("hosts", [])]
    return {
        "id": run_dir.name,
        "input": meta.get("input", "?"),
        "profile": meta.get("profile", "?"),
        "timestamp": meta.get("timestamp", ""),
        "molecule": meta.get("molecule", ""),
        "lang": meta.get("lang", "tr"),
        "hosts": hosts,
        "n_input": meta.get("n_input"),
        "n_candidates": len(cands),
        "has_pdf": (run_dir / "report.pdf").exists(),
        "has_html": (run_dir / "report.html").exists(),
        "has_xlsx": (run_dir / "candidates_full.xlsx").exists(),
    }


def list_runs() -> list[dict]:
    """outputs/ altındaki tüm run'ları run.json'dan özetler, en yeni önce."""
    runs = []
    if not OUTPUTS.exists():
        return runs
    for d in OUTPUTS.iterdir():
        if not d.is_dir():
            continue
        rj = d / "run.json"
        if not rj.exists():
            continue
        try:
            meta = json.loads(rj.read_text(encoding="utf-8"))
        except Exception:
            continue
        runs.append(_run_summary(d, meta))
    runs.sort(key=lambda r: r.get("timestamp", ""), reverse=True)
    return runs


def _clean_nan(obj: Any) -> Any:
    """NaN/Infinity → None (JSON uyumlu). pandas aday kayıtları NaN içerebilir;
    Starlette JSONResponse allow_nan=False ile serialize eder → aksi halde 500."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _clean_nan(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_nan(v) for v in obj]
    return obj


def get_run(run_id: str) -> dict | None:
    """Tek bir run'ın tam run.json'u (adaylar + meta gömülü)."""
    rj = _safe_run_dir(run_id)
    if rj is None:
        return None
    rj = rj / "run.json"
    if not rj.exists():
        return None
    try:
        # json.loads NaN'ı float('nan') olarak parse eder; API'ye vermeden temizle.
        return _clean_nan(json.loads(rj.read_text(encoding="utf-8")))
    except Exception:
        return None


def _safe_run_dir(run_id: str) -> Path | None:
    """Path-traversal'a karşı güvenli run dizini çözümü."""
    if not run_id or "/" in run_id or ".." in run_id:
        return None
    d = (OUTPUTS / run_id).resolve()
    if OUTPUTS.resolve() not in d.parents:
        return None
    return d if d.is_dir() else None


def delete_run(run_id: str) -> bool:
    """Bir koşu çıktısını (outputs/run_xxx) güvenle siler. Başarı → True."""
    d = _safe_run_dir(run_id)
    if d is None:
        return False
    import shutil as _sh
    _sh.rmtree(d, ignore_errors=True)
    return not d.exists()


_ALLOWED_FILES = {
    "report.html", "report.pdf", "candidates.csv",
    "candidates_full.xlsx", "top_peptides.fasta", "run.json",
}


def get_run_file(run_id: str, name: str) -> Path | None:
    """İndirme/serve için güvenli dosya yolu (yalnız beyaz-listedeki adlar)."""
    if name not in _ALLOWED_FILES:
        return None
    d = _safe_run_dir(run_id)
    if d is None:
        return None
    f = d / name
    return f if f.exists() else None


def run_pipeline(input_path: str, filename: str, profile: str,
                 host_names: list[str] | None, gram: str | None,
                 lang: str = "tr", adjuvant: str = "beta_defensin") -> Iterator[dict]:
    """pipeline.run() generator'ını sarar; her event'i dict olarak yield eder.

    Streamlit app.py:319-322 ile aynı çağrı. Son event `__result__` yerine
    frontend'e sade bir `{phase:"done", run_id, summary}` gönderilir.
    """
    det = detect(input_path)
    det.filename = filename

    # Girdi tanıma özetini ilk event olarak gönder (arayüz "Tanıma" kartı için).
    yield {
        "phase": "__detect__", "status": "info", "msg": det.summary,
        "data": {
            "fmt": det.fmt, "seq_type": det.seq_type, "molecule": det.molecule,
            "num_seqs": det.num_seqs, "avg_len": round(det.avg_len, 1),
            "min_len": det.min_len, "max_len": det.max_len,
            "is_gzipped": det.is_gzipped, "confident": det.confident,
            "notes": det.notes,
        },
    }

    # Planlanan adımlar (ilerleme çubuğu toplamı için).
    steps = build_plan(det, has_gpu=False)
    yield {
        "phase": "__plan__", "status": "info", "msg": f"{len(steps)} adım planlandı",
        "data": {"steps": plan_table(steps), "total": len(steps)},
    }

    gram_val = gram if profile == "bacteria" else None
    result_meta = None
    result_run_id = None

    for ev in pipeline.run(
        input_path, det, CFG, profile,
        host_names=host_names or None, overrides={},
        has_gpu=False, outdir=str(OUTPUTS), host_registry=HOSTS,
        organism_taxon=None, gram=gram_val, lang=lang, adjuvant=adjuvant,
    ):
        phase = ev["phase"]
        if phase == "__result__":
            data = ev.get("data") or {}
            result_meta = data.get("meta")
            paths = data.get("paths") or {}
            jpath = paths.get("json")
            if jpath:
                result_run_id = Path(str(jpath)).parent.name
            continue
        if phase == "__error__":
            yield {"phase": "__error__", "status": "error",
                   "msg": ev.get("msg", "Bilinmeyen hata"), "data": None}
            return
        # normal ilerleme event'i — data'yı JSON-güvenli tut
        yield {
            "phase": phase, "status": ev.get("status", ""),
            "msg": ev.get("msg", ""), "data": _jsonable(ev.get("data")),
        }

    # Bitiş: run_id + özet
    summary = None
    if result_run_id:
        run_dir = OUTPUTS / result_run_id
        if result_meta is not None:
            summary = _run_summary(run_dir, result_meta)
    yield {
        "phase": "__done__", "status": "done", "msg": "Tamamlandı",
        "data": {"run_id": result_run_id, "summary": summary},
    }


def _jsonable(obj: Any) -> Any:
    """Event data'sını JSON'a çevrilebilir hale getir (tuple/set/obj güvenliği)."""
    try:
        json.dumps(obj)
        return obj
    except (TypeError, ValueError):
        if isinstance(obj, dict):
            return {str(k): _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple, set)):
            return [_jsonable(v) for v in obj]
        return str(obj)
