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

# Araç durumu için modül listesi (app.py:288-297 deseni).
_TOOL_MODULES = [
    ("discovery", "DIAMOND+VFDB"), ("deeploc", "DeepLoc"),
    ("tmhmm_local", "TMHMM"), ("signalp", "SignalP"),
    ("iapred", "IApred"), ("bepipred", "BepiPred"),
    ("netmhc_local", "NetMHCpan"), ("toxinpred", "ToxinPred2"),
]


def tool_status() -> list[dict]:
    """Her gerçek aracın kurulu olup olmadığını döndürür (available())."""
    out = []
    for mod, label in _TOOL_MODULES:
        ok = False
        try:
            m = __import__(f"vaxforge.{mod}", fromlist=["available"])
            ok = bool(m.available())
        except Exception:
            ok = False
        out.append({"module": mod, "label": label, "available": ok})
    return out


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
    return {
        "profiles": CFG.profiles,
        "default_profile": CFG.default_profile,
        "hosts": hosts,
        "default_hosts": list(HOSTS.default_hosts),
        "tools": tool_status(),
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
                 lang: str = "tr") -> Iterator[dict]:
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
        organism_taxon=None, gram=gram_val, lang=lang,
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
