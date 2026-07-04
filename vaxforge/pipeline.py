"""Orkestratör — tüm adımları sırayla koşar ve ilerleme olayları yield eder.

Kullanım:
    for ev in run(path, det, cfg, "bacteria"):
        if ev["phase"] == "__result__":
            results = ev["data"]
        else:
            print(ev["phase"], ev["status"], ev["msg"])

Ağır/deferred adımlar (structure, docking_md) GPU yoksa atlanır ve öyle işaretlenir.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from . import (citations, discovery, epitope, funnel, ingest, mrna, report,
               scoring, survival)
from .config_loader import ThresholdConfig, flatten_for_report
from .detect import Detection
from .hosts import HostRegistry
from .plan import build_plan, plan_table


def _ev(phase, status, msg, data=None):
    return {"phase": phase, "status": status, "msg": msg, "data": data}


def run(path, det: Detection, cfg: ThresholdConfig, profile: str,
        host_names: list[str] | None = None,
        overrides: dict | None = None, has_gpu: bool = False,
        outdir: str | Path = "outputs", host_registry: HostRegistry | None = None):
    resolved = cfg.resolve(profile, overrides)
    steps = build_plan(det, has_gpu=has_gpu)
    reg = host_registry or HostRegistry.load()
    host_names = host_names or reg.default_hosts
    hosts = [reg.get(h) for h in host_names if h in reg.hosts]
    meta = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input": det.filename, "profile": profile,
        "hosts": [{"name": h.name, "label": h.label, "source": h.source,
                   "n_mhc_i": len(h.mhc_i), "n_mhc_ii": len(h.mhc_ii),
                   "predictors": h.predictors} for h in hosts],
        "thresholds": flatten_for_report(resolved),
        "plan": plan_table(steps),
        "citations": citations.for_report(),
    }

    # 1) Ingest (prep dahil)
    yield _ev("ingest", "running", "Proteom hazırlanıyor…")
    proteins = ingest.load_proteins(path, det)
    meta["n_input"] = len(proteins)
    yield _ev("ingest", "done", f"{len(proteins)} aday protein elde edildi", {"n": len(proteins)})
    if not proteins:
        yield _ev("__error__", "error", "Hiç protein elde edilemedi.")
        return

    # 2) Discovery
    yield _ev("discovery", "running", "Küratörlü DB / anahtar-kelime taraması…")
    proteins = discovery.run(proteins, resolved["discovery_vfdb"])
    meta["n_discovery"] = len(proteins)
    yield _ev("discovery", "done", f"{len(proteins)} virülans/hedef adayı kaldı", {"n": len(proteins)})
    if not proteins:
        yield _ev("__error__", "error", "Keşif adımından aday çıkmadı (eşikleri gevşetin).")
        return

    # 3) Funnel
    yield _ev("funnel", "running", "Antijen eleme hunisi (melez)…")
    proteins, fsum = funnel.run(proteins, resolved, profile=profile)
    meta["n_funnel"] = len(proteins)
    yield _ev("funnel", "done", f"{len(proteins)} antijen huniden geçti", fsum)
    if not proteins:
        yield _ev("__error__", "error", "Huniden aday çıkmadı.")
        return

    # 4) Epitope + çok-KONAKLI MHC (I + II)
    host_lbls = ", ".join(h.label for h in hosts) or "—"
    yield _ev("epitope", "running", f"B/T-hücre epitop tahmini · konaklar: {host_lbls}")
    peptides = epitope.run(proteins, resolved, hosts)
    meta["n_epitope"] = len(peptides)
    real_i = any("GERÇEK" in p.methods.get("mhc_score", "") for p in peptides if p.kind == "MHC-I")
    real_ii = any("GERÇEK" in p.methods.get("mhc_score", "") for p in peptides if p.kind == "MHC-II")
    tag = " · " + " ".join(["MHC-I=GERÇEK" if real_i else "MHC-I=proxy",
                            "MHC-II=GERÇEK" if real_ii else "MHC-II=proxy"])
    yield _ev("mhc_panel", "done",
              f"{len(peptides)} epitop (B + MHC-I + MHC-II) · {len(hosts)} konak" + tag,
              {"n": len(peptides), "hosts": [h.name for h in hosts]})

    # 5) Survival
    yield _ev("survival", "running", "Sağ kalım elemesi (alerjenite/toksisite)…")
    peptides, ssum = survival.run(peptides, resolved)
    yield _ev("survival", "done", f"{len(peptides)} peptit hayatta kaldı", ssum)

    # 6) Deferred ağır adımlar
    for hid, title in [("structure", "Peptit-MHC yapısı (AlphaFold)"),
                       ("docking_md", "Docking + MD")]:
        st = next((s for s in steps if s.id == hid), None)
        if st and st.status == "deferred":
            yield _ev(hid, "deferred", f"{title}: GPU yok → ertelendi (uzak worker)")

    # 6b) Scoring
    yield _ev("scoring", "running", "Adaylık puanı hesaplanıyor…")
    peptides = scoring.score(peptides, cfg.candidacy_weights())
    meta["n_survivors"] = len(peptides)
    yield _ev("scoring", "done", f"{len(peptides)} aday sıralandı",
              {"top": [(p.seq, p.kind, p.candidacy) for p in peptides[:5]]})

    # 7) mRNA konstrükt
    yield _ev("mrna", "running", "Çok-epitoplu mRNA konstrüktü kuruluyor…")
    construct = mrna.build(peptides)
    yield _ev("mrna", "done", f"Konstrükt hazır: {construct.metrics['mrna_len']} nt, "
              f"GC {construct.metrics['gc_percent']}%", construct.metrics)

    # 8) Rapor
    yield _ev("report", "running", "Rapor ve veri paketi üretiliyor…")
    run_dir = Path(outdir) / f"run_{meta['timestamp'].replace(':', '-')}"
    paths = report.write_package(run_dir, peptides, construct, meta)
    yield _ev("report", "done", f"Çıktılar: {run_dir}", {k: str(v) for k, v in paths.items()})

    yield _ev("__result__", "done", "Tamamlandı", {
        "peptides": peptides, "construct": construct, "paths": paths, "meta": meta,
    })
