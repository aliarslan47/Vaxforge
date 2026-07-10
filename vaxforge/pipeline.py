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

from . import (citations, discovery, epitope, funnel, iedb_match, ingest,
               population, report, scoring, survival)
from .config_loader import ThresholdConfig, flatten_for_report
from .detect import Detection
from .hosts import HostRegistry
from .plan import build_plan, plan_table


def _ev(phase, status, msg, data=None):
    return {"phase": phase, "status": status, "msg": msg, "data": data}


def run(path, det: Detection, cfg: ThresholdConfig, profile: str,
        host_names: list[str] | None = None,
        overrides: dict | None = None, has_gpu: bool = False,
        outdir: str | Path = "outputs", host_registry: HostRegistry | None = None,
        organism_taxon: str | None = None):
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
        "n_raw": det.num_seqs,          # dosyadaki ham dizi/CDS sayısı
        "molecule": det.molecule,
    }

    # 1) Ingest (prep dahil): CDS -> protein çevirisi / ORF vb.
    unit = "CDS" if det.molecule == "cds" else ("okuma" if det.molecule == "reads" else "dizi")
    yield _ev("ingest", "running", f"Girdi: {det.num_seqs} {unit}. Proteom hazırlanıyor…")
    proteins = ingest.load_proteins(path, det)
    meta["n_input"] = len(proteins)
    yield _ev("ingest", "done",
              f"{det.num_seqs} {unit} → {len(proteins)} protein elde edildi", {"n": len(proteins)})
    if not proteins:
        yield _ev("__error__", "error", "Hiç protein elde edilemedi.")
        return

    # 2) Discovery — SKORLAMA (NERVE2 gibi), eleme DEĞİL: tüm proteinler geçer,
    # her birine virülans skoru atanır (adaylık puanına bileşen). Gerçek eleme
    # funnel'da (lokalizasyon/TM/homoloji). VFDB'de olmayan yüzey antijenleri kaybolmaz.
    yield _ev("discovery", "running", "Virülans skoru (küratörlü VFDB / anahtar-kelime)…")
    proteins = discovery.run(proteins, resolved["discovery_vfdb"], profile=profile)
    meta["n_discovery"] = len(proteins)
    n_vf = sum(1 for pr in proteins if pr.annotations.get("vf_hit"))
    yield _ev("discovery", "done",
              f"{len(proteins)} proteine virülans skoru atandı "
              f"({n_vf} VFDB kanıtlı) — eleme funnel'da", {"n": len(proteins), "vf_hit": n_vf})
    if not proteins:
        yield _ev("__error__", "error", "İşlenecek protein yok.")
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
    nB = sum(1 for p in peptides if p.kind == "B")
    nI = sum(1 for p in peptides if p.kind == "MHC-I")
    nII = sum(1 for p in peptides if p.kind == "MHC-II")
    yield _ev("mhc_panel", "done",
              f"Sliding-window ile {len(peptides)} epitop üretildi "
              f"(B={nB}, MHC-I={nI}, MHC-II={nII}) · {len(hosts)} konak" + tag,
              {"n": len(peptides), "hosts": [h.name for h in hosts]})

    # 5) Survival — KADEMELİ eleme (sayı adım adım azalır)
    yield _ev("survival", "running", f"{len(peptides)} peptit sağ kalım elemesine giriyor…")
    peptides, ssum = survival.run(peptides, resolved)
    meta["n_after_allergen"] = ssum["alerjenite_sonrasi"]
    meta["n_after_toxicity"] = ssum["toksisite_sonrasi"]
    yield _ev("survival_allergen", "done",
              f"Alerjenite ({ssum['alerjenite']}): {ssum['girdi_peptit']} → "
              f"{ssum['alerjenite_sonrasi']} peptit ({ssum['alerjen_elenen']} alerjen elendi)", ssum)
    yield _ev("survival_toxicity", "done",
              f"Toksisite ({ssum['toksisite']}): {ssum['alerjenite_sonrasi']} → "
              f"{ssum['toksisite_sonrasi']} peptit ({ssum['toksik_elenen']} toksik elendi)", ssum)

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

    # 6c) Popülasyon kapsamı (IEDB) — aday setinin konak/bölge kapsamı.
    yield _ev("population", "running", "Popülasyon kapsamı hesaplanıyor (IEDB HLA frekansları)…")
    popcov = population.for_candidates(peptides, hosts)
    meta["population_coverage"] = popcov
    hcov = popcov.get("hosts", {}).get("human", {})
    world = (hcov.get("mhc_i") or {}).get("World") if isinstance(hcov.get("mhc_i"), dict) else None
    msg = (f"İnsan MHC-I dünya kapsamı: {world['coverage']}%" if world
           else "İnsan HLA verisi yok/insan konak seçilmedi" if popcov.get("available")
           else "IEDB Population Coverage aracı kurulu değil")
    yield _ev("population", "done", msg, {"hosts": list(popcov.get("hosts", {}))})

    # 6d) IEDB literatür/bilinen-epitop taraması (SALT YORUM — skoru değiştirmez).
    # Adaylar deneysel doğrulanmış epitoplarla eşleştirilir; organizma taxon'u
    # verilirse o setle (recall benchmark dahil), yoksa canlı substring yedeğiyle.
    yield _ev("iedb", "running", "Aday peptitler IEDB'de (bilinen epitoplar/literatür) taranıyor…")
    im_summary = iedb_match.annotate_candidates(peptides, taxon_iri=organism_taxon)
    meta["iedb_match"] = im_summary
    if im_summary.get("available"):
        bm = im_summary.get("benchmark") or {}
        rec = f" · recall={bm['recall']}" if bm.get("recall") is not None else ""
        msg = (f"{im_summary['n_matched']}/{im_summary['n_candidates']} aday "
               f"bilinen IEDB epitobuyla eşleşti ({im_summary.get('source','')})" + rec)
    else:
        msg = im_summary.get("note", "IEDB taraması yapılamadı")
    yield _ev("iedb", "done", msg, {"n_matched": im_summary.get("n_matched", 0)})

    # 7) Rapor
    yield _ev("report", "running", "Rapor ve veri paketi üretiliyor…")
    run_dir = Path(outdir) / f"run_{meta['timestamp'].replace(':', '-')}"
    paths = report.write_package(run_dir, peptides, meta)
    yield _ev("report", "done", f"Çıktılar: {run_dir}", {k: str(v) for k, v in paths.items()})

    yield _ev("__result__", "done", "Tamamlandı", {
        "peptides": peptides, "paths": paths, "meta": meta,
    })
