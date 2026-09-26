"""Orkestratör — tüm adımları sırayla koşar ve ilerleme olayları yield eder.

Kullanım:
    for ev in run(path, det, cfg, "bacteria"):
        if ev["phase"] == "__result__":
            results = ev["data"]
        else:
            print(ev["phase"], ev["status"], ev["msg"])

(Yapısal doğrulama + docking/MD adımları şimdilik çıkarıldı; odak aday-belirleme.)
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from . import (citations, conservation, discovery, epitope, funnel, iedb_match,
               ingest, mev, mimicry, population, report, scoring, survival)
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
        organism_taxon: str | None = None, gram: str | None = None,
        lang: str = "tr", adjuvant: str = "beta_defensin",
        strain_paths: list[str] | None = None):
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
        "gram": gram if profile == "bacteria" else None,
        "lang": lang,
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
    proteins, fsum = funnel.run(proteins, resolved, profile=profile,
                                hosts=hosts, gram=gram)
    meta["n_funnel"] = len(proteins)
    meta["funnel_methods"] = fsum   # rapor/UI "Kullanılan yöntemler" bölümü (Faz 1)
    yield _ev("funnel", "done", f"{len(proteins)} antijen huniden geçti", fsum)
    if not proteins:
        yield _ev("__error__", "error", "Huniden aday çıkmadı.")
        return

    # 3b) Konservasyon — suşlar arası korunmuşluk.
    # Yalnız funnel-survivor adaylarda; diamond ortolog + BLOSUM62 hizalama (hafif).
    cons_min = float(resolved["conservation"].params["min_percent"].value) \
        if "conservation" in resolved else 80.0
    # FAZ 8: kullanıcı suş vermediyse ama organizma taxon'u varsa NCBI'den otomatik suş çek.
    auto_strain_note = None
    if not strain_paths and organism_taxon:
        try:
            from . import strain_fetch
            if strain_fetch.available(organism_taxon):
                yield _ev("conservation", "running",
                          "Suş verisi yok — NCBI Datasets'ten otomatik suş çekiliyor…")
                fetched = strain_fetch.fetch_strains(organism_taxon, n=3)
                if fetched.get("ok"):
                    strain_paths = fetched["paths"]
                    auto_strain_note = fetched["note"]
                    meta["conservation_auto_strains"] = {
                        "accessions": fetched["accessions"], "note": fetched["note"]}
                else:
                    meta["conservation_auto_strains"] = {"accessions": [], "note": fetched.get("note")}
        except Exception as e:
            meta["conservation_auto_strains"] = {"accessions": [], "note": f"otomatik çekme atlandı ({type(e).__name__})"}
    if conservation.available(strain_paths):
        yield _ev("conservation", "running",
                  f"Suş konservasyonu ({len(strain_paths)} suş) hesaplanıyor…")
        proteins, csum = conservation.run(proteins, strain_paths, min_percent=cons_min)
        if auto_strain_note:
            csum["auto_fetched"] = True
            csum["auto_note"] = auto_strain_note
        meta["conservation"] = csum
        src = f" · {auto_strain_note}" if auto_strain_note else ""
        yield _ev("conservation", "done",
                  f"{csum.get('n_conserved', 0)}/{len(proteins)} aday ≥%{cons_min} korunmuş "
                  f"({csum.get('n_strains', 0)} suş){src}", csum)
    else:
        proteins, csum = conservation.run(proteins, strain_paths, min_percent=cons_min)
        meta["conservation"] = csum
        yield _ev("conservation", "done",
                  "Konservasyon hesaplanmadı (suş verisi verilmedi)", csum)

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

    # 5b) Moleküler mimikri — epitop-seviyesi konak self-eşleşmesi (YUMUŞAK: elemez,
    # işaretler + skorda hafif ceza). Protein-seviyesi homoloji funnel'da zaten sert.
    mim_k = int(resolved["molecular_mimicry"].params["k"].value) \
        if "molecular_mimicry" in resolved else 9
    mim = mimicry.predict([p.seq for p in peptides], [h.name for h in hosts], k=mim_k)
    n_mim = 0       # tam eşleşme (sert mimikri)
    n_near = 0      # kısmi eşleşme (tek konservatif ikame — BLOSUM62)
    for p in peptides:
        d = mim.get(p.seq)
        if d:
            p.metrics["self_mimicry"] = d["self_match"]
            p.metrics["self_mimicry_host"] = d["host"]
            p.metrics["self_mimicry_match"] = d["match"]
            p.metrics["self_mimicry_near"] = d["near_match"]
            p.metrics["self_mimicry_near_sub"] = d["near_sub"]
            p.metrics["self_mimicry_near_blosum"] = d["near_blosum"]
            p.methods["self_mimicry"] = (f"konak self {mim_k}-mer tam eşleşme + "
                                         "BLOSUM62 Hamming-1 kısmi (konservatif ikame)")
            if d["self_match"]:
                n_mim += 1
            elif d["near_match"]:
                n_near += 1
    meta["mimicry"] = {
        "available": True, "k": mim_k,
        "n_candidates": len(peptides), "n_exact": n_mim, "n_near": n_near,
        "hosts": [h.name for h in hosts],
        "method": ("tam k-mer eşleşme (sert mimikri) + BLOSUM62 Hamming-1 konservatif "
                   "ikame (kısmi mimikri) — YUMUŞAK: işaretler, elemez"),
    }
    yield _ev("mimicry", "done",
              f"Moleküler mimikri: {n_mim} tam + {n_near} kısmi (BLOSUM62) / {len(peptides)} "
              f"peptit konak self {mim_k}-mer eşleşmesi (yumuşak uyarı)",
              {"n_exact": n_mim, "n_near": n_near, "k": mim_k})

    # 5c) Glikozilasyon — glikan-maskeli B-epitop filtresi (YUMUŞAK: işaretler, elemez).
    # N-X-S/T sequon ile örtüşen B-epitobu in vivo glikanla kaplı olabilir → antikor erişemez.
    yield _ev("glycosylation", "running",
              "Glikozilasyon taraması (N-X-S/T sequon → glikan-maskeli epitop)…")
    try:
        from . import glyco
        gsum = glyco.annotate_peptides(peptides)
        meta["glycosylation"] = gsum
        yield _ev("glycosylation", "done",
                  f"Glikozilasyon: {gsum['n_bcell_masked']}/{gsum['n_bcell']} B-epitobu "
                  f"glikan-maskeli sequon içeriyor (yumuşak uyarı)", gsum)
    except Exception as e:  # opsiyonel: pipeline'ı düşürmesin
        meta["glycosylation"] = None
        yield _ev("glycosylation", "done", f"Glikozilasyon atlandı ({type(e).__name__})", {})

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

    # 6e) MEV — çok-epitoplu aşı konstruktu inşası + karakterizasyonu.
    # Top epitoplar linker (AAY/GPGPG/KK/EAAAK) + β-defensin adjuvan ile
    # birleştirilir; konstrukt-seviyesi fizikokimyasal/immünoinformatik özellikler.
    yield _ev("mev", "running", "Çok-epitoplu aşı konstruktu (MEV) inşa ediliyor…")
    try:
        mev_out = mev.run(peptides, adjuvant=adjuvant)
        meta["mev"] = mev_out
        mp = mev_out["properties"]
        yield _ev("mev", "done",
                  f"MEV: {len(mev_out['sequence'])} aa · MHC-I={mev_out['n_by_kind']['MHC-I']}, "
                  f"MHC-II={mev_out['n_by_kind']['MHC-II']}, B={mev_out['n_by_kind']['B']} · "
                  f"MW {mp['mw_kda']} kDa · antijenite {mp['antigenicity']['score']}",
                  {"length": len(mev_out["sequence"]), "adjuvant": mev_out["adjuvant"]})
    except Exception as e:  # MEV opsiyonel: hata pipeline'ı düşürmesin
        meta["mev"] = None
        yield _ev("mev", "done", f"MEV atlandı ({type(e).__name__})", {})

    # 6e-2) YAPISAL DOĞRULAMA — ESMFold (ESMAtlas API) → pLDDT (katlanma güveni) +
    # RSA (yüzey erişilebilirliği) + konformasyonel/yüzey B-epitop yamaları. MEV
    # konstruktun DÜZGÜN KATLANIP KATLANMADIĞINI ve B-epitopların yüzeye açık olup
    # olmadığını sayısal doğrular. Ağ yok / >400 aa → dürüst ARAÇ_YOK (uydurma yok).
    mev_seq_st = (meta.get("mev") or {}).get("sequence")
    if mev_seq_st:
        yield _ev("structure", "running",
                  "Yapısal doğrulama (ESMFold · ESMAtlas → pLDDT / RSA / konformasyonel B)…")
        try:
            from . import bcell_conf, rsa as rsa_mod, structure
            st = None
            note = ""
            if not structure.available():
                note = "ESMAtlas erişilemedi (çevrimdışı) — yapı öngörülmedi (ARAÇ_YOK)"
            else:
                pred = structure.predict(mev_seq_st)
                if not pred:
                    note = (f"MEV {len(mev_seq_st)} aa > {structure.MAX_LEN} aa (ESMAtlas sınırı) "
                            "— yapı öngörülmedi" if len(mev_seq_st) > structure.MAX_LEN
                            else "ESMFold öngörüsü alınamadı (API hatası)")
                else:
                    rsa_list = rsa_mod.rsa_from_pdb(pred["pdb"])
                    surf = rsa_mod.surface_summary(rsa_list) if rsa_list else {}
                    conf = bcell_conf.predict(mev_seq_st, rsa_list) if rsa_list else None
                    st = {
                        "available": True, "source": "MEV konstrukt",
                        "tool": "ESMFold (ESMAtlas API) — öngörü, deneysel değil",
                        "n_res": pred["n_res"], "plddt_mean": pred["plddt_mean"],
                        "plddt_min": pred["plddt_min"], "surface": surf,
                        "conformational_b": conf,
                    }
            meta["structure"] = st or {"available": False, "note": note}
            if st:
                cb = st.get("conformational_b") or {}
                yield _ev("structure", "done",
                          f"MEV yapısı: pLDDT {st['plddt_mean']} (min {st['plddt_min']}) · "
                          f"yüzey %{st['surface'].get('percent_surface', '—')} · "
                          f"konformasyonel B yaması {cb.get('n_patches', 0)}", st)
            else:
                yield _ev("structure", "done", f"Yapısal doğrulama atlandı — {note}",
                          {"available": False})
        except Exception as e:  # opsiyonel: pipeline'ı düşürmesin
            meta["structure"] = None
            yield _ev("structure", "done", f"Yapısal doğrulama atlandı ({type(e).__name__})", {})

    # 6e-3) IN SILICO KLONLAMA — MEV proteini E. coli kodonlarına uyarlanır (CAI),
    # GC% + iç restriksiyon bölgeleri + pET vektör/etiket önerisi. Üretilebilirlik halkası.
    mev_seq_cl = (meta.get("mev") or {}).get("sequence")
    if mev_seq_cl:
        yield _ev("cloning", "running",
                  "In silico klonlama (kodon-optimizasyon · CAI · restriksiyon · vektör)…")
        try:
            from . import cloning
            cl = cloning.run(mev_seq_cl, host="ecoli")
            meta["cloning"] = cl
            if cl:
                yield _ev("cloning", "done",
                          f"Klonlama: {cl['dna_len_bp']} bp · CAI {cl['cai']} · GC %{cl['gc_percent']} · "
                          f"iç kesim {cl['n_internal_sites']} · {cl['suggested_vector']}", cl)
            else:
                yield _ev("cloning", "done", "Klonlama atlandı (dizi yok)", {})
        except Exception as e:  # opsiyonel: pipeline'ı düşürmesin
            meta["cloning"] = None
            yield _ev("cloning", "done", f"Klonlama atlandı ({type(e).__name__})", {})

    # 6f) İMMÜNİZASYON TESTİ — ImmForge: aday MEV in silico immünize edilir → OLUR/OLMAZ.
    # Ters aşı bilişimi döngü kapanışı: konak MHC allelleri için GERÇEK NetMHCpan → motor → karar.
    mev_seq = (meta.get("mev") or {}).get("sequence")
    if mev_seq:
        yield _ev("immunization", "running",
                  "In silico immünizasyon testi (ImmForge · NetMHCpan → motor)…")
        try:
            from . import immforge_score
            verdict = immforge_score.score_construct(mev_seq, hosts)
            meta["immunogenicity"] = verdict
            hl = verdict.get("headline") or {"verdict": verdict.get("verdict", "—")}
            sc = f" (skor {hl['score']:.0f}/100)" if hl.get("score") is not None else ""
            yield _ev("immunization", "done",
                      f"İmmünizasyon kararı: {hl.get('verdict', '—')}{sc}", verdict)
        except Exception as e:  # opsiyonel: pipeline'ı düşürmesin
            meta["immunogenicity"] = None
            yield _ev("immunization", "done",
                      f"İmmünizasyon testi atlandı ({type(e).__name__})", {})

    # 7) Rapor
    yield _ev("report", "running", "Rapor ve veri paketi üretiliyor…")
    run_dir = Path(outdir) / f"run_{meta['timestamp'].replace(':', '-')}"
    paths = report.write_package(run_dir, peptides, meta)
    yield _ev("report", "done", f"Çıktılar: {run_dir}", {k: str(v) for k, v in paths.items()})

    yield _ev("__result__", "done", "Tamamlandı", {
        "peptides": peptides, "paths": paths, "meta": meta,
    })
