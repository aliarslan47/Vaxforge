"""Adım 3 — Antijen eleme hunisi (melez: sert filtreler + skorlama).

Sert filtreler (geçemeyen elenir): lokalizasyon, transmembran sayısı,
insan homolojisi (ZORUNLU — ama bu makinede insan proteom DB'si yoksa
'çalıştırılmadı' olarak işaretlenir, sessizce elemez).
Yumuşak (skorlama): sinyal peptidi, antijenite, konservasyon, adezin.

Tüm eşikler config'ten gelir. Heuristik yöntemler 'method_*' ile etiketlenir;
gerçek araçlar (DeepLoc/SignalP/DeepTMHMM/VaxiJen/blastp) takılınca değişir.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from . import antigen_acc, iapred, psortb, signalp, tmhmm_local
from .config_loader import ResolvedTool
from .hosts import Host
from .models import ProteinRecord
from .sequtils import (antigenicity_proxy, count_tm_helices,
                       predict_localization, sanitize, signal_peptide_prob)

_TOOLS = Path(__file__).resolve().parent.parent / "tools"
_DIAMOND = _TOOLS / "diamond"
_HUMAN_DB = _TOOLS / "db" / "human.dmnd"


def human_db_available() -> bool:
    return _DIAMOND.exists() and _HUMAN_DB.exists()


def _host_db_path(host: Host) -> Path | None:
    """Konağın homology_db yolunu çözer (tools/ altına göre). Yoksa None."""
    if not host.homology_db:
        return None
    p = _TOOLS / host.homology_db
    return p if p.exists() else None


def _homology_vs_db(proteins: list[ProteinRecord], db_path: Path,
                    evalue: float) -> dict[str, float]:
    """qseqid -> verilen proteom DB'sine en iyi %kimlik (diamond blastp). Hata -> boş."""
    if not _DIAMOND.exists():
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for pr in proteins:
            fh.write(f">{pr.id}\n{sanitize(pr.seq)}\n")
        qpath = fh.name
    try:
        cmd = [str(_DIAMOND), "blastp", "-q", qpath, "-d", str(db_path),
               "--outfmt", "6", "qseqid", "pident", "--max-target-seqs", "1",
               "--evalue", str(evalue), "--quiet"]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception:
        Path(qpath).unlink(missing_ok=True)
        return {}
    Path(qpath).unlink(missing_ok=True)
    best: dict[str, float] = {}
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            qid, pid = p[0], float(p[1])
            best[qid] = max(best.get(qid, 0.0), pid)
    return best


def _self_homology(proteins: list[ProteinRecord], hosts: list[Host],
                   evalue: float) -> tuple[dict[str, dict], list[str], list[str]]:
    """Host-driven otoimmünite taraması.

    Seçili HER konağın kendi proteomuna karşı diamond blastp koşar. Bir protein
    HERHANGİ bir seçili konağın proteomuna eşik üstü benzerse otoimmünite riski.

    Döner: (per_protein, taranan_konaklar, DB'siz_konaklar)
      per_protein[pid] = {'max_pident': float, 'by_host': {host_label: pident}}
    """
    scanned: list[str] = []
    missing: list[str] = []
    per_host: dict[str, dict[str, float]] = {}
    for h in hosts:
        db = _host_db_path(h)
        if db is None:
            missing.append(h.label)
            continue
        scanned.append(h.label)
        per_host[h.label] = _homology_vs_db(proteins, db, evalue)
    per_protein: dict[str, dict] = {}
    for pr in proteins:
        by_host = {lbl: best.get(pr.id, 0.0) for lbl, best in per_host.items()}
        per_protein[pr.id] = {
            "max_pident": max(by_host.values()) if by_host else 0.0,
            "by_host": by_host,
        }
    return per_protein, scanned, missing

# heuristik lokalizasyon çıktısını config sözlüğüne eşle
_LOC_MAP = {"secreted": "extracellular", "outer_membrane": "outer_membrane",
            "membrane": "membrane", "cytoplasm": "cytoplasm"}


def run(proteins: list[ProteinRecord], tools: dict[str, ResolvedTool],
        profile: str = "bacteria", hosts: list[Host] | None = None,
        gram: str | None = None) -> tuple[list[ProteinRecord], dict]:
    """Huniyi uygula. (hayatta kalanlar, özet) döndürür.

    Lokalizasyon (bakteride sert filtre): önce PSORTb (prokaryota-özel, Gram-tipine
    göre; en doğru) → yoksa DeepLoc (ökaryot-yanlı, yedek) → yoksa heuristik.
    Virüs/parazit profilinde lokalizasyon eleme yapmaz ('çalıştırılmadı').

    Otoimmünite (self-homoloji): HOST-DRIVEN — seçili her konağın kendi proteomuna
    karşı taranır (yalnız insan değil). DB'si olmayan konak için dürüst 'atlandı'.
    """
    hosts = hosts or []
    loc_hard = profile == "bacteria"
    loc_allowed = set(tools["localization"].params["allowed"].value)
    max_tm = int(tools["transmembrane"].params["max_helices"].value)
    vaxijen_thr = float(tools["antigenicity_vaxijen"].params["threshold"].value)
    sp_min = float(tools["signal_peptide"].params["min_probability"].value)
    hh_id = float(tools["human_homology"].params["min_identity"].value)
    hh_ev = float(tools["human_homology"].params["evalue"].value)

    survivors: list[ProteinRecord] = []
    # -- host-driven self-homoloji (otoimmünite) — seçili her konak
    hh_data, hh_scanned, hh_missing = _self_homology(proteins, hosts, hh_ev)
    hh_available = bool(hh_scanned)
    seq_pairs = [(pr.id, sanitize(pr.seq)) for pr in proteins]
    tm_available = tmhmm_local.available()
    tm_real = tmhmm_local.predict(seq_pairs) if tm_available else {}
    # -- lokalizasyon: bakteride PSORTb (Gram) > heuristik (Kyte-Doolittle + sinyal)
    psortb_ok = loc_hard and gram in ("positive", "negative") and psortb.available()
    psortb_real = psortb.predict(seq_pairs, gram) if psortb_ok else {}
    sp_available = signalp.available()
    sp_real = signalp.predict(seq_pairs, profile=profile) if sp_available else {}
    ia_available = iapred.available()
    ia_res = iapred.predict(seq_pairs) if ia_available else {}
    ag_available = antigen_acc.available()
    for pr in proteins:
        # -- lokalizasyon: PSORTb (prokaryota-özel) > heuristik (KD+sinyal)
        if psortb_ok and pr.id in psortb_real:
            loc = psortb_real[pr.id]["localization"]
            pr.annotations["localization_raw"] = psortb_real[pr.id]["raw"]
            pr.annotations["localization_score"] = psortb_real[pr.id]["score"]
            pr.annotations["method_localization"] = f"GERÇEK (PSORTb, Gram{'+' if gram=='positive' else '−'})"
        else:
            loc = _LOC_MAP.get(predict_localization(pr.seq), "cytoplasm")
            pr.annotations["method_localization"] = "heuristik (KD+sinyal)"
        # -- transmembran (gerçek TMHMM, yoksa Kyte-Doolittle)
        if tm_available and pr.id in tm_real:
            tm = tm_real[pr.id]
            pr.annotations["method_tm"] = "GERÇEK (TMHMM-2.0, yerel)"
        else:
            tm = count_tm_helices(pr.seq)
            pr.annotations["method_tm"] = "Kyte-Doolittle hidropati (klasik yöntem)"
        # -- sinyal peptidi (gerçek SignalP, yoksa heuristik)
        if sp_available and pr.id in sp_real:
            sp = sp_real[pr.id]["sp_prob"]
            pr.annotations["method_signalp"] = "GERÇEK (SignalP-5.0)"
        else:
            sp = signal_peptide_prob(pr.seq)
            pr.annotations["method_signalp"] = "heuristik"
        if ia_available and pr.id in ia_res:
            ag = ia_res[pr.id]["norm"]
            pr.annotations["antigenicity_raw"] = ia_res[pr.id]["score"]
            pr.annotations["antigenicity_category"] = ia_res[pr.id]["category"]
            ag_method = "GERÇEK (IApred — VaxiJen'den iyi)"
        elif ag_available:
            ag = antigen_acc.predict(pr.seq)
            ag_method = "z-ACC + lojistik (VaxiJen yöntemi, resmi model değil)"
        else:
            ag = antigenicity_proxy(pr.seq)
            ag_method = "heuristik proxy (VaxiJen DEĞİL)"

        pr.annotations.update({
            "localization": loc, "tm_helices": tm, "signalp": sp,
            "antigenicity": ag,
            "method_antigenicity": ag_method,
        })

        # -- sert (yalnız bakteride): lokalizasyon
        if loc_hard and loc == "unknown":
            # PSORTb belirsiz — muhafazakâr; gerçek yüzey proteinini kaçırmamak için
            # ELEME YAPMA, dürüst 'belirsiz' notuyla geçir.
            pr.trace.append({"step": "funnel:localization", "ok": True,
                             "reason": "PSORTb 'Unknown' — lokalizasyon belirsiz, elenmedi (muhafazakârlık)",
                             "warning": True})
        elif loc_hard:
            ok_loc = loc in loc_allowed
            pr.mark("funnel:localization", ok_loc, f"{loc} {'∈' if ok_loc else '∉'} izinli", localization=loc)
        else:
            pr.annotations["method_localization"] = f"ÇALIŞTIRILMADI ({profile}: bakteri-dışı, DeepLoc gerekir)"
            pr.trace.append({"step": "funnel:localization", "ok": True,
                             "reason": f"ATLANDI — lokalizasyon heuristiği yalnız bakteride geçerli ({profile})",
                             "warning": True})
        # -- sert: TM sayısı
        ok_tm = tm <= max_tm
        pr.mark("funnel:transmembrane", ok_tm, f"{tm} TM ≤ {max_tm}" if ok_tm else f"{tm} TM > {max_tm}", tm=tm)
        # -- sert: self-homoloji (ZORUNLU otoimmünite güvenliği) — HOST-DRIVEN
        if hh_available:
            hd = hh_data.get(pr.id, {"max_pident": 0.0, "by_host": {}})
            pid = hd["max_pident"]
            # en yüksek benzerliği veren konağı bul (raporlama için)
            worst_host = max(hd["by_host"], key=hd["by_host"].get, default="—") if hd["by_host"] else "—"
            pr.annotations["self_homology"] = pid
            pr.annotations["self_homology_by_host"] = hd["by_host"]
            pr.annotations["human_homology"] = pid   # geriye-uyum (rapor/evaluate)
            pr.annotations["method_human_homology"] = f"GERÇEK (diamond → {', '.join(hh_scanned)} proteom)"
            ok_hh = pid < hh_id   # eşik ÜSTÜ konağa benzer -> elenir
            pr.mark("funnel:human_homology", ok_hh,
                    f"konak-öz %{pid} < {hh_id} (güvenli)" if ok_hh
                    else f"{worst_host} %{pid} ≥ {hh_id} → OTOİMMÜNİTE riski, elendi", human_pident=pid)
        else:
            pr.annotations["human_homology"] = "not_run"
            pr.trace.append({"step": "funnel:human_homology", "ok": True,
                             "reason": "ATLANDI — seçili konakların hiçbirinde proteom DB'si yok (gerçek kullanımda ZORUNLU)",
                             "warning": True})

        # -- yumuşak skorlama (0-1): sinyal + antijenite + adezin proxy
        sp_score = min(1.0, sp / max(sp_min, 1e-6)) * 0.5 if sp < sp_min else 0.5 + 0.5 * min(1.0, sp)
        funnel_score = round(0.5 * ag + 0.3 * min(1.0, sp_score) + 0.2 * (1 if ag >= vaxijen_thr else 0.3), 3)
        pr.annotations["funnel_score"] = funnel_score

        if pr.passed:
            survivors.append(pr)

    if psortb_ok:
        loc_method = f"GERÇEK (PSORTb, Gram{'+' if gram=='positive' else '−'})"
    else:
        loc_method = "heuristik (Kyte-Doolittle + sinyal)"
    hh_method = (f"GERÇEK (diamond → {', '.join(hh_scanned)})" if hh_available
                 else "ÇALIŞTIRILMADI (seçili konaklarda DB yok)")
    if hh_missing:
        hh_method += f" · DB'siz: {', '.join(hh_missing)}"
    summary = {
        "girdi": len(proteins), "hayatta": len(survivors),
        "lokalizasyon": loc_method + (" [sert:bakteri]" if loc_hard else f" [atlandı:{profile}]"),
        "transmembran": "GERÇEK (TMHMM-2.0)" if tm_available else "Kyte-Doolittle (klasik)",
        "sinyal_peptidi": "GERÇEK (SignalP-5.0)" if sp_available else "heuristik",
        "antijenite": "GERÇEK (IApred)" if ia_available else ("z-ACC model" if ag_available else "heuristik"),
        "self_homoloji": hh_method,
    }
    return survivors, summary
