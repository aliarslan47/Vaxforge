"""Adım 4 — B/T-hücre epitop tahmini (sliding-window) + ÇOK-KONAKLI MHC.

B-hücre: Parker hidrofilisite kayan-pencere (gerçek: BepiPred-3.0).
MHC-I ve MHC-II: SEÇİLEN KONAK(LAR)ın allelleri taranır (predictors katmanı;
insana kilitli değil). Her peptit için hangi konaklarda sunulduğu (host_coverage)
hesaplanır — sizin 'hangi organizmada işe yarar' haritanız.

Yordayıcı gerçekliği predictors.py'de: MHCflurry (insan MHC-I gerçek), diğerleri
NetMHCpan/NetMHCIIpan kurulunca gerçek, yoksa proxy (dürüst etiketli).
"""

from __future__ import annotations

from . import bepipred, bepipred3, immunogenicity, netctl, predictors
from .config_loader import ResolvedTool
from .hosts import Host
from .models import Peptide, ProteinRecord
from .sequtils import kolaskar_mean, parker_mean, sanitize


def _bcell(pr: ProteinRecord, tool: ResolvedTool, top: int = 5,
           bp3: list[float] | None = None) -> list[Peptide]:
    """Lineer B-hücre epitopları.

    Öncelik: GERÇEK BepiPred-3.0 (ESM-2, batch — `bp3` per-residue skorları) >
    GERÇEK BepiPred-1.0 (yerel, per-residue) > GERÇEK klasik yöntemler
    (Kolaskar-Tongaonkar antijenite + Parker hidrofilisite). Hepsi
    yayınlanmış/atıflanabilir — proxy değil.
    """
    s = sanitize(pr.seq)
    win = int(tool.params["min_length"].value)
    if len(s) < win:
        return []

    # Per-residue skor kaynağı seçimi (dürüst 'method' etiketiyle).
    if bp3 is not None and len(bp3) == len(s):
        bp, use_bp, thr = bp3, True, bepipred3.DEFAULT_THRESHOLD
        method = "GERÇEK (BepiPred-3.0, ESM-2)"
    else:
        bp = bepipred.predict_residues(s)
        use_bp = len(bp) == len(s)
        thr = bepipred.DEFAULT_THRESHOLD
        method = ("GERÇEK (BepiPred-1.0, yerel)" if use_bp
                  else "GERÇEK (Kolaskar-Tongaonkar + Parker, klasik yöntemler)")

    scored = []
    for i in range(len(s) - win + 1):
        pep = s[i:i + win]
        kt = kolaskar_mean(pep)
        pk = parker_mean(pep)
        if use_bp:
            bpm = sum(bp[i:i + win]) / win           # BepiPred pencere ortalaması
            # BepiPred-3.0 skoru zaten 0-1 olasılık; 1.0 ham skoru sigmoide sokulur.
            score = round(bpm, 3) if bp3 is not None else round(1 / (1 + pow(2.718, -bpm)), 3)
            passed = bpm >= thr
            scored.append((score, i, pep, kt, pk, round(bpm, 3), passed))
        else:
            score = round(0.6 * min(1.0, max(0.0, (kt - 0.9) / 0.5))
                          + 0.4 * (1 / (1 + pow(2.718, -pk / 3))), 3)
            passed = kt >= max(1.0, kolaskar_mean(s))
            scored.append((score, i, pep, kt, pk, None, passed))
    scored.sort(reverse=True, key=lambda x: x[0])

    out = []
    for score, i, pep, kt, pk, bpm, passed in scored[:top]:
        p = Peptide(seq=pep, parent=pr.id, kind="B", start=i)
        p.metrics["bcell_score"] = score
        p.metrics["kolaskar"] = round(kt, 3)
        p.metrics["parker"] = round(pk, 2)
        if bpm is not None:
            p.metrics["bepipred"] = bpm
        p.methods["bcell_score"] = method
        p.passed = passed
        out.append(p)
    return out


def _mhc(pr: ProteinRecord, hosts: list[Host], mhc_class: str, kind_label: str,
         lengths: list[int], rank_strong: float, rank_weak: float, top: int = 6,
         proc: dict[int, dict[str, float]] | None = None) -> list[Peptide]:
    s = sanitize(pr.seq)
    pos = {}
    for L in lengths:
        for i in range(len(s) - L + 1):
            pos.setdefault(s[i:i + L], i)
    peps = list(pos)
    if not peps:
        return []
    per_host = {h.name: predictors.predict(peps, h, mhc_class, rank_weak) for h in hosts}

    ranked = []
    for pep in peps:
        best_score, best_rank, best_allele, best_host = 0.0, 99.0, "", ""
        presented, detail, methods = [], {}, set()
        for h in hosts:
            d = per_host[h.name].get(pep)
            if not d:
                continue
            detail[h.name] = {"rank": d["rank"], "n_alleles": d["n_alleles"], "allele": d["best_allele"],
                              "bound_alleles": d.get("bound_alleles", [])}
            methods.add(d["method"])
            if d["rank"] <= rank_weak:
                presented.append(h.name)
            if d["score"] > best_score:
                best_score = d["score"]
            if d["rank"] < best_rank:
                best_rank, best_allele, best_host = d["rank"], d["best_allele"], h.name
        cov = len(presented) / max(1, len(hosts))
        ranked.append((best_score, pep, best_rank, best_allele, best_host, presented, cov, detail, methods))

    ranked.sort(reverse=True)
    out = []
    for best_score, pep, best_rank, best_allele, best_host, presented, cov, detail, methods in ranked:
        p = Peptide(seq=pep, parent=pr.id, kind=kind_label, start=pos[pep])
        p.metrics.update({
            "mhc_score": best_score, "pseudo_rank": best_rank,
            "best_allele": best_allele, "best_host": best_host,
            "hosts_presented": presented, "host_coverage": len(presented),
            "coverage_frac": cov, "per_host": detail,
        })
        p.methods["mhc_score"] = " | ".join(sorted(methods)) if methods else "—"
        # MHC-I: bağlanmadan SONRA CD8+ T-hücre immünojenitesi (en iyi allel maskesiyle)
        if mhc_class == "mhc_i":
            raw = immunogenicity.score(pep, best_allele)
            p.metrics["immunogenicity"] = raw
            p.metrics["immunogenicity_norm"] = immunogenicity.score_norm(raw)
            p.methods["immunogenicity"] = immunogenicity.METHOD
            # Antijen işleme: peptidin C-terminaline denk gelen 9-mer'in kesim+TAP'ı
            if proc:
                ct = proc.get(pos[pep] + len(pep) - 1)
                if ct:
                    p.metrics["cleavage"] = ct["cle"]
                    p.metrics["tap"] = ct["tap"]
                    p.metrics["processing_norm"] = netctl.processing_norm(ct["cle"], ct["tap"])
                    p.methods["processing"] = netctl.METHOD
        p.passed = best_rank <= rank_weak
        p.notes.append("güçlü" if best_rank <= rank_strong else "zayıf")
        out.append(p)
        if len(out) >= top:
            break
    return out


def run(proteins: list[ProteinRecord], tools: dict[str, ResolvedTool],
        hosts: list[Host]) -> list[Peptide]:
    peptides: list[Peptide] = []
    i_tool, ii_tool = tools["mhc_class_i"], tools["mhc_class_ii"]
    # B-hücre: GERÇEK BepiPred-3.0 (ESM-2) TÜM proteinler için TEK batch (ESM-2
    # bir kez yüklenir). Kuruluysa gerçek per-residue skor; yoksa _bcell içinde
    # BepiPred-1.0 / klasik yönteme düşer.
    bp3_scores = bepipred3.predict_batch([(pr.id, sanitize(pr.seq)) for pr in proteins])
    # Antijen işleme (proteozomal kesim + TAP): funnel-survivor proteinleri TEK
    # batch NetCTL. Allel-bağımsız cle/tap protein-başına C-terminal pozisyona
    # göre haritalanır; kuruluysa gerçek, yoksa boş (MHC-I metrikleri atlanır).
    proc_scores = netctl.predict([(pr.id, sanitize(pr.seq)) for pr in proteins])
    for pr in proteins:
        made = _bcell(pr, tools["bcell_epitope"], bp3=bp3_scores.get(pr.id))
        made += _mhc(pr, hosts, "mhc_i", "MHC-I",
                     list(i_tool.params["peptide_lengths"].value),
                     float(i_tool.params["rank_strong"].value),
                     float(i_tool.params["rank_weak"].value),
                     proc=proc_scores.get(pr.id))
        made += _mhc(pr, hosts, "mhc_ii", "MHC-II",
                     [int(ii_tool.params["peptide_length"].value)],
                     float(ii_tool.params["rank_strong"].value),
                     float(ii_tool.params["rank_weak"].value))
        parent_ag = pr.annotations.get("antigenicity")
        for p in made:
            p.metrics["parent_antigenicity"] = parent_ag
            p.metrics["gene"] = pr.annotations.get("gene")
            p.metrics["locus_tag"] = pr.annotations.get("locus_tag")
            p.metrics["location"] = pr.annotations.get("location")
            p.metrics["protein_id"] = pr.annotations.get("protein_id")
        peptides += made
    return peptides
