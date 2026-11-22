"""Adım 4 — B/T-hücre epitop tahmini (sliding-window) + ÇOK-KONAKLI MHC.

B-hücre: Parker hidrofilisite kayan-pencere (gerçek: BepiPred-3.0).
MHC-I ve MHC-II: SEÇİLEN KONAK(LAR)ın allelleri taranır (predictors katmanı;
insana kilitli değil). Her peptit için hangi konaklarda sunulduğu (host_coverage)
hesaplanır — sizin 'hangi organizmada işe yarar' haritanız.

Yordayıcı gerçekliği predictors.py'de: yerel NetMHCpan/NetMHCIIpan (gerçek),
yoksa IEDB API, o da yoksa proxy (dürüst etiketli).
"""

from __future__ import annotations

from . import (bepipred, bepipred3, complexity, ifngamma, immunogenicity,
               mhc_motif, predictors, processing)
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


def _host_label(hosts: list[Host], name: str) -> str:
    """Konak adını okunur TÜR etiketine çevirir ('bovine' -> 'Sığır (Bos taurus, BoLA)').

    Anchor pozisyonları allele, alleller de konağa özgü olduğundan rapor
    tablosunun altında hangi TÜRE göre bakıldığı yazılmak zorundadır.
    """
    for h in hosts or []:
        if h.name == name:
            return getattr(h, "label", None) or name
    return name or "?"


def _mhc(pr: ProteinRecord, hosts: list[Host], mhc_class: str, kind_label: str,
         lengths: list[int], rank_strong: float, rank_weak: float,
         per_host: dict[str, dict[str, dict]], top: int = 6,
         proc: dict[int, dict[str, float]] | None = None,
         proc_backend=None) -> list[Peptide]:
    s = sanitize(pr.seq)
    pos = {}
    for L in lengths:
        for i in range(len(s) - L + 1):
            pos.setdefault(s[i:i + L], i)
    peps = list(pos)
    if not peps:
        return []
    # per_host: TÜM proteinler için ÖNCEDEN hesaplanmış GLOBAL tahmin (protein
    # başına ayrı netMHCpan çağrısı YERİNE tek batch — bkz. run()). Bu proteinin
    # peptitleri global sözlükten okunur.

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
            # MHC yarığı anchor motifi (SALT YORUMLAMA — sıralamayı etkilemez).
            # Anchor pozisyonları ALLELE ÖZGÜDÜR; hangi türün hangi allelinde
            # ölçüldüğü anchor_basis alanına yazılır ve raporda tablonun altına
            # basılır (2026-11-22).
            mhc_motif.annotate(p, best_allele, "mhc_i", host=_host_label(hosts, best_host))
            # Antijen işleme: peptidin C-terminaline denk gelen 9-mer'in kesim+TAP'ı
            # Normalizasyon omurganın KENDİSİNDEN gelir: bir omurganın cle'si 0-1
            # olasılık, IEDB'nin proteasome_score'u log-ölçekli — sabit bir
            # modülden normalize etmek ikinci durumda her şeyi 1.0'a doyurur.
            if proc and proc_backend is not None:
                ct = proc.get(pos[pep] + len(pep) - 1)
                if ct:
                    p.metrics["cleavage"] = ct["cle"]
                    p.metrics["tap"] = ct["tap"]
                    p.metrics["processing_norm"] = proc_backend.processing_norm(
                        ct["cle"], ct["tap"])
                    p.methods["processing"] = proc_backend.METHOD
        else:
            # MHC-II: anchor pozisyonları da allele özgüdür (P1/P4/P6/P9 çekirdek
            # cepleri ama çekirdek ofseti allele ve peptide göre kayar). Ölçüm
            # NetMHCIIpan varsa yapılır; yoksa annotate() hiçbir şey yazmaz
            # (varsayılan pozisyon UYDURULMAZ).
            mhc_motif.annotate(p, best_allele, "mhc_ii", host=_host_label(hosts, best_host))
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
    # batch. Omurga IEDB processing API (bkz. processing.backend) — hangisi
    # gerçekten koşuyorsa (bkz. processing.backend). Allel-bağımsız cle/tap
    # protein-başına C-terminal pozisyona göre haritalanır; ikisi de yoksa boş
    # (MHC-I işleme metrikleri dürüstçe ölçülmemiş kalır).
    proc_backend = processing.backend()
    proc_scores = (proc_backend.predict([(pr.id, sanitize(pr.seq)) for pr in proteins])
                   if proc_backend is not None else {})

    # --- MHC bağlanma: TÜM proteinlerin peptitleri TEK GLOBAL BATCH ---
    # Eskiden _mhc her protein için ayrı netMHCpan çağırıyordu (580 protein =
    # 580+ süreç başlatma → seri darboğaz). netMHCpan zaten tek çağrıda çok peptit
    # işler; bu yüzden peptitleri bir kez toplayıp TEK batch tahmin ediyoruz,
    # sonra protein bazında sözlükten okuyoruz. Sonuç: birkaç süreç, dakikalar.
    i_lengths = list(i_tool.params["peptide_lengths"].value)
    ii_length = int(ii_tool.params["peptide_length"].value)
    i_rank_weak = float(i_tool.params["rank_weak"].value)
    ii_rank_weak = float(ii_tool.params["rank_weak"].value)

    def _all_peptides(lengths: list[int]) -> list[str]:
        allp: set[str] = set()
        for pr in proteins:
            s = sanitize(pr.seq)
            for L in lengths:
                for i in range(len(s) - L + 1):
                    allp.add(s[i:i + L])
        return sorted(allp)

    all_i = _all_peptides(i_lengths)
    all_ii = _all_peptides([ii_length])
    # Güçlü eşik de omurgaya GEÇİRİLİR: `score` artık %rank'tan sınıfın kendi
    # eşiklerine demirli log-ölçekle türetiliyor (bkz. mhc_rank).
    i_rank_strong = float(i_tool.params["rank_strong"].value)
    ii_rank_strong = float(ii_tool.params["rank_strong"].value)
    pred_i = {h.name: predictors.predict(all_i, h, "mhc_i", i_rank_weak,
                                         rank_strong=i_rank_strong) for h in hosts}
    pred_ii = {h.name: predictors.predict(all_ii, h, "mhc_ii", ii_rank_weak,
                                          rank_strong=ii_rank_strong) for h in hosts}

    for pr in proteins:
        made = _bcell(pr, tools["bcell_epitope"], bp3=bp3_scores.get(pr.id))
        made += _mhc(pr, hosts, "mhc_i", "MHC-I", i_lengths,
                     float(i_tool.params["rank_strong"].value), i_rank_weak,
                     pred_i, proc=proc_scores.get(pr.id),
                     proc_backend=proc_backend)
        made += _mhc(pr, hosts, "mhc_ii", "MHC-II", [ii_length],
                     float(ii_tool.params["rank_strong"].value), ii_rank_weak,
                     pred_ii)
        parent_ag = pr.annotations.get("antigenicity")
        for p in made:
            # Dizi karmaşıklığı: homopolimerik/düşük karmaşıklıklı peptit
            # bağlanma eşiğini geçse bile aşı adayı OLAMAZ (gerekçe + kalibrasyon:
            # complexity.py). `passed`'a dokunmaz — bağlanma kararı ayrı ölçüt.
            complexity.annotate(p)
            p.metrics["parent_antigenicity"] = parent_ag
            p.metrics["gene"] = pr.annotations.get("gene")
            p.metrics["locus_tag"] = pr.annotations.get("locus_tag")
            p.metrics["location"] = pr.annotations.get("location")
            p.metrics["protein_id"] = pr.annotations.get("protein_id")
            # kaynak proteinin funnel (araç) sonuçlarını adaya taşı (rapor için)
            for k in ("localization", "tm_helices", "signalp", "surface_glyco", "human_homology",
                      "antigenicity_raw", "antigenicity_category", "virulence",
                      "vf_identity", "vf_keyword", "vf_hit",
                      "conservation_percent", "conservation_entropy_mean",
                      "conservation_n_orthologs",
                      "method_localization", "method_tm", "method_signalp",
                      "method_human_homology", "method_antigenicity", "method_discovery",
                      "method_conservation"):
                if k in pr.annotations:
                    p.metrics[k] = pr.annotations[k]
            # Epitop-penceresi konservasyonu: kaynak profilinin [start, start+len]
            # ortalaması (0-1). Skorlama bileşeni + rapor. Profil yoksa (suş yok) None.
            prof = pr.annotations.get("conservation_profile")
            if prof:
                a, b = p.start, min(len(prof), p.start + len(p.seq))
                if b > a:
                    win = prof[a:b]
                    p.metrics["epitope_conservation"] = round(sum(win) / len(win), 3)
        peptides += made

    # IFN-γ (HTL/Th1): MHC-II peptitleri için TEK batch (model bir kez yüklenir).
    # Yalnız insan/fare modeli var → seçili konaklardan uygun olanı kullanılır;
    # yoksa (ör. sığır/domuz) atlanır (dürüstçe metrik eklenmez).
    ifn_host = ifngamma.model_host([h.name for h in hosts])
    ii_peps = [p for p in peptides if p.kind == "MHC-II"]
    if ifn_host and ii_peps:
        ifn_tool = tools.get("ifn_gamma")
        thr = (float(ifn_tool.params["threshold"].value) if ifn_tool
               else ifngamma.DEFAULT_THRESHOLD)
        res = ifngamma.predict([p.seq for p in ii_peps], ifn_host, threshold=thr)
        for p in ii_peps:
            d = res.get(p.seq)
            if not d:
                continue
            p.metrics["ifn_gamma"] = d["ml"]
            p.metrics["ifn_gamma_blast"] = d["blast"]
            p.metrics["ifn_gamma_total"] = d["total"]
            p.metrics["ifn_gamma_inducer"] = d["inducer"]
            p.metrics["ifn_gamma_norm"] = d["norm"]
            p.metrics["ifn_gamma_host"] = ifn_host
            p.methods["ifn_gamma"] = ifngamma.METHOD
    return peptides
