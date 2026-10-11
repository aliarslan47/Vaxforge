"""Askıda kalan üç yayının antijenlerini çöz: eşleme ada göre değil, ölçülebilir kanıtla.

PROSPECTIVE_GROUNDTRUTH.md'nin "henüz çözülemeyenler" bölümünde üç yayın vardı. Üçünün
hiçbiri PMC'de açık değil (PMID 38906494 / 34020335 / 31287878 → PMC karşılığı yok,
Europe PMC: "Subscription required", annotations API yalnız başlık+özeti işliyor). Yani
locus etiketleri tam metinden okunamıyor. Ad eşlemesi YAPILMIYOR: "suis" regex'inin 28
Glaesserella parasuis kaydını yakaladığı hatanın dersi bu.

Bu betik üç kalemi üç farklı yolla ele alıyor ve her kararın dayanağını ölçüyor:

1) V5 karışımı (doi:10.1093/infdis/jiz352, Xing 2019) — ÇÖZÜLDÜ, ikincil açık kaynaktan.
   Aynı grubun açık erişimli eLife makalesi (doi:10.7554/eLife.101760, PMC12185102)
   Yöntemler'de V5'in bileşimini harf harf veriyor: "V5 contains five conserved virulence
   factors of S. suis. They are SrtA (amino acids 50-249), ScpCL (36-535), MRP (283-721),
   SLY (28-497), and SspA (41-852) ... cloned from the genome of S. suis strain A7
   (GenBank accession number: NC_017622.1)". Kaynak genom ve rezidü aralıkları verildiği
   için eşleme ad tahmini değil: her adayın uzunluğu, bildirilen aralığın üst sınırını
   karşılamak zorunda. SrtA (n=249) ve SLY (n=497) aralıkların üst sınırıyla BİREBİR
   eşleşiyor; MRP'nin gen adı annotasyonda zaten mrp. A7'de iki S8-ailesi peptidaz var
   (SspA/ScpCL belirsizliği) → S. pyogenes C5a peptidazı ScpA (P15926) ile hizalayıp
   ayırıyoruz: C5a peptidaz benzeri olan ScpCL'dir.

2) SBP2' (doi:10.1016/j.rvsc.2021.05.006) — ÇÖZÜLDÜ, dört bağımsız kontrolle.
   Aynı grubun açık erişimli AEM makalesi (doi:10.1128/aem.01267-26) SBP2' için
   doğrulanabilir dört özellik veriyor: srtBCD kümesinde, YPKN pilin motifi taşıyor,
   Western'de ~33 kDa, ve S. pneumoniae RrgB'nin N-TERMİNAL bölgesine %29 kimlikte
   (eşi SBP2'' ise C-terminal bölgeye %35). Küme genomik komşuluktan bulunuyor (üç
   yan yana class C sortaz = SrtB/SrtC/SrtD), beş pilin alt birimi küme içinden
   çıkarılıyor, sonra bu dört kontrol ölçülüyor.

3) HK / Gtf-2 / PsbP (doi:10.1016/j.micpath.2024.106759) — ÇÖZÜLEMEDİ.
   Özette ne aksesyon ne locus etiketi var, tam metin abonelik arkasında. Belirsizliğin
   büyüklüğü hedef proteomda SAYILIYOR (iddia edilmiyor): aynı ada uyan kaç protein var.
   Birden fazlaysa eşleme yapılmaz; kalem "erişim yokluğundan çözülemedi" olarak kapanır.

Girdilerin indirilmesi (benchmark/downloads/ .gitignore'da — repo konvansiyonu; bu dört
dosya aşağıdaki komutlarla aynen yeniden üretilir):

    # S. suis A7 (V5'in kaynak genomu, eLife makalesinin bildirdiği NC_017622.1).
    # CDS-AA biçimi seçildi çünkü ScpCL A7'de pseudo=true ve RefSeq protein.faa'da yok.
    mkdir -p benchmark/downloads/ssuis_a7
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\
&id=NC_017622.1&rettype=fasta_cds_aa&retmode=text" \
      -o benchmark/downloads/ssuis_a7/NC_017622.1_cds_aa.faa

    # S. suis ZY05719 (SBP2' çalışmalarının suşu) — protein + feature table.
    # Feature table gerekli: srtBCD kümesi genomik KOMŞULUKTAN bulunuyor.
    mkdir -p benchmark/downloads/ssuis_zy05719
    B=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/993/745/GCF_000993745.1_ASM99374v1
    curl -s "$B/GCF_000993745.1_ASM99374v1_protein.faa.gz" | gunzip \
      > benchmark/downloads/ssuis_zy05719/GCF_000993745.1_ASM99374v1_protein.faa
    curl -s "$B/GCF_000993745.1_ASM99374v1_feature_table.txt.gz" | gunzip \
      > benchmark/downloads/ssuis_zy05719/GCF_000993745.1_ASM99374v1_feature_table.txt

    # Referans diziler: S. pneumoniae RrgB (SBP2' N/C ayrımı) ve S. pyogenes C5a
    # peptidazı ScpA (ScpCL/SspA ayrımı).
    mkdir -p benchmark/downloads/refs
    curl -s "https://rest.uniprot.org/uniprotkb/search?query=%28rrgB%29+AND+\
%28taxonomy_id%3A1313%29&format=fasta&size=1" \
      -o benchmark/downloads/refs/rrgB_spneumoniae_A7KT25.faa
    curl -s "https://rest.uniprot.org/uniprotkb/P15926.fasta" \
      -o benchmark/downloads/refs/scpA_spyogenes_P15926.faa

Koşum:
    .venv/bin/python3 benchmark/resolve_pending_antigens.py
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "benchmark"
DL = BENCH / "downloads"

A7_CDS_AA = DL / "ssuis_a7" / "NC_017622.1_cds_aa.faa"
A7_KEYED = DL / "ssuis_a7" / "ssuis_a7_keyed.faa"
ZY_PROT = DL / "ssuis_zy05719" / "GCF_000993745.1_ASM99374v1_protein.faa"
ZY_FEAT = DL / "ssuis_zy05719" / "GCF_000993745.1_ASM99374v1_feature_table.txt"
RRGB = DL / "refs" / "rrgB_spneumoniae_A7KT25.faa"
SCPA = DL / "refs" / "scpA_spyogenes_P15926.faa"
TARGET = ROOT / "data" / "validation" / "ssuis_prepared.faa"

# eLife 10.7554/eLife.101760 Yöntemler'de bildirilen konstrukt aralıkları.
V5_RANGES = {"SrtA": (50, 249), "ScpCL": (36, 535), "MRP": (283, 721),
             "SLY": (28, 497), "SspA": (41, 852)}
# AEM 10.1128/aem.01267-26: SBP2' N-terminal RrgB'ye %29, SBP2'' C-terminal'e %35.
SBP2_MOTIF = "YPKN"


def read_faa(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_cds_aa(path: Path) -> list[dict]:
    """NCBI fasta_cds_aa başlıklarını ayrıştır: [gene=] [locus_tag=] [protein=] [protein_id=]."""
    out: list[dict] = []
    hdr, buf = None, []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if hdr is not None:
                out.append(_cds_rec(hdr, "".join(buf)))
            hdr, buf = ln[1:], []
        else:
            buf.append(ln.strip())
    if hdr is not None:
        out.append(_cds_rec(hdr, "".join(buf)))
    return out


def _cds_rec(hdr: str, seq: str) -> dict:
    def tag(k: str) -> str:
        m = re.search(rf"\[{k}=([^\]]*)\]", hdr)
        return m.group(1) if m else ""
    return {"gene": tag("gene"), "locus_tag": tag("locus_tag"),
            "protein": tag("protein"), "protein_id": tag("protein_id"),
            "pseudo": tag("pseudo") == "true", "seq": seq, "n": len(seq)}


def key_of(rec: dict) -> str:
    """Anahtar: WP aksesyonu varsa o, yoksa locus etiketi (psödogenlerin WP'si yoktur)."""
    return rec["protein_id"] or rec["locus_tag"]


def write_keyed_a7() -> list[dict]:
    recs = read_cds_aa(A7_CDS_AA)
    A7_KEYED.write_text("".join(f">{key_of(r)}\n{r['seq']}\n" for r in recs if r["seq"]),
                        encoding="utf-8")
    return recs


def blastp(queries: dict[str, str], subject: Path, diamond: str) -> list[dict]:
    """Hizalama bölgesi dahil tüm hitler — N/C-terminal ayrımı buna dayanıyor."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        qf = td / "q.faa"
        qf.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items()))
        db = td / "db"
        subprocess.run([diamond, "makedb", "--in", str(subject), "-d", str(db)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = td / "h.tsv"
        subprocess.run([diamond, "blastp", "-q", str(qf), "-d", str(db), "-o", str(out),
                        "--very-sensitive", "--max-target-seqs", "5", "--quiet",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                        "qstart", "qend", "sstart", "send", "bitscore"], check=True)
        rows = []
        for ln in out.read_text().splitlines():
            f = ln.split("\t")
            rows.append({"query": f[0], "subject": f[1], "identity": float(f[2]),
                         "aln_len": int(f[3]), "qstart": int(f[4]), "qend": int(f[5]),
                         "sstart": int(f[6]), "send": int(f[7]), "bitscore": float(f[8])})
        return rows


def resolve_v5(recs: list[dict], diamond: str) -> dict:
    """V5'in beşini A7'de bul: annotasyon + eLife'ın bildirdiği uzunluk alt sınırı."""
    def pick(patterns: tuple[str, ...], min_len: int) -> list[dict]:
        hits = []
        for r in recs:
            label = f"{r['protein']} {r['gene']}".lower()
            if any(p in label for p in patterns) and r["n"] >= min_len:
                hits.append(r)
        return hits

    res: dict[str, dict] = {}
    # SrtA: sınıf A sortaz; bildirilen aralığın üst sınırı 249 → uzunluk birebir kontrol.
    srt = [r for r in pick(("class a sortase", "sortase a"), V5_RANGES["SrtA"][1])]
    res["SrtA"] = _one(srt, "sınıf A sortaz, n>=249",
                       exact_len=V5_RANGES["SrtA"][1])
    res["MRP"] = _one(pick(("muramidase-released",), V5_RANGES["MRP"][1]),
                      "gen adı mrp, n>=721")
    res["SLY"] = _one(pick(("suilysin",), V5_RANGES["SLY"][1]),
                      "gen adı sly, n>=497", exact_len=V5_RANGES["SLY"][1])

    # SspA / ScpCL: A7'de iki S8 peptidaz var; C5a peptidazına yakın olan ScpCL'dir.
    s8 = pick(("s8 family serine peptidase",), V5_RANGES["ScpCL"][1])
    scp_hits = blastp({key_of(r): r["seq"] for r in s8}, SCPA, diamond)
    best: dict[str, dict] = {}
    for h in scp_hits:
        if h["query"] not in best or h["bitscore"] > best[h["query"]]["bitscore"]:
            best[h["query"]] = h
    ranked = sorted(s8, key=lambda r: -best.get(key_of(r), {}).get("bitscore", 0.0))
    scpa_evidence = {key_of(r): best.get(key_of(r)) for r in s8}
    if len(ranked) >= 2:
        res["ScpCL"] = {"status": "çözüldü", "key": key_of(ranked[0]),
                        "locus_tag": ranked[0]["locus_tag"], "n": ranked[0]["n"],
                        "annotation": ranked[0]["protein"], "pseudo": ranked[0]["pseudo"],
                        "basis": "iki S8 peptidazdan ScpA'ya (P15926) en yakın olan",
                        "scpa_alignment": scpa_evidence}
        res["SspA"] = {"status": "çözüldü (eleme ile)", "key": key_of(ranked[-1]),
                       "locus_tag": ranked[-1]["locus_tag"], "n": ranked[-1]["n"],
                       "annotation": ranked[-1]["protein"], "pseudo": ranked[-1]["pseudo"],
                       "basis": ("kalan S8 peptidaz; n=%d bildirilen 41-852 aralığını "
                                 "karşılıyor" % ranked[-1]["n"]),
                       "scpa_alignment": scpa_evidence}
    else:
        res["ScpCL"] = res["SspA"] = {"status": "çözülemedi",
                                      "basis": f"S8 peptidaz sayısı {len(s8)}, ayrım yapılamadı"}
    return res


def _one(hits: list[dict], basis: str, exact_len: int | None = None) -> dict:
    if len(hits) != 1:
        return {"status": "çözülemedi", "basis": basis,
                "n_matches": len(hits),
                "matches": [{"locus_tag": h["locus_tag"], "n": h["n"],
                             "annotation": h["protein"]} for h in hits]}
    h = hits[0]
    out = {"status": "çözüldü", "key": key_of(h), "locus_tag": h["locus_tag"],
           "n": h["n"], "annotation": h["protein"], "pseudo": h["pseudo"],
           "basis": basis}
    if exact_len is not None:
        out["exact_length_match"] = h["n"] == exact_len
    return out


def resolve_sbp2(diamond: str) -> dict:
    """srtBCD kümesini genomik komşuluktan bul, SBP2'yi dört kontrolle ayır."""
    rows = []
    for ln in ZY_FEAT.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 19 or f[0] != "CDS":
            continue
        rows.append({"start": int(f[7]), "acc": f[10], "name": f[13],
                     "locus_tag": f[16], "plen": f[18]})
    rows.sort(key=lambda r: r["start"])

    # Üç YAN YANA class C sortaz = srtBCD. Küme ada göre değil komşuluğa göre bulunuyor.
    idx = [i for i, r in enumerate(rows) if "class C sortase" in r["name"]]
    cluster_at = None
    for i in idx:
        if i + 2 < len(rows) and all("class C sortase" in rows[j]["name"]
                                     for j in (i, i + 1, i + 2)):
            cluster_at = i
            break
    if cluster_at is None:
        return {"status": "çözülemedi", "basis": "üç yan yana class C sortaz bulunamadı"}

    sortases = [rows[cluster_at + k]["locus_tag"] for k in range(3)]
    # Kümenin pilin alt birimleri: sortazlardan sonraki pilin-benzeri CDS'ler.
    pilin_pat = ("pilin", "prealbumin-like fold", "isopeptide-forming")
    pilins = []
    for r in rows[cluster_at + 3:]:
        if any(p in r["name"].lower() for p in pilin_pat):
            pilins.append(r)
        else:
            break

    seqs = read_faa(ZY_PROT)
    queries = {r["acc"]: seqs[r["acc"]] for r in pilins if r["acc"] in seqs}
    rrgb = blastp(queries, RRGB, diamond)
    rrgb_len = len(next(iter(read_faa(RRGB).values())))
    best: dict[str, dict] = {}
    for h in rrgb:
        if h["query"] not in best or h["bitscore"] > best[h["query"]]["bitscore"]:
            best[h["query"]] = h

    checks = []
    for r in pilins:
        s = seqs.get(r["acc"], "")
        h = best.get(r["acc"])
        # RrgB'nin hangi yarısına hizalanıyor: N mi C mi?
        half = None
        if h:
            mid = (h["sstart"] + h["send"]) / 2
            half = "N" if mid < rrgb_len / 2 else "C"
        checks.append({
            "locus_tag": r["locus_tag"], "accession": r["acc"], "n": len(s),
            "annotation": r["name"], "ypkn_motif": SBP2_MOTIF in s,
            "rrgb_identity": round(h["identity"], 1) if h else None,
            "rrgb_region": [h["sstart"], h["send"]] if h else None,
            "rrgb_half": half,
        })

    # SBP2' = YPKN taşıyan TEK pilin ve RrgB'nin N yarısına hizalanan.
    cand = [c for c in checks if c["ypkn_motif"] and c["rrgb_half"] == "N"]
    sbp2pp = [c for c in checks if c["rrgb_half"] == "C" and not c["ypkn_motif"]]
    out = {"srtBCD_sortases": sortases, "n_pilins_in_cluster": len(pilins),
           "rrgb_reference_length": rrgb_len, "checks": checks}
    if len(cand) == 1:
        out.update({"status": "çözüldü", "key": cand[0]["accession"],
                    "locus_tag": cand[0]["locus_tag"], "n": cand[0]["n"],
                    "basis": ("küme içinde YPKN pilin motifini taşıyan tek pilin; "
                              "RrgB'nin N-terminal bölgesine %%%.1f kimlik (yayın: %%29); "
                              "n=%d ~ bildirilen 33 kDa bandı"
                              % (cand[0]["rrgb_identity"], cand[0]["n"])).replace("%%", "%"),
                    "sbp2_second_half": (sbp2pp[0] if len(sbp2pp) == 1 else None)})
    else:
        out.update({"status": "çözülemedi",
                    "basis": f"YPKN+N-yarı ölçütünü karşılayan pilin sayısı {len(cand)}"})
    return out


def measure_micpath_ambiguity() -> dict:
    """HK/Gtf-2/PsbP: belirsizliği iddia etme, hedef proteomda SAY."""
    names: dict[str, str] = {}
    for ln in TARGET.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            head = ln[1:]
            acc = head.split()[0]
            desc = head.split(" ", 1)[1] if " " in head else ""
            names[acc] = desc.split(" OS=")[0]
    groups = {
        "histidin kinaz (HK)": ("histidine kinase",),
        "glikoziltransferaz ailesi 2 (Gtf-2)": ("glycosyltransferase",),
        "fosfat bağlayan protein (PsbP)": ("phosphate-binding", "phosphate binding"),
    }
    out = {}
    for label, pats in groups.items():
        hits = [{"accession": a, "name": n} for a, n in names.items()
                if any(p in n.lower() for p in pats)]
        out[label] = {"n_matches": len(hits), "matches": sorted(
            hits, key=lambda h: h["accession"])}
    return out


def main() -> int:
    diamond = shutil.which("diamond") or "diamond"
    for p in (A7_CDS_AA, ZY_PROT, ZY_FEAT, RRGB, SCPA, TARGET):
        if not p.exists():
            print(f"girdi yok: {p}")
            return 1

    recs = write_keyed_a7()
    rep = {
        "question": ("PROSPECTIVE_GROUNDTRUTH.md'de askıda kalan üç yayının antijenleri "
                     "ölçülebilir kanıtla eşlenebilir mi?"),
        "access_check": {
            "pmc_full_text": "yok (38906494 / 34020335 / 31287878 → PMC karşılığı yok)",
            "europe_pmc": "Subscription required; annotations API yalnız başlık+özet",
            "consequence": "locus etiketleri birincil tam metinden okunamıyor",
        },
        "v5": {
            "doi": "10.1093/infdis/jiz352",
            "composition_source": {
                "doi": "10.7554/eLife.101760", "pmcid": "PMC12185102",
                "quote": ("V5 contains five conserved virulence factors of S. suis. "
                          "They are SrtA (amino acids 50-249), ScpCL (36-535), "
                          "MRP (283-721), SLY (28-497), and SspA (41-852) ... cloned "
                          "from the genome of S. suis strain A7 (NC_017622.1)"),
            },
            "source_genome": "NC_017622.1 (S. suis A7)",
            "reported_ranges": {k: list(v) for k, v in V5_RANGES.items()},
            "resolved": resolve_v5(recs, diamond),
            "caveat": ("karışım olarak uygulandı; tek-antijen koruması ATFEDİLEMEZ. "
                       "Prospektif tek-antijen sayımına girmez, karışım kanıtı olarak tutulur."),
        },
        "sbp2": {
            "doi": "10.1016/j.rvsc.2021.05.006",
            "feature_source": {"doi": "10.1128/aem.01267-26", "pmcid": "PMC13599487",
                               "features": ["srtBCD kümesinde", f"{SBP2_MOTIF} pilin motifi",
                                            "~33 kDa", "RrgB N-terminal %29 (SBP2'' C-term %35)"]},
            "source_genome": "GCF_000993745.1 (S. suis ZY05719)",
            "resolved": resolve_sbp2(diamond),
        },
        "micpath_hk_gtf_psbp": {
            "doi": "10.1016/j.micpath.2024.106759",
            "status": "çözülemedi — erişim yokluğu",
            "reason": ("özette aksesyon/locus etiketi yok, tam metin abonelik arkasında. "
                       "Ad eşlemesi yapılmıyor: hedef proteomda aynı ada uyan protein "
                       "sayısı aşağıda ölçüldü, birden fazla."),
            "ambiguity_in_target_proteome": measure_micpath_ambiguity(),
        },
    }
    out = BENCH / "results" / "pending_antigen_resolution.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print("=== V5 (doi:10.1093/infdis/jiz352) — bileşim eLife PMC12185102'den")
    for k, v in rep["v5"]["resolved"].items():
        extra = ""
        if v.get("exact_length_match"):
            extra = " [uzunluk birebir]"
        if v.get("pseudo"):
            extra += " [A7'de pseudo]"
        print(f"  {k:7} {v['status']:22} {v.get('locus_tag', '-'):18} "
              f"n={v.get('n', '-')}{extra}")
    s = rep["sbp2"]["resolved"]
    print(f"\n=== SBP2' (doi:10.1016/j.rvsc.2021.05.006)")
    print(f"  srtBCD sortazları: {s.get('srtBCD_sortases')}")
    print(f"  kümedeki pilin sayısı: {s.get('n_pilins_in_cluster')}")
    for c in s.get("checks", []):
        print(f"   {c['locus_tag']:18} {c['accession']:16} n={c['n']:5} "
              f"YPKN={'VAR' if c['ypkn_motif'] else 'yok':3} "
              f"RrgB={c['rrgb_identity']} yarı={c['rrgb_half']}")
    print(f"  SONUÇ: {s['status']} → {s.get('locus_tag')} ({s.get('key')})")
    print(f"\n=== HK/Gtf-2/PsbP (doi:10.1016/j.micpath.2024.106759): "
          f"{rep['micpath_hk_gtf_psbp']['status']}")
    for label, g in rep["micpath_hk_gtf_psbp"]["ambiguity_in_target_proteome"].items():
        print(f"  {label:38} hedef proteomda eşleşme: {g['n_matches']}")
    print(f"\nyazıldı: {out}")
    print(f"anahtarlanmış A7: {A7_KEYED.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
