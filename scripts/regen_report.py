"""Bir koşunun çıktılarından (candidates_full.xlsx + run.json) Peptide nesnelerini
yeniden kurup raporu GÜNCEL kodla yeniden üretir — pipeline'ı yeniden koşmadan.
Kullanım: python scripts/regen_report.py <run_dir> [tr|en]
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from vaxforge import report
from vaxforge.models import Peptide

run_dir = Path(sys.argv[1])
lang = sys.argv[2] if len(sys.argv) > 2 else "tr"
meta = json.loads((run_dir / "run.json").read_text())
meta["lang"] = lang

xl = run_dir / "candidates_full.xlsx"
adaylar = pd.read_excel(xl, "Adaylar")
tools = pd.read_excel(xl, "Arac_sonuclari")


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _num_prefix(s):
    m = re.match(r"\s*(-?\d+\.?\d*)", str(s))
    return float(m.group(1)) if m else None


# peptit -> {araç: sonuç, ...} (+ durum)
tmap, tstatus = {}, {}
for _, r in tools.iterrows():
    tmap.setdefault(r["peptit"], {})[str(r["araç"])] = str(r["sonuç"])
    tstatus.setdefault(r["peptit"], {})[str(r["araç"])] = str(r.get("durum", ""))


def _find(d, key):
    for k, v in d.items():
        if key in k:
            return v
    return None


def _parse_kv(s, sep="="):
    out = {}
    for part in str(s).split(";"):
        if sep in part:
            k, v = part.split(sep, 1)
            out[k.strip()] = v.strip()
    return out


def _parse_motif(s):
    # "HLA-B*08:01∈{L,M,I}" -> {"HLA-B*08:01": ["L","M","I"]}
    out = {}
    for part in str(s).split(";"):
        m = re.match(r"(.+?)∈\{(.+?)\}", part.strip())
        if m:
            out[m.group(1)] = [x.strip() for x in m.group(2).split(",") if x.strip()]
    return out


# MHC güçlü eşiği (notes için)
thr = {(t.get("tool"), t.get("param")): t.get("value") for t in meta.get("thresholds", [])}
rs_i = _f(thr.get(("mhc_class_i", "rank_strong"))) or 0.5
rs_ii = _f(thr.get(("mhc_class_ii", "rank_strong"))) or 2.0

peptides = []
for _, r in adaylar.iterrows():
    p = Peptide(seq=str(r["peptide"]), parent=str(r["cds_source"]),
                kind=str(r["kind"]), start=0)
    p.candidacy = _f(r["candidacy"])
    p.passed = bool(r["passed"])
    td = tmap.get(r["peptide"], {})
    ts = tstatus.get(r["peptide"], {})
    ag = _num_prefix(_find(td, "Antijenite"))
    rank = _f(r["best_rank"])
    strong = rank is not None and rank <= (rs_i if r["kind"] == "MHC-I" else rs_ii)
    p.notes = ["güçlü"] if strong else ["zayıf"]
    ifn_row = _find(td, "IFN-γ")
    ifn_ok = ifn_row is not None and "PASS" in (_find(ts, "IFN-γ") or "").upper()
    m = {
        "parent_antigenicity": ag, "pseudo_rank": rank,
        "best_allele": (None if pd.isna(r["best_allele"]) else str(r["best_allele"])),
        "allergen": bool(r["allergen"]) if not pd.isna(r["allergen"]) else False,
        "toxicity": _f(r["toxicity"]), "immunogenicity": _f(r["immunogenicity"]),
        "processing_norm": _f(r["processing_norm"]), "bcell_score": _f(r["bcell_score"]),
        "ifn_gamma_inducer": ifn_ok, "host_coverage": _f(r["host_coverage"]),
        "gene": (None if pd.isna(r["gene"]) else str(r["gene"])),
        "locus_tag": (None if pd.isna(r["locus_tag"]) else str(r["locus_tag"])),
        "location": (None if pd.isna(r["location"]) else str(r["location"])),
        "mhc_score": _f(r["mhc_score"]),
    }
    if not pd.isna(r.get("anchor_residues")) and str(r["anchor_residues"]).strip():
        m["anchor_residues"] = _parse_kv(r["anchor_residues"])
        m["allele_anchor_motif"] = _parse_motif(r["allele_anchor_motif"])
        m["anchor_match"] = (None if pd.isna(r["anchor_match"]) else str(r["anchor_match"]))
    # IEDB
    iem = str(r.get("iedb_match", "")).strip()
    matched = bool(iem) and iem.lower() not in ("nan", "false", "no", "yok", "-", "")
    if matched:
        pmids = [x.strip() for x in str(r.get("iedb_pmids", "")).split(",") if x.strip() and x.strip().lower() != "nan"]
        m["iedb"] = {"matched": True, "match_type": iem,
                     "epitope_seq": ("" if pd.isna(r.get("iedb_epitope")) else str(r["iedb_epitope"])),
                     "organisms": ([] if pd.isna(r.get("iedb_organism")) else [str(r["iedb_organism"])]),
                     "pmids": pmids}
    else:
        m["iedb"] = {"matched": False}
    p.metrics.update(m)
    peptides.append(p)

paths = report.write_package(run_dir, peptides, meta)
print(f"yeniden üretildi ({lang}):", len(peptides), "aday")
for k, v in paths.items():
    print(f"  {k}: {v}")
