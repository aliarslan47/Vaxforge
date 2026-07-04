"""Tanınan dosyayı aday antijen proteinlerine (proteom) dönüştürür.

  proteome  -> protein dizileri doğrudan
  cds       -> çeviri (frame 0)
  genome    -> 6-çerçeve ORF bulma (fallback gen tahmini)
  reads     -> (prototip) 6-çerçeve ORF bulma; gerçekte assembly+Prodigal gerekir

Gerçek araçlar (SPAdes, Prodigal) takıldığında bu fallback devre dışı kalır.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .detect import Detection
from .models import ProteinRecord

MIN_ORF_AA = 40

# NCBI CDS/protein FASTA başlığındaki alanlar: [gene=..] [locus_tag=..] [location=..]
_LOCUS_RE = {
    "gene": re.compile(r"\[gene=([^\]]+)\]"),
    "locus_tag": re.compile(r"\[locus_tag=([^\]]+)\]"),
    "protein": re.compile(r"\[protein=([^\]]+)\]"),
    "protein_id": re.compile(r"\[protein_id=([^\]]+)\]"),
    "location": re.compile(r"\[location=([^\]]+)\]"),
}


def _parse_locus(desc: str) -> dict:
    """CDS başlığından gen/lokus/konum bilgisini çıkarır (varsa)."""
    out = {}
    for key, rx in _LOCUS_RE.items():
        m = rx.search(desc or "")
        if m:
            out[key] = m.group(1).strip()
    return out


def _open(path: Path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")


def _fmt(det: Detection) -> str:
    return "fastq" if det.fmt == "fastq" else "fasta"


def load_proteins(path: str | Path, det: Detection) -> list[ProteinRecord]:
    path = Path(path)
    if det.molecule == "genbank":
        return _read_genbank(path)
    if det.molecule == "proteome":
        return _read_proteins(path, det)
    if det.molecule == "cds":
        return _translate_cds(path, det)
    # genome / reads -> ORF bulma
    return _find_orfs(path, det)


def _read_genbank(path) -> list[ProteinRecord]:
    """GenBank CDS'lerini protein olarak alır; gen/lokus/konum doğrudan gelir."""
    out = []
    with _open(path) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            idx = 0
            for f in rec.features:
                if f.type != "CDS":
                    continue
                q = f.qualifiers
                prot = (q.get("translation", [None])[0])
                if not prot:
                    try:
                        prot = str(f.extract(rec.seq).translate(to_stop=True))
                    except Exception:
                        continue
                if not prot or len(prot) < 20:
                    continue
                idx += 1
                locus = q.get("locus_tag", [None])[0]
                pid = q.get("protein_id", [None])[0] or locus or f"{rec.id}_cds{idx}"
                pr = ProteinRecord(id=pid, seq=prot, source="genbank")
                pr.annotations.update({
                    "gene": q.get("gene", [None])[0],
                    "locus_tag": locus,
                    "protein_id": q.get("protein_id", [None])[0],
                    "location": str(f.location),
                    "desc": q.get("product", [""])[0],
                })
                out.append(pr)
    return out


def _read_proteins(path, det) -> list[ProteinRecord]:
    out = []
    with _open(path) as fh:
        for rec in SeqIO.parse(fh, _fmt(det)):
            seq = str(rec.seq).replace("*", "")
            if len(seq) >= 20:
                pr = ProteinRecord(id=rec.id, seq=seq, source="proteome")
                pr.annotations["desc"] = rec.description
                pr.annotations.update(_parse_locus(rec.description))
                out.append(pr)
    return out


def _translate_cds(path, det) -> list[ProteinRecord]:
    out = []
    with _open(path) as fh:
        for rec in SeqIO.parse(fh, _fmt(det)):
            nt = str(rec.seq)
            nt = nt[: len(nt) - (len(nt) % 3)]
            if not nt:
                continue
            prot = str(Seq(nt).translate(to_stop=True))
            if len(prot) >= 20:
                pr = ProteinRecord(id=rec.id, seq=prot, source="translate")
                pr.annotations["desc"] = rec.description
                pr.annotations.update(_parse_locus(rec.description))
                out.append(pr)
    return out


def _find_orfs(path, det) -> list[ProteinRecord]:
    out = []
    with _open(path) as fh:
        for rec in SeqIO.parse(fh, _fmt(det)):
            nt = str(rec.seq).upper()
            idx = 0
            for strand, s in ((+1, nt), (-1, str(Seq(nt).reverse_complement()))):
                for frame in range(3):
                    sub = s[frame:]
                    sub = sub[: len(sub) - (len(sub) % 3)]
                    prot = str(Seq(sub).translate())
                    for orf in prot.split("*"):
                        if len(orf) >= MIN_ORF_AA:
                            m = orf.find("M")
                            orf2 = orf[m:] if m >= 0 else orf
                            if len(orf2) >= MIN_ORF_AA:
                                idx += 1
                                out.append(ProteinRecord(
                                    id=f"{rec.id}_orf{idx}({'+' if strand>0 else '-'}{frame})",
                                    seq=orf2, source="orf"))
    return out
