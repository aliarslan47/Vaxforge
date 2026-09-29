"""7 aday virüsün RefSeq proteomunu indir + boyuta göre sırala (koşmadan önce gösterim).

Kural: aksesyonlar NCBI'dan doğrulanır (uydurma yok). Poliprotein virüslerde
mat_peptide çıkarılır (EMCV kalıbı); segmentli/çok-CDS virüslerde tüm CDS translation'ları.
Çıktı: data/validation/<slug>_proteome.faa + özet tablo (protein sayısı, residü, en küçük→büyük).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

Entrez.email = "aliarslan47@gmail.com"
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "validation"

# (slug, gösterim adı, taxid, aile, poliprotein mi)
CANDIDATES = [
    ("reovirus",  "Mammalian orthoreovirus (Reovirus T3)",      351073, "Reoviridae (segmentli)",  False),
    ("cowpox",    "Cowpox virus",                                10243,  "Poxviridae",              False),
    ("vaccinia",  "Vaccinia virus",                              10245,  "Poxviridae",              False),
    ("rvfv",      "Rift Valley fever virus",                     11588,  "Phenuiviridae (segmentli)", False),
    ("veev",      "Venezuelan equine encephalitis virus",        11036,  "Togaviridae (poliprotein)", True),
    ("jev",       "Japanese encephalitis virus",                 11072,  "Flaviviridae (poliprotein)", True),
    ("getah",     "Getah virus",                                 59300,  "Togaviridae (poliprotein)", True),
]


def find_refseq_accessions(taxid: int) -> list[str]:
    """Taxon için RefSeq NC_ nuccore aksesyonlarını bul."""
    term = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP]"
    h = Entrez.esearch(db="nuccore", term=term, retmax=50)
    rec = Entrez.read(h); h.close()
    ids = rec["IdList"]
    if not ids:
        return []
    h = Entrez.esummary(db="nuccore", id=",".join(ids))
    summ = Entrez.read(h); h.close()
    accs = []
    for d in summ:
        acc = d.get("AccessionVersion", "")
        title = d.get("Title", "")
        if acc.startswith("NC_"):
            accs.append((acc, title, int(d.get("Length", 0))))
    return accs


def extract_proteins(accs: list, polyprotein: bool):
    """RefSeq genom(lar)dan protein sekanslarını çıkar."""
    prots = []  # (id, name, seq)
    for acc, title, _ln in accs:
        h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
        gb = SeqIO.read(h, "genbank"); h.close()
        mats = [f for f in gb.features if f.type == "mat_peptide"]
        cds = [f for f in gb.features if f.type == "CDS"]
        feats = mats if (polyprotein and mats) else cds
        for f in feats:
            q = f.qualifiers
            name = (q.get("product", q.get("gene", ["?"]))[0])
            if "translation" in q:
                seq = q["translation"][0]
            else:
                try:
                    seq = str(f.extract(gb.seq).translate(to_stop=True))
                except Exception:
                    continue
            pid = q.get("protein_id", [f"{acc}:{f.location.start}"])[0]
            if len(seq) >= 8:
                prots.append((pid, name, seq))
        time.sleep(0.34)
    return prots


def main():
    rows = []
    for slug, name, taxid, family, poly in CANDIDATES:
        try:
            accs = find_refseq_accessions(taxid)
            if not accs:
                print(f"[!] {name}: RefSeq bulunamadı (taxid {taxid})", file=sys.stderr)
                continue
            prots = extract_proteins(accs, poly)
            residues = sum(len(s) for _, _, s in prots)
            outp = OUT / f"{slug}_proteome.faa"
            with open(outp, "w") as fh:
                for pid, pname, seq in prots:
                    fh.write(f">{pid} {pname}\n{seq}\n")
            acc_str = ",".join(a for a, _, _ in accs)
            rows.append((slug, name, family, len(prots), residues, acc_str))
            print(f"[ok] {name}: {len(prots)} protein, {residues} residü  ({acc_str})", flush=True)
        except Exception as e:
            print(f"[ERR] {name}: {e}", file=sys.stderr, flush=True)

    rows.sort(key=lambda r: r[4])  # residüye göre en küçük -> büyük
    print("\n=== SIRALI (en küçük proteom -> en büyük) ===")
    print(f"{'#':<3}{'slug':<10}{'virüs':<42}{'protein':>8}{'residü':>9}")
    for i, (slug, name, fam, np_, res, acc) in enumerate(rows, 1):
        print(f"{i:<3}{slug:<10}{name[:40]:<42}{np_:>8}{res:>9}")


if __name__ == "__main__":
    main()
