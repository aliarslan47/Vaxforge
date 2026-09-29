"""prep_new_viruses.py sonrası 3 kalan sorunu çöz (getah + pox L1R ortolog).

- getah NC_006558: mat_peptide yok (yalnız 1 transframe) → 3 CDS'i yaz, ground-truth =
  yapısal poliprotein 'C-P62-6K-E1' (E2 içerir; RefSeq'te ayrı anote DEĞİL — dürüst not).
- vaccinia/cowpox: L1R ürün-adıyla etiketsiz → kanonik vaccinia L1R proteini (RefSeq protein)
  çekilir, hem vaccinia hem cowpox proteomuna ORTOLOG (en yüksek kimlik) eşleştirilir.
new_virus_groundtruth.json güncellenir. Bkz [[vaxforge-literature-rule]].
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

Entrez.email = "aliarslan47@gmail.com"
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "validation"
GTJSON = ROOT / "benchmark" / "new_virus_groundtruth.json"

gt_map = json.loads(GTJSON.read_text())


def read_faa(slug):
    prots = []
    pid = name = None
    seq = []
    for line in (OUT / f"{slug}_proteome.faa").read_text().splitlines():
        if line.startswith(">"):
            if pid:
                prots.append((pid, name, "".join(seq)))
            h = line[1:].split(" ", 1)
            pid, name, seq = h[0], (h[1] if len(h) > 1 else ""), []
        else:
            seq.append(line.strip())
    if pid:
        prots.append((pid, name, "".join(seq)))
    return prots


# ---- 1) GETAH: CDS ile yeniden yaz, gt = yapısal poliprotein ----
print("=== getah (NC_006558.1) CDS ===", flush=True)
gb = SeqIO.read(Entrez.efetch(db="nuccore", id="NC_006558.1", rettype="gb", retmode="text"), "genbank")
prots = []
gt_id = None
for f in [x for x in gb.features if x.type == "CDS"]:
    q = f.qualifiers
    name = q.get("product", ["?"])[0]
    seq = q["translation"][0] if "translation" in q else str(f.extract(gb.seq).translate(to_stop=True))
    pid = q.get("protein_id", ["?"])[0]
    if len(seq) >= 8:
        prots.append((pid, name, seq))
        if "E1" in name or "P62" in name or "structural" in name.lower():
            gt_id = pid
            gt_name = name
with open(OUT / "getah_proteome.faa", "w") as fh:
    for pid, name, seq in prots:
        fh.write(f">{pid} {name}\n{seq}\n")
gt_map["getah"].update({
    "n_protein": len(prots), "gt_id": gt_id, "gt_acc": gt_id or "",
    "gt_name": "structural polyprotein (C-P62-6K-E1, E2 içerir)",
    "note": "RefSeq'te mat_peptide yok → E2 ayrı değil; yapısal poliprotein ground-truth (dürüst). "
            "Zayıf-kanıt virüs (subunit patenti yok).",
})
print(f"  {len(prots)} CDS yazıldı, gt={gt_id}", flush=True)

# ---- 2) Kanonik vaccinia L1R proteini çek ----
print("\n=== kanonik vaccinia L1R ===", flush=True)
h = Entrez.esearch(db="protein", term="L1R[Gene] AND Vaccinia virus[Organism] AND refseq[filter]", retmax=5)
ids = Entrez.read(h)["IdList"]; h.close()
l1r_seq = None
if ids:
    rec = SeqIO.read(Entrez.efetch(db="protein", id=ids[0], rettype="fasta", retmode="text"), "fasta")
    l1r_seq = str(rec.seq)
    print(f"  L1R referans: {rec.id} · {len(l1r_seq)} aa · {rec.description[:60]}", flush=True)
else:
    print("  [!] L1R referansı bulunamadı", flush=True)

# ---- 3) Ortolog eşleştir (vaccinia + cowpox) ----
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1


def best_ortholog(prots, query):
    self_sc = aligner.score(query, query)
    best = (None, None, 0.0)
    for pid, name, seq in prots:
        try:
            r = aligner.score(query, seq) / self_sc  # kendine-normalize benzerlik
        except Exception:
            continue
        if r > best[2]:
            best = (pid, name, round(r, 3))
    return best


if l1r_seq:
    for slug in ("vaccinia", "cowpox"):
        prots = read_faa(slug)
        pid, name, r = best_ortholog(prots, l1r_seq)
        ok = r >= 0.40
        gt_map[slug].update({
            "gt_id": pid if ok else None, "gt_acc": pid if ok else "",
            "gt_name": f"L1R ortoloğu ({name})" if ok else "L1R (çözülemedi)",
            "note": f"kanonik vaccinia-L1R'ye ortolog benzerlik {r}" + ("" if ok else " <0.40 → unresolved; çok-konak kolu yine koşar"),
        })
        print(f"  {slug}: {pid} · {name} · benzerlik {r} · {'OK' if ok else 'DÜŞÜK'}", flush=True)

GTJSON.write_text(json.dumps(gt_map, indent=2, ensure_ascii=False), encoding="utf-8")
print("\ngüncellendi:", GTJSON)
