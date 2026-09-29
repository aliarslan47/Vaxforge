"""Sorunlu 4 virüs proteomunu TEMİZDEN hazırla (tek-referans + doğru ground-truth).

Problem (2026-09-29 denetimi):
- VEEV & getah: diskteki faa yalnız POLİPROTEİN içeriyor → E2 ayrı dizi değil. mat_peptide çıkarımı gerek.
- vaccinia & cowpox: çok-referans birleşik (408/233) → tek RefSeq'e indirilecek.
- Ground-truth L1R (pox): vaccinia'da ürün etiketiyle; cowpox'ta jenerik CPXV adları → vaccinia L1R'ye
  ORTOLOG eşleştirmeyle (en yüksek kimlik) çöz; güvenli değilse 'unresolved' (çok-konak kolu yine koşar).

Çıktı: temiz data/validation/{veev,getah,vaccinia,cowpox}_proteome.faa (üzerine yaz)
       + benchmark/new_virus_groundtruth.json (slug -> gt_id/gt_name/gt_acc/taxon)
Kural: aksesyon NCBI'dan (uydurma yok). Bkz [[vaxforge-competitor-benchmark]] [[vaxforge-literature-rule]].
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

Entrez.email = "aliarslan47@gmail.com"
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "validation"
GTJSON = ROOT / "benchmark" / "new_virus_groundtruth.json"

# Pinlenen tek RefSeq aksesyonları (deterministik — çok-referans karışmasın)
PINNED = {
    "veev":     ("NC_001449.1", 11036, True,  "E2"),          # Togaviridae poliprotein → E2
    "getah":    ("NC_006558.1", 59300, True,  "E2"),          # Togaviridae poliprotein → E2
    "vaccinia": ("NC_006998.1", 10245, False, "L1R"),         # Poxviridae WR → L1R
    "cowpox":   ("NC_003663.2", 10243, False, "L1R"),         # Poxviridae Brighton Red → L1R (ortolog)
}


def fetch_gb(acc: str):
    for attempt in range(4):
        try:
            h = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
            gb = SeqIO.read(h, "genbank"); h.close()
            return gb
        except Exception as e:
            print(f"  [retry {attempt+1}] {acc}: {e}", flush=True)
            time.sleep(3)
    raise RuntimeError(f"efetch başarısız: {acc}")


def extract(gb, polyprotein: bool):
    """(pid, name, seq) listesi. Poliprotein → mat_peptide; değilse CDS translation."""
    mats = [f for f in gb.features if f.type == "mat_peptide"]
    cds = [f for f in gb.features if f.type == "CDS"]
    feats = mats if (polyprotein and mats) else cds
    prots, seen = [], set()
    for i, f in enumerate(feats):
        q = f.qualifiers
        name = q.get("product", q.get("gene", ["?"]))[0]
        if "translation" in q:
            seq = q["translation"][0]
        else:
            try:
                seq = str(f.extract(gb.seq).translate(to_stop=True))
            except Exception:
                continue
        pid = q.get("protein_id", [f"{gb.id}:mp{i}"])[0]
        if len(seq) >= 8 and pid not in seen:
            seen.add(pid)
            prots.append((pid, name, seq))
    return prots


def resolve_by_product(prots, patterns):
    """Ürün adında pattern eşleşen ilk proteini döndür (pid, name)."""
    for pid, name, _ in prots:
        for pat in patterns:
            if re.search(pat, name, re.I):
                return pid, name
    return None, None


def resolve_ortholog(prots, query_seq, min_identity=0.40):
    """query_seq'e en yüksek kimlikli proteini bul (cowpox L1R için)."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = None
    best = (None, None, 0.0)
    for pid, name, seq in prots:
        try:
            sc = aligner.score(query_seq, seq)
            ident = sc / max(len(query_seq), 1)  # kaba benzerlik oranı
        except Exception:
            continue
        if ident > best[2]:
            best = (pid, name, ident)
    if best[0] and best[2] >= min_identity:
        return best[0], best[1], round(best[2], 3)
    return None, None, round(best[2], 3)


def main():
    gt_map = {}
    vaccinia_l1r_seq = None

    # vaccinia'yı ÖNCE işle (cowpox ortolog sorgusu için L1R dizisi lazım)
    order = ["veev", "getah", "vaccinia", "cowpox"]
    for slug in order:
        acc, taxid, poly, gt_kind = PINNED[slug]
        print(f"\n=== {slug} ({acc}) ===", flush=True)
        try:
            gb = fetch_gb(acc)
            prots = extract(gb, poly)
            residues = sum(len(s) for _, _, s in prots)
            outp = OUT / f"{slug}_proteome.faa"
            with open(outp, "w") as fh:
                for pid, name, seq in prots:
                    fh.write(f">{pid} {name}\n{seq}\n")
            print(f"  yazıldı: {len(prots)} protein, {residues} residü -> {outp.name}", flush=True)

            # ground-truth çöz
            gt_id = gt_name = None
            gt_note = ""
            if gt_kind == "E2":
                gt_id, gt_name = resolve_by_product(prots, [r"\bE2\b", r"E2 (envelope|glyco)"])
                gt_note = "mat_peptide E2 (envelope glikoprotein)"
            elif gt_kind == "L1R":
                if slug == "vaccinia":
                    gt_id, gt_name = resolve_by_product(prots, [r"\bL1R\b", r"L1 "])
                    if gt_id:
                        vaccinia_l1r_seq = next(s for p, n, s in prots if p == gt_id)
                    gt_note = "ürün-etiketi L1R (IMV zar proteini)"
                else:  # cowpox — ortolog
                    if vaccinia_l1r_seq:
                        gt_id, gt_name, ident = resolve_ortholog(prots, vaccinia_l1r_seq)
                        gt_note = f"vaccinia-L1R ortoloğu (kaba kimlik {ident})"
                        if gt_id is None:
                            gt_note = f"L1R ortoloğu güvenli değil (max {ident}) → unresolved; çok-konak kolu yine koşar"
                    else:
                        gt_note = "vaccinia L1R çözülemedi → cowpox ortoloğu atlandı"
            gt_map[slug] = {
                "proteome": str(outp), "taxon": f"NCBITaxon:{taxid}",
                "n_protein": len(prots), "accession": acc,
                "gt_id": gt_id, "gt_name": gt_name or gt_kind, "gt_acc": gt_id or "",
                "note": gt_note,
            }
            print(f"  ground-truth: {gt_id} · {gt_name} · {gt_note}", flush=True)
        except Exception as e:
            print(f"  [ERR] {slug}: {e}", file=sys.stderr, flush=True)
            gt_map[slug] = {"error": str(e)}
        time.sleep(0.5)

    GTJSON.write_text(json.dumps(gt_map, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nground-truth haritası -> {GTJSON}")


if __name__ == "__main__":
    main()
