"""İleriye dönük (post-2019-07-31) ground-truth setini kur ve sırasını hesapla.

Gerekçe: kıyas ground-truth'umuzun %59'u Vaxign-ML'in eğitim etiketi
(benchmark/results/vaxignml_circularity.json). Bu dairesellik, karşılaştırmanın
genelleme hakkında söyleyebileceğini sınırlıyor. Çözüm, Vaxign-ML'in eğitim kesiti
olan 31-07-2019'dan SONRA yayımlanmış ve Protegen'e uzak antijenlerle ikinci,
bağımsız bir set kurmak.

Üç adım, hiçbiri pipeline koşusu gerektirmiyor:
  1) Küratörlük kaydındaki antijenin dizisi kaynak proteomdan çekilir ve diamond ile
     HEDEF proteomda (bizim koştuğumuz suş) karşılığı bulunur. Ad eşlemesi yapılmaz:
     'suis' regex'inin 28 Glaesserella parasuis kaydını yakaladığı hata tekrarlanmasın.
  2) Aynı dizi Protegen v4.0'a karşı hizalanır; >=%90 kimlik = Vaxign-ML'in eğitim
     etiketi, %30-90 = homolog (onun negatif-seçim eşiğinin üstü), <%30 = bağımsız.
  3) Sıra, diskteki outputs/validation_<bug>/run.json aday listesinden yeniden
     hesaplanır: protein başına en yüksek candidacy, azalan sırada (run_bacterium.py
     ile özdeş tanım).

Koşum:
    .venv/bin/python3 benchmark/build_prospective_groundtruth.py ssuis \
        --protegen benchmark/downloads/protegen/protegen-bacterium-4.0.faa
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "benchmark"

# Protegen örtüşme sınıfları: Vaxign-ML'in kendi negatif-seçim eşiği %30 olduğu için
# %30-90 bandı "negatif olarak da seçilemeyecek kadar yakın" demek, yani bilgi sızdırır.
TRAIN_LABEL_ID = 90.0
HOMOLOGUE_ID = 30.0
# Hedef proteomda "aynı protein" sayılma eşiği: aynı türün başka suşu, ortolog değil.
SAME_PROTEIN_ID = 90.0


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


def _unwrap(acc: str) -> str:
    """UniProt başlığından çıplak aksesyon: tr|A0A0H3MU02|A0A0H3MU02_STRS4 -> A0A0H3MU02."""
    parts = acc.split("|")
    return parts[1] if len(parts) >= 3 else acc


def blast_best(queries: dict[str, str], subject: Path, diamond: str,
               extra: tuple[str, ...] = ()) -> dict[str, tuple[str, float, float]]:
    """Her sorgu için en iyi hit: {query: (subject_id, identity, bitscore)}."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        qf = td / "query.faa"
        qf.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items()))
        db = td / "subject"
        subprocess.run([diamond, "makedb", "--in", str(subject), "-d", str(db)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = td / "hits.tsv"
        subprocess.run([diamond, "blastp", "-q", str(qf), "-d", str(db), "-o", str(out),
                        "--very-sensitive", "--max-target-seqs", "5", "--quiet",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore",
                        *extra], check=True)
        best: dict[str, tuple[str, float, float]] = {}
        for ln in out.read_text().splitlines():
            q, s, pid, bs = ln.split("\t")[:4]
            pid, bs = float(pid), float(bs)
            if q not in best or bs > best[q][2]:
                best[q] = (s, pid, bs)
        return best


def protein_ranks(run_json: Path) -> tuple[dict[str, int], int]:
    """run_bacterium.py ile özdeş tanım: protein başına en yüksek candidacy, azalan."""
    run = json.loads(run_json.read_text(encoding="utf-8"))
    best: dict[str, float] = {}
    for c in run.get("candidates") or []:
        src = c.get("cds_source")
        if not src:
            continue
        sc = c.get("candidacy")
        if sc is None:
            continue
        best[src] = max(best.get(src, 0.0), float(sc))
    ordered = sorted(best, key=lambda a: -best[a])
    return {a: i + 1 for i, a in enumerate(ordered)}, len(ordered)


def classify(pid: float | None) -> str:
    if pid is None:
        return "bağımsız (<%30)"
    if pid >= TRAIN_LABEL_ID:
        return "eğitim etiketi (>=%90)"
    if pid >= HOMOLOGUE_ID:
        return "homolog (%30-90)"
    return "bağımsız (<%30)"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("bug")
    ap.add_argument("--protegen", required=True)
    ap.add_argument("--diamond", default=shutil.which("diamond") or "diamond")
    args = ap.parse_args()

    cur = json.loads((BENCH / "prospective" / f"{args.bug}.json").read_text(encoding="utf-8"))
    protegen = Path(args.protegen)
    if not protegen.exists():
        print(f"Protegen FASTA yok: {protegen}")
        return 1

    # 1) küratörlük kaydındaki dizileri kaynak proteomlardan çek
    sources = {k: read_faa(ROOT / v) for k, v in cur["source_proteomes"].items()}
    queries: dict[str, str] = {}
    missing = []
    for a in cur["antigens"]:
        pool = sources.get(a["strain"], {})
        seq = pool.get(a["protein_id"])
        if not seq:
            missing.append(a["protein_id"])
            continue
        queries[a["locus_tag"]] = seq
    if missing:
        print(f"kaynak proteomda bulunamadı: {missing}")

    # 2) hedef proteomda karşılığı + Protegen örtüşmesi
    target = ROOT / cur["target_proteome"]
    tgt_hits = blast_best(queries, target, args.diamond)
    pro_hits = blast_best(queries, protegen, args.diamond)

    # 3) sıra: diskteki aday listesinden
    run_json = ROOT / "outputs" / f"validation_{args.bug}" / "run.json"
    ranks, n_cand_proteins = protein_ranks(run_json) if run_json.exists() else ({}, 0)

    rows = []
    for a in cur["antigens"]:
        lt = a["locus_tag"]
        th = tgt_hits.get(lt)
        ph = pro_hits.get(lt)
        acc = _unwrap(th[0]) if th and th[1] >= SAME_PROTEIN_ID else None
        rows.append({
            **{k: a[k] for k in ("name", "strain", "protein_id", "locus_tag", "doi",
                                 "pmid", "year", "evidence", "protection", "grade",
                                 "prior_art", "caveat")},
            # Karışım olarak uygulanan antijenlerde koruma tek bileşene atfedilemez.
            "mixture": bool(a.get("mixture")),
            "identification": a.get("identification"),
            "target_accession": acc,
            "target_identity": round(th[1], 1) if th else None,
            "in_target_proteome": acc is not None,
            "protegen_hit": ph[0] if ph else None,
            "protegen_identity": round(ph[1], 1) if ph else None,
            "overlap_class": classify(ph[1] if ph else None),
            "protein_rank": ranks.get(acc) if acc else None,
            "produced_candidate": bool(acc and acc in ranks),
        })

    usable = [r for r in rows
              if r["in_target_proteome"] and not r["prior_art"] and not r["mixture"]
              and r["overlap_class"].startswith("bağımsız")]

    # Aday üretmeyen antijenlerin nedenini ölç, tahmin etme: lokalizasyon kapısı
    # kayıpların bağlayıcı kısıtı çıktıysa (bkz. MISS_ATTRIBUTION §3) bağımsız bir
    # sette de aynı şekilde çıkmalı. Ölçüm, atfı doğrulanabilir kılıyor.
    missed = [r for r in rows if r["in_target_proteome"] and not r["produced_candidate"]]
    if missed:
        sys.path.insert(0, str(ROOT))
        from vaxforge import psortb  # noqa: PLC0415
        if psortb.available():
            proteome = read_faa(target)
            pairs = [(r["target_accession"], proteome[r["target_accession"]])
                     for r in missed if r["target_accession"] in proteome]
            loc = psortb.predict(pairs, "positive")
            for r in rows:
                acc = r["target_accession"]
                if acc and acc in loc:
                    r["psortb_localization"] = loc[acc].get("localization")
                    r["psortb_score"] = loc[acc].get("score")
    rep = {
        "bug": args.bug,
        "question": ("Vaxign-ML'in eğitim kesitinden SONRA yayımlanmış ve Protegen'e "
                     "uzak koruyucu antijenlerde VaxForge nasıl sıralıyor?"),
        "cutoff": cur["cutoff"],
        "cutoff_rationale": cur["cutoff_rationale"],
        "inclusion_criteria": cur["inclusion_criteria"],
        "n_curated": len(rows),
        "n_in_target_proteome": sum(1 for r in rows if r["in_target_proteome"]),
        "n_prospective_usable": len(usable),
        "n_mixture_excluded": sum(1 for r in rows if r["mixture"]),
        "n_candidate_proteins_in_run": n_cand_proteins,
        "run_json": str(run_json.relative_to(ROOT)) if run_json.exists() else None,
        "thresholds": {"same_protein_identity": SAME_PROTEIN_ID,
                       "train_label_identity": TRAIN_LABEL_ID,
                       "homologue_identity": HOMOLOGUE_ID},
        "antigens": rows,
        "negative_evidence": cur.get("negative_evidence", []),
        "pending_full_text": cur.get("pending_full_text", []),
        "resolved_from_secondary_sources": cur.get("resolved_from_secondary_sources", []),
    }
    out = BENCH / "results" / f"prospective_groundtruth_{args.bug}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"{'antijen':44} {'hedef':11} {'kimlik':6} {'Protegen':9} {'sınıf':22} sıra")
    for r in rows:
        tail = (str(r["protein_rank"]) if r["protein_rank"]
                else f"aday üretmedi [{r.get('psortb_localization')} {r.get('psortb_score')}]")
        print(f"{r['name'][:44]:44} {str(r['target_accession']):11} "
              f"{str(r['target_identity']):6} {str(r['protegen_identity']):9} "
              f"{r['overlap_class'][:22]:22} {tail}")
    print(f"\nküratörlük: {len(rows)} | hedef proteomda: {rep['n_in_target_proteome']} | "
          f"karışım (sayıma girmez): {rep['n_mixture_excluded']} | "
          f"prospektif kullanılabilir (yeni + bağımsız): {len(usable)}")
    print(f"koşudaki aday protein sayısı: {n_cand_proteins}")
    print(f"yazıldı: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
