"""Adım 3b — Suşlar arası konservasyon (korunmuşluk).

Bir aday antijen suşlar arasında ne kadar KORUNMUŞ? Geniş koruma için yüksek
istenir (varyabl antijen → dar-suş koruması; ör. MenB fHbp varyant-zengin,
Vaccine Design Ch25). Tek proteom tek suştur → karşılaştıracak bir şey yoktur;
bu yüzden konservasyon YALNIZ kullanıcı ek suş proteom(lar)ı verirse hesaplanır.

YÖNTEM (hafif, GPU'suz — yalnız funnel-survivor adaylarda koşar):
  1. Her ek suş proteomundan bir kez diamond DB kur (önbellek: tools/db/.strain_cache/).
  2. Her aday proteini her suşa 'diamond blastp' ile tara → en iyi ortolog (sseqid).
  3. Aday ↔ ortolog BLOSUM62 lokal hizalama → aday kalıntılarına korunmuşluk yansıt.
  4. Kalıntı-başı konservasyon = (kendi suşu + ortologların aynı olan oranı);
     protein skoru = kalıntı ortalaması. Ek olarak kolon Shannon entropisi (düşük=korunmuş).

Suş yoksa / diamond yoksa → DÜRÜST 'hesaplanmadı' (uydurma yok). Yedek deseni
allergen/discovery ile aynı; 'method_conservation' etiketlidir.
"""

from __future__ import annotations

import hashlib
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import Align, SeqIO

from .models import ProteinRecord
from .sequtils import sanitize

_TOOLS = Path(__file__).resolve().parent.parent / "tools"
_CACHE_DIR = _TOOLS / "db" / ".strain_cache"

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _diamond_bin() -> str | None:
    local = _TOOLS / "diamond"
    if local.exists():
        return str(local)
    return shutil.which("diamond")


def available(strain_paths: list[str] | None) -> bool:
    """Konservasyon bu koşuda hesaplanabilir mi? (suş dosyası + diamond)."""
    return bool(strain_paths) and _diamond_bin() is not None


def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    return a


def _make_db(strain_path: str) -> tuple[str, dict[str, str]] | None:
    """Suş proteomundan diamond DB kur (önbellekli) + {id: seq} sözlüğü döndür.

    Dönüş: (dmnd_yolu, id->seq). Hata olursa None.
    """
    p = Path(strain_path)
    if not p.exists():
        return None
    try:
        seqs = {r.id: sanitize(str(r.seq)) for r in SeqIO.parse(str(p), "fasta")}
    except Exception:
        return None
    if not seqs:
        return None
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key = hashlib.md5(f"{p.resolve()}:{p.stat().st_size}".encode()).hexdigest()[:16]
    dmnd = _CACHE_DIR / f"{key}.dmnd"
    if not dmnd.exists():
        with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
            for sid, s in seqs.items():
                fh.write(f">{sid}\n{s}\n")
            faa = fh.name
        try:
            r = subprocess.run(
                [str(_diamond_bin()), "makedb", "--in", faa, "-d", str(dmnd), "--quiet"],
                capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                Path(faa).unlink(missing_ok=True)
                return None
        except Exception:
            Path(faa).unlink(missing_ok=True)
            return None
        Path(faa).unlink(missing_ok=True)
    return str(dmnd), seqs


def _best_orthologs(proteins: list[ProteinRecord], dmnd: str) -> dict[str, str]:
    """qseqid -> en iyi ortolog sseqid (bir suş DB'sinde)."""
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for pr in proteins:
            fh.write(f">{pr.id}\n{sanitize(pr.seq)}\n")
        qpath = fh.name
    try:
        cmd = [str(_diamond_bin()), "blastp", "-q", qpath, "-d", dmnd,
               "--outfmt", "6", "qseqid", "sseqid", "bitscore",
               "--max-target-seqs", "1", "--evalue", "1e-5", "--quiet"]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception:
        Path(qpath).unlink(missing_ok=True)
        return {}
    Path(qpath).unlink(missing_ok=True)
    best: dict[str, tuple[str, float]] = {}
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 3:
            continue
        qid, sid, bits = p[0], p[1], float(p[2])
        if qid not in best or bits > best[qid][1]:
            best[qid] = (sid, bits)
    return {q: v[0] for q, v in best.items()}


def _identity_vector(cand: str, ortholog: str, aligner) -> list[int] | None:
    """Aday her kalıntısı ortologda aynı mı? (1/0). Hizalama başarısızsa None."""
    s = sanitize(cand)
    try:
        aln = aligner.align(s, ortholog)[0]
    except Exception:
        return None
    a, b = aln[0], aln[1]  # aday-hizalı, ortolog-hizalı
    vec = [0] * len(s)
    ci = 0
    for x, y in zip(a, b):
        if x != "-":
            if y != "-" and x == y:
                vec[ci] = 1
            ci += 1
    return vec


def _col_entropy(cand: str, aligned_cols: list[list[str]]) -> list[float]:
    """Her aday kolonundaki amino asit dağılımının Shannon entropisi (bit). Düşük=korunmuş."""
    out = []
    for i in range(len(cand)):
        counts: dict[str, int] = {}
        col = [cand[i]] + [c[i] for c in aligned_cols if c[i] != "-"]
        for aa in col:
            if aa in _AA:
                counts[aa] = counts.get(aa, 0) + 1
        n = sum(counts.values()) or 1
        h = -sum((c / n) * math.log2(c / n) for c in counts.values())
        out.append(round(h, 3))
    return out


def run(proteins: list[ProteinRecord], strain_paths: list[str] | None,
        min_percent: float = 80.0) -> tuple[list[ProteinRecord], dict]:
    """Funnel-survivor adaylara suş-konservasyonu yaz. Dönüş: (proteins, summary).

    Suş yoksa dürüstçe conservation_percent=None ('hesaplanmadı') bırakır.
    """
    if not available(strain_paths):
        for pr in proteins:
            pr.annotations["conservation_percent"] = None
            pr.annotations["conservation_profile"] = None
            pr.annotations["method_conservation"] = "hesaplanmadı (suş verisi yok)"
        return list(proteins), {"computed": False, "n_strains": 0,
                                "note": "suş verisi verilmedi"}

    # 1) Suş DB'leri + sözlükleri hazırla (önbellekli)
    strains: list[tuple[str, dict[str, str]]] = []
    for sp in strain_paths:
        db = _make_db(sp)
        if db:
            strains.append(db)
    if not strains:
        for pr in proteins:
            pr.annotations["conservation_percent"] = None
            pr.annotations["conservation_profile"] = None
            pr.annotations["method_conservation"] = "hesaplanmadı (suş DB kurulamadı)"
        return list(proteins), {"computed": False, "n_strains": 0,
                                "note": "suş DB kurulamadı"}

    aligner = _aligner()
    n_strains = len(strains)

    # 2) Her suşta en iyi ortologları bul (aday listesi tek batch)
    per_strain_hits = []  # [(dmnd_seqs, {qid: sseqid})]
    for dmnd, seqs in strains:
        hits = _best_orthologs(proteins, dmnd)
        per_strain_hits.append((seqs, hits))

    conserved_counts = 0
    for pr in proteins:
        cand = sanitize(pr.seq)
        L = len(cand)
        # kendi suşu her kalıntıda 'aynı' (taban 1); ortologlar eklenir
        agree = [1] * L      # her pozisyonda kaç suşta aynı (kendi dahil)
        total = [1] * L      # her pozisyonda kaç suş kıyaslandı (kendi dahil)
        aligned_cols: list[list[str]] = []
        n_orthologs = 0
        for seqs, hits in per_strain_hits:
            sid = hits.get(pr.id)
            if not sid or sid not in seqs:
                continue  # bu suşta ortolog yok → o suş bu pozisyonlara katkı vermez
            vec = _identity_vector(cand, seqs[sid], aligner)
            if vec is None:
                continue
            n_orthologs += 1
            for i in range(L):
                total[i] += 1
                agree[i] += vec[i]
            # entropi için hizalı ortolog karakterlerini kolonlara yaz
            aligned_cols.append(_aligned_ortholog_chars(cand, seqs[sid], aligner))

        profile = [round(agree[i] / total[i], 3) for i in range(L)]
        pct = round(100 * sum(profile) / L, 1) if L else None
        entropy = _col_entropy(cand, aligned_cols) if aligned_cols else None
        mean_ent = round(sum(entropy) / len(entropy), 3) if entropy else None

        pr.annotations["conservation_percent"] = pct
        pr.annotations["conservation_profile"] = profile
        pr.annotations["conservation_entropy_mean"] = mean_ent
        pr.annotations["conservation_n_orthologs"] = n_orthologs
        pr.annotations["method_conservation"] = (
            f"GERÇEK (diamond ortolog + BLOSUM62 hizalama, {n_orthologs}/{n_strains} suş)"
            if n_orthologs else f"ortolog bulunamadı ({n_strains} suşta)")
        if pct is not None and pct >= min_percent:
            conserved_counts += 1

    return list(proteins), {
        "computed": True, "n_strains": n_strains,
        "n_conserved": conserved_counts, "min_percent": min_percent,
    }


def _aligned_ortholog_chars(cand: str, ortholog: str, aligner) -> list[str]:
    """Aday her kalıntısına denk gelen ortolog karakteri ('-' = boşluk/ortolog yok)."""
    s = sanitize(cand)
    out = ["-"] * len(s)
    try:
        aln = aligner.align(s, ortholog)[0]
    except Exception:
        return out
    a, b = aln[0], aln[1]
    ci = 0
    for x, y in zip(a, b):
        if x != "-":
            if ci < len(out):
                out[ci] = y
            ci += 1
    return out
