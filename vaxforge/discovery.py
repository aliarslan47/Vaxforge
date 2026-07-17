"""Adım 2 — Virülans faktörü / aşı hedefi keşfi.

GERÇEK yol: proteomu 'diamond blastp' ile küratörlü VFDB'ye tarar
(tools/diamond + tools/db/vfdb.dmnd varsa). Eşikler config'ten (evalue/identity/
coverage). Ek olarak protein açıklamasında virülans anahtar-kelime taraması.

diamond/VFDB yoksa FALLBACK: paket-içi mini referansa yerel hizalama + anahtar-
kelime (heuristik). Sonuçlar 'method' ile etiketlenir. Sert filtre.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import Align, SeqIO

from .config_loader import ResolvedTool
from .models import ProteinRecord
from .sequtils import sanitize

_TOOLS = Path(__file__).resolve().parent.parent / "tools"
_VFDB = _TOOLS / "db" / "vfdb.dmnd"


def _diamond_bin() -> str | None:
    """DIAMOND binary yolu: önce paket-içi tools/diamond, yoksa sistem PATH'i.
    (Eskiden yalnız tools/diamond'a bakılıyordu; sistemde kurulu diamond
    görülmediği için gerçek araç boşuna 'yedek' sanılıyordu.)"""
    local = _TOOLS / "diamond"
    if local.exists():
        return str(local)
    return shutil.which("diamond")
_REF_PATH = Path(__file__).resolve().parent / "refdata" / "virulence_ref.faa"

_KEYWORDS = [
    "virulence", "adhesin", "invasin", "toxin", "hemolysin", "haemolysin",
    "outer membrane", "secreted", "pilus", "fimbria", "flagell", "protease",
    "iron", "siderophore", "capsul", "effector", "porin", "lipoprotein",
]

# Viral yüzey/yapısal protein anahtar-kelimeleri (yalnız açıklama etiketlemek için;
# virüs profilinde keşif zaten baypas edilir).
_VIRAL_KEYWORDS = [
    "envelope", "glycoprotein", "spike", "capsid", "nucleoprotein", "nucleocapsid",
    "surface", "fusion", "hemagglutinin", "haemagglutinin", "neuraminidase",
    "matrix", "membrane protein", "polyprotein", "receptor-binding",
]


def diamond_available() -> bool:
    return _diamond_bin() is not None and _VFDB.exists()


# runner/tool_status uyumu için standart isim (eskiden yalnız diamond_available
# vardı; runner available() çağırıp bulamayınca araç yanlışlıkla 'yedek' oluyordu).
def available() -> bool:
    return diamond_available()


def _virulence_score(vf_hit: bool, pident: float, kw: str | None) -> float:
    """Virülans kanıtını 0-1 skora indirger (sert kapı DEĞİL, adaylık bileşeni).

    Güçlü VFDB eşleşmesi → yüksek; zayıf/kısmi benzerlik → orta; yalnız küratörlü
    anahtar-kelime → 0.5; hiç kanıt yok → 0.15 (kanıt yokluğu, virülans-değil kanıtı
    değildir — bu yüzden taban sıfır değil). NERVE2'nin virülans-skoru mantığı.
    """
    score = 0.0
    if vf_hit:
        score = 0.7 + 0.3 * min(1.0, pident / 100.0)
    elif pident > 0:
        score = 0.2 + 0.5 * min(1.0, pident / 100.0)
    if kw:
        score = max(score, 0.5)
    return round(score if score > 0 else 0.15, 3)


def _run_diamond(proteins: list[ProteinRecord], evalue: float, min_id: float,
                 min_cov: float) -> dict[str, dict]:
    """qseqid -> {pident, coverage, hit, stitle}. Hata olursa boş döner."""
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for pr in proteins:
            fh.write(f">{pr.id}\n{sanitize(pr.seq)}\n")
        qpath = fh.name
    lengths = {pr.id: len(sanitize(pr.seq)) for pr in proteins}
    try:
        cmd = [str(_diamond_bin()), "blastp", "-q", qpath, "-d", str(_VFDB),
               "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "stitle",
               "--max-target-seqs", "1", "--evalue", str(evalue), "--quiet"]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception:
        Path(qpath).unlink(missing_ok=True)
        return {}
    Path(qpath).unlink(missing_ok=True)
    best: dict[str, dict] = {}
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 6:
            continue
        qid, sid, pident, length, ev, stitle = p[0], p[1], float(p[2]), int(p[3]), float(p[4]), p[5]
        cov = 100 * length / max(1, lengths.get(qid, length))
        if qid in best and best[qid]["pident"] >= pident:
            continue
        best[qid] = {"pident": round(pident, 1), "coverage": round(cov, 1),
                     "hit": pident >= min_id and cov >= min_cov, "stitle": stitle[:80], "sseqid": sid}
    return best


# --- fallback: mini-ref hizalama --------------------------------------------
def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    return a


def _mini_ref_identity(seq: str, refs, aligner) -> tuple[float, str]:
    s = sanitize(seq)
    best_pct, best_ref = 0.0, ""
    for rid, rseq in refs:
        try:
            aln = aligner.align(s, rseq)[0]
        except Exception:
            continue
        a, b = aln[0], aln[1]
        matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
        pct = 100 * matches / max(1, min(len(s), len(rseq)))
        if pct > best_pct:
            best_pct, best_ref = pct, rid
    return round(best_pct, 1), best_ref


def run(proteins: list[ProteinRecord], tool: ResolvedTool,
        profile: str = "bacteria") -> list[ProteinRecord]:
    # VFDB bakteriyel virülans faktörleri DB'sidir; viral/parazit proteinleri
    # buraya vurmaz. Bu profillerde keşfi SERT filtre olarak uygulamak her şeyi
    # eler. Bu yüzden bakteri-dışında keşif baypas edilir: tüm proteinler aday
    # geçer (funnel + epitop adımları eler), viral anahtar-kelimeler etiketlenir.
    if profile != "bacteria":
        for pr in proteins:
            desc = str(pr.annotations.get("desc", "")).lower()
            kw = next((k for k in _VIRAL_KEYWORDS + _KEYWORDS if k in desc), None)
            # Virüs/parazitte VFDB anlamsız; virülans skoru anahtar-kelimeden
            # (yüzey/yapısal protein etiketi) yaklaşık — 0.15 taban, kw varsa 0.5.
            pr.annotations.update({
                "vf_identity": 0.0, "vf_ref": "", "vf_keyword": kw or "",
                "virulence": 0.5 if kw else 0.15,
                "method_discovery": f"ATLANDI ({profile}: VFDB bakteriye özgü) — tüm proteinler aday",
            })
            pr.mark("discovery", True,
                    (f"keşif atlandı ({profile}); anahtar-kelime: {kw}" if kw
                     else f"keşif atlandı ({profile}); tüm proteinler aday"),
                    vf_identity=0.0, keyword=kw)
        return list(proteins)

    min_id = float(tool.params["min_identity"].value)
    min_cov = float(tool.params["min_coverage"].value)
    evalue = float(tool.params["evalue"].value)

    use_diamond = diamond_available()
    dia = _run_diamond(proteins, evalue, min_id, min_cov) if use_diamond else {}
    if use_diamond and not dia:
        use_diamond = False  # diamond çalışmadı -> fallback
    refs = None if use_diamond else [(r.id, sanitize(str(r.seq))) for r in SeqIO.parse(_REF_PATH, "fasta")]
    aligner = None if use_diamond else _aligner()

    # SKORLAMA adımı (NERVE2 gibi) — ELEME DEĞİL: hiçbir protein burada düşmez.
    # Her proteine virülans skoru atanır; gerçek eleme funnel'da (lokalizasyon/TM/
    # insan-homoloji sert filtreleri) yapılır. Böylece VFDB'de olmayan gerçek yüzey
    # antijenleri (ör. N. meningitidis NHBA) kaybolmaz.
    for pr in proteins:
        desc = str(pr.annotations.get("desc", "")).lower()
        kw = next((k for k in _KEYWORDS if k in desc), None)

        if use_diamond:
            d = dia.get(pr.id, {"pident": 0.0, "coverage": 0.0, "hit": False, "stitle": ""})
            vf_hit, pct, ref = d["hit"], d["pident"], d.get("stitle", "")
            method = "GERÇEK (diamond blastp → VFDB) + anahtar-kelime"
            reason = (f"VFDB %{pct} kimlik ({ref})" if vf_hit
                      else f"anahtar-kelime: {kw}" if kw else f"VFDB hit yok (%{pct})")
        else:
            pct, ref = _mini_ref_identity(pr.seq, refs, aligner)
            vf_hit = pct >= min_id
            method = "heuristik (mini-ref + anahtar-kelime) — VFDB yok"
            reason = (f"mini-ref %{pct} ({ref})" if vf_hit
                      else f"anahtar-kelime: {kw}" if kw else f"mini-ref %{pct}, anahtar-kelime yok")

        vir = _virulence_score(vf_hit, pct, kw)
        pr.annotations.update({"vf_identity": pct, "vf_ref": ref, "vf_keyword": kw or "",
                               "virulence": vir, "vf_hit": vf_hit,
                               "method_discovery": method})
        # Soft: her zaman geçer (sert kapı değil); kanıt 'reason'da + virülans skorunda.
        pr.mark("discovery", True, f"{reason} → virülans skoru {vir}",
                vf_identity=pct, keyword=kw, virulence=vir)
    return list(proteins)
