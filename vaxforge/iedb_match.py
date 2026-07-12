"""Adım 6d — Aday peptidlerin IEDB'de (bilinen/literatürde) taranması.

AMAÇ: Pipeline'ın en sonunda kalan aday peptidleri, deneysel olarak DOĞRULANMIŞ
epitoplarla (IEDB — Immune Epitope Database) dizi düzeyinde eşleştirir. Bir aday
bilinen bir epitopla örtüşüyorsa, bu güçlü bir POZİTİF KONTROL sinyalidir ve
raporda kaynak makale (PubMed) ile birlikte gösterilir.

DÜRÜSTLÜK / TASARIM KARARLARI:
  * LLM literatür madenciliği YOK. Eşleştirme tamamen deterministik dizi
    karşılaştırmasıdır (exact + örtüşme); kaynak = IEDB'nin küratörlü kayıtları.
  * Bu adım ADAYLIK PUANINI DEĞİŞTİRMEZ (salt yorumlama, anchor-motifi gibi).
    Sebep: aynı IEDB verisini validasyonda 'ground truth' olarak kullanıyoruz;
    skora katarsak recall ölçümü döngüsel/taraflı olurdu.
  * Veri IQ-API'den (query-api.iedb.org) çekilir ve YEREL önbelleğe alınır
    (tools/db/iedb/<taxon>.json) — VFDB'yi bir kez indirmek gibi. İlk çekimden
    sonra çevrimdışı çalışır. Ağ yoksa ve önbellek yoksa dürüstçe 'çalıştırılmadı'.

İki mod:
  (1) Hedefli (organizma taxon'u verilirse): o kaynak-organizmanın tüm lineer
      pozitif epitopları çekilir, adaylar YERELDE eşleştirilir. Bu set aynı zamanda
      recall benchmark'ının paydasıdır (benchmark()).
  (2) Yedek (taxon yok): her aday, IQ-API'de canlı substring (ilike) sorgusuyla
      TÜM IEDB'ye karşı aranır.

Referans: Vita R, ve ark. The Immune Epitope Database (IEDB): 2024 update.
Nucleic Acids Res. 2025;53(D1):D436-D443.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

_API = "https://query-api.iedb.org"
_CACHE = Path(__file__).resolve().parent.parent / "tools" / "db" / "iedb"
_TIMEOUT = 60
# taxon yokken aday-başına canlı IEDB sorgusu üst sınırı (fazlası ağ-bağlı,
# pratik değil → atlanır; toplu yol için organizma taxon'u verilmeli).
LIVE_MAX = 150
_PAGE = 10000  # IQ-API sayfa başına maksimum satır

METHOD = "GERÇEK (IEDB IQ-API, deneysel epitop kaydı)"
DEFAULT_K = 8  # minimum örtüşme çekirdek uzunluğu (aa)


# --------------------------------------------------------------------------- #
# Düşük seviye API
# --------------------------------------------------------------------------- #
def _get(path: str, params: dict, headers: dict | None = None):
    """IQ-API GET → (json, content_range_str|None). Hata/timeout → (None, None)."""
    qs = urllib.parse.urlencode(params, safe="{}():,.*")
    url = f"{_API}/{path}?{qs}"
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
            body = json.loads(r.read().decode("utf-8"))
            cr = r.headers.get("Content-Range")
            return body, cr
    except Exception:
        return None, None


def online() -> bool:
    body, _ = _get("epitope_search", {"limit": 1})
    return body is not None


# --------------------------------------------------------------------------- #
# PMID çözümü
# --------------------------------------------------------------------------- #
def _resolve_pmids(ref_ids: list[int]) -> dict[int, str]:
    """IEDB reference_id -> PubMed ID (varsa). Toplu, parçalı sorgu."""
    out: dict[int, str] = {}
    uniq = sorted({int(r) for r in ref_ids if r is not None})
    for i in range(0, len(uniq), 100):
        chunk = uniq[i:i + 100]
        ids = ",".join(str(x) for x in chunk)
        body, _ = _get("reference_search",
                       {"reference_id": f"in.({ids})",
                        "select": "reference_id,pubmed_id"})
        if not body:
            continue
        for row in body:
            pm = row.get("pubmed_id")
            if pm:
                out[row["reference_id"]] = str(pm)
    return out


# --------------------------------------------------------------------------- #
# Organizma epitop setini çek + önbellekle
# --------------------------------------------------------------------------- #
def fetch_source_organism(taxon_iri: str, resolve_pmids: bool = True,
                          refresh: bool = False) -> list[dict] | None:
    """Bir kaynak-organizmanın tüm lineer POZİTİF epitoplarını çek + önbellekle.

    taxon_iri örn: 'NCBITaxon:2697049' (SARS-CoV-2). Dönen kayıt:
      {seq, length, organisms[], antigens[], reference_ids[], pmids[],
       assay_evidence[], iedb_id}
    Ağ yok ve önbellek yoksa None.
    """
    _CACHE.mkdir(parents=True, exist_ok=True)
    slug = taxon_iri.replace(":", "_")
    cache = _CACHE / f"{slug}.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text(encoding="utf-8"))

    select = ("structure_id,linear_sequence,source_organism_names,"
              "parent_source_antigen_names,reference_ids,mhc_allele_evidences,"
              "qualitative_measures")
    base = {
        "structure_type": "eq.Linear peptide",
        "qualitative_measure": "neq.Negative",
        "source_organism_iris": f"cs.{{{taxon_iri}}}",
        "select": select,
        "order": "structure_id",  # offset ile tutarlı sayfalama için zorunlu
    }
    records: list[dict] = []
    offset = 0
    while True:
        body, _cr = _get("epitope_search",
                         {**base, "limit": _PAGE, "offset": offset})
        if body is None:
            return None  # ağ hatası (kısmi veriyi önbelleğe yazma)
        for x in body:
            seq = (x.get("linear_sequence") or "").strip().upper()
            if not seq or not seq.isalpha():
                continue
            records.append({
                "seq": seq,
                "length": len(seq),
                "organisms": x.get("source_organism_names") or [],
                "antigens": x.get("parent_source_antigen_names") or [],
                "reference_ids": x.get("reference_ids") or [],
                "assay_evidence": x.get("mhc_allele_evidences") or [],
                "iedb_id": x.get("structure_id"),
            })
        if len(body) < _PAGE:
            break
        offset += _PAGE

    if resolve_pmids:
        all_refs = [r for rec in records for r in rec["reference_ids"]]
        pmap = _resolve_pmids(all_refs)
        for rec in records:
            rec["pmids"] = sorted({pmap[r] for r in rec["reference_ids"]
                                   if r in pmap})
    else:
        for rec in records:
            rec["pmids"] = []

    cache.write_text(json.dumps(records, ensure_ascii=False), encoding="utf-8")
    return records


# --------------------------------------------------------------------------- #
# Yerel eşleştirme indeksi
# --------------------------------------------------------------------------- #
class Index:
    """Epitop kayıtlarından exact + k-mer örtüşme indeksi."""

    def __init__(self, records: list[dict], k: int = DEFAULT_K):
        self.k = k
        self.records = records
        self.exact: dict[str, list[int]] = {}
        self.kmer: dict[str, set[int]] = {}
        for idx, rec in enumerate(records):
            s = rec["seq"]
            self.exact.setdefault(s, []).append(idx)
            for i in range(len(s) - k + 1):
                self.kmer.setdefault(s[i:i + k], set()).add(idx)

    def match(self, seq: str) -> dict | None:
        """En iyi eşleşmeyi döndür (exact > içerme > çekirdek). Yoksa None."""
        seq = seq.strip().upper()
        k = self.k
        # 1) exact
        if seq in self.exact:
            return self._pack(self.exact[seq][0], "exact", seq)
        # 2/3) k-mer örtüşmesi olan tüm kayıt adayları
        cand: set[int] = set()
        for i in range(len(seq) - k + 1):
            cand |= self.kmer.get(seq[i:i + k], set())
        if not cand:
            return None
        best = None  # (öncelik, örtüşme_uzunluğu, idx, tip)
        for idx in cand:
            es = self.records[idx]["seq"]
            if seq in es:
                pr, ov, tip = 3, len(seq), "aday_epitop_içinde"
            elif es in seq:
                pr, ov, tip = 2, len(es), "epitop_aday_içinde"
            else:
                pr, ov, tip = 1, _shared(seq, es, k), f"ortak_{k}mer"
            key = (pr, ov)
            if best is None or key > best[0]:
                best = (key, idx, tip, ov)
        if best is None:
            return None
        _key, idx, tip, ov = best
        return self._pack(idx, tip, seq, overlap=ov)

    def _pack(self, idx: int, match_type: str, query: str,
              overlap: int | None = None) -> dict:
        rec = self.records[idx]
        return {
            "matched": True,
            "match_type": match_type,
            "overlap": overlap if overlap is not None else len(query),
            "epitope_seq": rec["seq"],
            "iedb_id": rec.get("iedb_id"),
            "organisms": rec.get("organisms", []),
            "antigens": rec.get("antigens", []),
            "pmids": rec.get("pmids", []),
            "reference_ids": rec.get("reference_ids", []),
            "assay_evidence": rec.get("assay_evidence", []),
        }


def _shared(a: str, b: str, k: int) -> int:
    """a ile b arasındaki en uzun ortak alt-dizinin (contiguous) uzunluğu (>=k varsa)."""
    bset = set()
    for i in range(len(b) - k + 1):
        bset.add(b[i:i + k])
    best = 0
    # en uzun ardışık ortak parça: kaba ama yeterli (küçük peptitler)
    for i in range(len(a)):
        for j in range(i + k, len(a) + 1):
            if a[i:j] in b:
                best = max(best, j - i)
    return best or (k if any(a[i:i + k] in bset for i in range(len(a) - k + 1)) else 0)


# --------------------------------------------------------------------------- #
# Canlı yedek: tek peptit, tüm IEDB'ye karşı substring
# --------------------------------------------------------------------------- #
def _live_match(seq: str) -> dict | None:
    seq = seq.strip().upper()
    body, _ = _get("epitope_search", {
        "linear_sequence": f"ilike.*{seq}*",
        "structure_type": "eq.Linear peptide",
        "qualitative_measure": "neq.Negative",
        "select": ("structure_id,linear_sequence,source_organism_names,"
                   "parent_source_antigen_names,reference_ids"),
        "limit": 5,
    })
    if not body:
        return None
    # exact eşleşmeyi tercih et; yoksa en kısa (en spesifik) epitop.
    x = min(body, key=lambda r: (0 if (r.get("linear_sequence") or "").upper() == seq
                                 else 1, len(r.get("linear_sequence") or "")))
    es = (x.get("linear_sequence") or "").upper()
    pmap = _resolve_pmids(x.get("reference_ids") or [])
    return {
        "matched": True,
        "match_type": "exact" if es == seq else "aday_epitop_içinde",
        "overlap": len(seq),
        "epitope_seq": es,
        "iedb_id": x.get("structure_id"),
        "organisms": x.get("source_organism_names") or [],
        "antigens": x.get("parent_source_antigen_names") or [],
        "pmids": sorted(set(pmap.values())),
        "assay_evidence": [],
    }


# --------------------------------------------------------------------------- #
# Ana giriş: adayları IEDB'ye karşı işaretle
# --------------------------------------------------------------------------- #
def annotate_candidates(peptides, taxon_iri: str | None = None,
                        k: int = DEFAULT_K, live_fallback: bool = True) -> dict:
    """Her aday peptidin metrics['iedb'] alanını doldurur. Skorlamayı DEĞİŞTİRMEZ.

    Dönüş: özet dict (kaynak, taxon, kaç aday eşleşti, benchmark varsa).
    """
    records = None
    source = None
    if taxon_iri:
        # Toplu seti PMID çözmeden çek (hızlı); PMID'leri yalnız EŞLEŞEN adaylar
        # için sonradan çözeceğiz.
        records = fetch_source_organism(taxon_iri, resolve_pmids=False)
        source = f"IEDB kaynak-organizma seti ({taxon_iri})"

    n_matched = 0
    if records is not None:
        index = Index(records, k=k)
        hits = []
        for p in peptides:
            hit = index.match(p.seq)
            _apply(p, hit)
            if hit:
                hits.append(hit)
                n_matched += 1
        # eşleşen adayların referanslarını -> PMID (tek toplu sorgu)
        ref_ids = [r for h in hits for r in (h.get("reference_ids") or [])]
        if ref_ids:
            pmap = _resolve_pmids(ref_ids)
            for h in hits:
                h["pmids"] = sorted({pmap[r] for r in (h.get("reference_ids") or [])
                                     if r in pmap})
        summary = {
            "available": True, "source": source, "taxon": taxon_iri,
            "n_records": len(records), "n_candidates": len(peptides),
            "n_matched": n_matched,
            "benchmark": benchmark([p.seq for p in peptides], records, k=k),
        }
        return summary

    # taxon yok → aday-başına CANLI substring sorgusu. Bu O(N) ağ isteğidir
    # (her biri _TIMEOUT'a kadar); çok adayda saatlerce sürer/askıda kalır.
    # Bu yüzden ÜST SINIR: LIVE_MAX'ı aşarsa atla (dürüst not). Toplu/hızlı yol
    # için organizma taxon'u verilmeli (fetch_source_organism → yerel Index).
    if live_fallback and len(peptides) > LIVE_MAX:
        for p in peptides:
            p.metrics["iedb"] = {"matched": None,
                                 "note": f"atlandı ({len(peptides)} aday > {LIVE_MAX}); "
                                         "organizma taxon'u verilirse toplu eşleştirilir"}
            p.methods["iedb"] = METHOD
        return {"available": False, "source": None, "taxon": None,
                "n_candidates": len(peptides), "n_matched": 0, "benchmark": None,
                "note": f"aday-başına canlı IEDB taraması atlandı ({len(peptides)} > "
                        f"{LIVE_MAX} aday, çok yavaş) — organizma taxon'u önerilir"}
    if live_fallback and online():
        for p in peptides:
            hit = _live_match(p.seq)
            _apply(p, hit)
            n_matched += 1 if hit else 0
        return {"available": True, "source": "IEDB canlı substring (tüm DB)",
                "taxon": None, "n_candidates": len(peptides),
                "n_matched": n_matched, "benchmark": None}

    # hiçbir şey yok — dürüst 'çalıştırılmadı'
    for p in peptides:
        p.metrics["iedb"] = {"matched": None, "note": "IEDB erişilemedi (ağ+önbellek yok)"}
    return {"available": False, "source": None, "taxon": taxon_iri,
            "n_candidates": len(peptides), "n_matched": 0, "benchmark": None,
            "note": "IEDB erişilemedi (ağ yok ve yerel önbellek yok)"}


def _apply(p, hit: dict | None):
    if hit:
        p.metrics["iedb"] = hit
        p.methods["iedb"] = METHOD
    else:
        p.metrics["iedb"] = {"matched": False}
        p.methods["iedb"] = METHOD


# --------------------------------------------------------------------------- #
# Recall benchmark (validasyon)
# --------------------------------------------------------------------------- #
def benchmark(predicted_seqs: list[str], records: list[dict],
              k: int = DEFAULT_K, min_len: int = 8, max_len: int = 25) -> dict:
    """Bilinen (IEDB doğrulanmış) epitopların ne kadarını tahminlerimiz yakaladı.

    recall = (örtüşülen benzersiz bilinen epitop) / (toplam benzersiz bilinen epitop)
    precision_benzeri = (bilinen bir epitopla eşleşen tahmin) / (toplam tahmin)

    Örtüşme: exact, içerme (iki yön), veya >=k ortak çekirdek.
    """
    known = sorted({r["seq"] for r in records
                    if min_len <= len(r["seq"]) <= max_len})
    preds = [s.strip().upper() for s in predicted_seqs if s]
    if not known:
        return {"n_known": 0, "recall": None, "note": "IEDB'de uygun epitop yok"}

    # tahminlerden k-mer indeksi
    pred_kmers: dict[str, set[int]] = {}
    for pi, s in enumerate(preds):
        for i in range(len(s) - k + 1):
            pred_kmers.setdefault(s[i:i + k], set()).add(pi)

    def overlaps(ep: str) -> bool:
        for i in range(len(ep) - k + 1):
            if ep[i:i + k] in pred_kmers:
                return True
        # ep, bir tahminin içinde tümüyle yer alıyor mu (kısa ep için)
        return any(ep in s or s in ep for s in preds)

    hit = sum(1 for ep in known if overlaps(ep))
    # tersi: kaç tahmin bilinen bir epitopla eşleşti
    idx = Index(records, k=k)
    pred_matched = sum(1 for s in preds if idx.match(s))
    return {
        "n_known": len(known),
        "n_known_hit": hit,
        "recall": round(hit / len(known), 4) if known else None,
        "n_pred": len(preds),
        "n_pred_matched": pred_matched,
        "precision_like": round(pred_matched / len(preds), 4) if preds else None,
        "k": k,
    }
