"""Yüklenen dosyayı otomatik tanıma.

Amaç: kullanıcı bir dosya bıraktığında pipeline'ın 'bu nedir ve ne yapmam
gerekir' sorusuna kendi cevap vermesi.

Tanınan girdiler:
  - FASTQ (ham okumalar)        -> QC + assembly + gen tahmini dalı
  - FASTA / nükleotid (genom)   -> gen tahmini -> proteom dalı
  - FASTA / nükleotid (CDS)     -> çeviri -> proteom dalı
  - FASTA / protein (proteom)   -> doğrudan antijen taraması

Harici bir araca ihtiyaç duymadan, dosyanın ilk kısmını okuyarak karar verir.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path

# Nükleotid alfabesi (belirsizlik kodları dahil)
_NT_CHARS = set("ACGTUNRYSWKMBDHV")
# Sadece proteinlerde görülen (nükleotid olmayan) amino asit harfleri
_PROTEIN_ONLY = set("EFILPQZ*")

_MAX_SCAN_SEQS = 200      # tip kararı için taranacak en fazla kayıt
_MAX_SCAN_BYTES = 2_000_000


@dataclass
class Detection:
    fmt: str = "unknown"            # fasta | fastq | unknown
    seq_type: str = "unknown"       # nucleotide | protein | unknown
    molecule: str = "unknown"       # reads | genome | cds | proteome | unknown
    num_seqs: int = 0
    total_len: int = 0
    avg_len: float = 0.0
    max_len: int = 0
    min_len: int = 0
    is_gzipped: bool = False
    filename: str = ""
    notes: list[str] = field(default_factory=list)
    confident: bool = True

    @property
    def summary(self) -> str:
        return (
            f"{self.fmt.upper()} / {self.seq_type} / {self.molecule} — "
            f"{self.num_seqs} kayıt, ort. {self.avg_len:.0f} bp/aa"
        )


def _open_text(path: Path) -> tuple[io.TextIOBase, bool]:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    gz = magic == b"\x1f\x8b"
    if gz:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace"), True
    return open(path, "rt", encoding="utf-8", errors="replace"), False


def _classify_seq_type(seq: str) -> str:
    """Bir dizinin nükleotid mi protein mi olduğunu içerik oranıyla belirler."""
    s = seq.strip().upper()
    if not s:
        return "unknown"
    if _PROTEIN_ONLY & set(s):
        return "protein"
    nt = sum(1 for c in s if c in _NT_CHARS)
    frac = nt / len(s)
    # Nükleotid dizileri neredeyse tümüyle ACGT(U/N)'dir; proteinler değil.
    return "nucleotide" if frac >= 0.9 else "protein"


def detect(path: str | Path) -> Detection:
    """Bir dosyayı tanır ve Detection döndürür."""
    path = Path(path)
    det = Detection(filename=path.name)
    handle, det.is_gzipped = _open_text(path)
    lengths: list[int] = []
    type_votes = {"nucleotide": 0, "protein": 0}
    try:
        first = ""
        for line in handle:
            if line.strip():
                first = line.strip()
                break
        if not first:
            det.notes.append("Boş dosya.")
            det.confident = False
            return det

        if first.startswith("LOCUS"):
            handle.seek(0)
            _scan_genbank(handle, det)
            return det
        if first.startswith("@"):
            det.fmt = "fastq"
        elif first.startswith(">"):
            det.fmt = "fasta"
        else:
            det.fmt = "unknown"
            det.confident = False
            det.notes.append("İlk satır ne '>' ne '@' ile başlıyor; format tanınamadı.")
            return det

        handle.seek(0)
        if det.fmt == "fasta":
            lengths, type_votes = _scan_fasta(handle)
        else:
            lengths, type_votes = _scan_fastq(handle)
    finally:
        handle.close()

    if not lengths:
        det.confident = False
        det.notes.append("Hiç kayıt okunamadı.")
        return det

    det.num_seqs = len(lengths)
    det.total_len = sum(lengths)
    det.avg_len = det.total_len / len(lengths)
    det.max_len = max(lengths)
    det.min_len = min(lengths)
    det.seq_type = "protein" if type_votes["protein"] > type_votes["nucleotide"] else "nucleotide"

    _infer_molecule(det)
    return det


def _scan_genbank(handle, det: Detection) -> None:
    """GenBank: kayıt (kromozom/kontig) ve CDS özniteliklerini sayar."""
    from Bio import SeqIO
    det.fmt = "genbank"
    det.seq_type = "nucleotide"
    det.molecule = "genbank"
    n_rec = 0
    cds_lens = []
    try:
        for rec in SeqIO.parse(handle, "genbank"):
            n_rec += 1
            for f in rec.features:
                if f.type == "CDS":
                    tr = f.qualifiers.get("translation", [None])[0]
                    cds_lens.append(len(tr) if tr else max(1, len(f.location) // 3))
    except Exception as e:
        det.confident = False
        det.notes.append(f"GenBank ayrıştırma sorunu: {e}")
    det.num_seqs = len(cds_lens)          # CDS sayısı (kullanıcı bunu görmek istiyor)
    if cds_lens:
        det.total_len = sum(cds_lens)
        det.avg_len = det.total_len / len(cds_lens)
        det.max_len = max(cds_lens)
        det.min_len = min(cds_lens)
    if cds_lens:
        det.notes.append(f"{n_rec} GenBank kaydı (kromozom/kontig), {len(cds_lens)} CDS bulundu — "
                         "CDS çevirileri doğrudan kullanılacak (gen/lokus/konum ile).")
    else:
        det.notes.append(f"{n_rec} GenBank kaydı bulundu ama CDS annotasyonu YOK "
                         "(çıplak genom) — ham dizide 6-çerçeve ORF taramasıyla protein çıkarılacak.")


def _scan_fasta(handle) -> tuple[list[int], dict]:
    lengths: list[int] = []
    votes = {"nucleotide": 0, "protein": 0}
    seq_parts: list[str] = []
    read_bytes = 0

    def flush():
        if seq_parts:
            seq = "".join(seq_parts)
            lengths.append(len(seq))
            t = _classify_seq_type(seq[:2000])
            if t in votes:
                votes[t] += 1

    for line in handle:
        read_bytes += len(line)
        if line.startswith(">"):
            flush()
            seq_parts = []
            if len(lengths) >= _MAX_SCAN_SEQS or read_bytes >= _MAX_SCAN_BYTES:
                break
        else:
            seq_parts.append(line.strip())
    flush()
    return lengths, votes


def _scan_fastq(handle) -> tuple[list[int], dict]:
    lengths: list[int] = []
    votes = {"nucleotide": 0, "protein": 0}
    read_bytes = 0
    while True:
        header = handle.readline()
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            break
        read_bytes += len(header) + len(seq) + len(plus) + len(qual)
        s = seq.strip()
        lengths.append(len(s))
        t = _classify_seq_type(s)
        if t in votes:
            votes[t] += 1
        if len(lengths) >= _MAX_SCAN_SEQS or read_bytes >= _MAX_SCAN_BYTES:
            break
    return lengths, votes


def _infer_molecule(det: Detection) -> None:
    """Format + tip + uzunluk profiline göre molekül türünü tahmin eder."""
    if det.fmt == "fastq":
        det.molecule = "reads"
        det.notes.append("Ham okumalar: assembly ve gen tahmini gerekir.")
        return
    if det.seq_type == "protein":
        det.molecule = "proteome"
        det.notes.append("Protein dizileri: doğrudan antijen taramasına girer.")
        return
    # nükleotid FASTA: genom mu, CDS/gen mi?
    if det.num_seqs <= 20 and det.avg_len >= 20000:
        det.molecule = "genome"
        det.notes.append("Az sayıda uzun kontig: genom gibi görünüyor, gen tahmini yapılacak.")
    elif det.avg_len <= 15000:
        det.molecule = "cds"
        det.notes.append("Çok sayıda kısa dizi: CDS/gen gibi, çevrilerek proteine dönüştürülecek.")
    else:
        det.molecule = "genome"
        det.confident = False
        det.notes.append("Belirsiz nükleotid girdisi; varsayılan olarak genom kabul edildi.")
