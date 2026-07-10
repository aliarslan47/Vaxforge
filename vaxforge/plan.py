"""Tanıma sonucundan pipeline planı üretir.

Her adım: id, başlık, hangi araçlar, GPU gerektirir mi, uzak-worker'a mı gider,
ve girdiye göre gerekli mi. (Yapısal doğrulama / docking-MD adımları şimdilik
çıkarıldı; odak aday-belirleme. GPU'lu 'deferred' altyapısı korunur, ileride
geri eklenebilir.)
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .detect import Detection


@dataclass
class Step:
    id: str
    title: str
    tools: list[str] = field(default_factory=list)
    gpu: bool = False
    deferred: bool = False        # bu makinede koşamaz (GPU/ağır) -> uzak worker
    status: str = "pending"       # pending | running | done | skipped | deferred
    note: str = ""


# Girdi-öncesi hazırlık adımları (dala göre değişir)
_PREP = {
    "reads": [
        Step("qc", "Okuma kalite kontrolü (QC)", ["fastp"]),
        Step("assembly", "De novo assembly", ["SPAdes"], note="12 GB RAM'de büyük genomda sıkışabilir"),
        Step("gene_pred", "Gen tahmini", ["Prodigal"]),
    ],
    "genome": [
        Step("gene_pred", "Gen tahmini", ["Prodigal"]),
    ],
    "cds": [
        Step("translate", "CDS -> protein çevirisi", ["Biopython"]),
    ],
    "proteome": [],
}

# Girdiden bağımsız ortak omurga
def _core(det: Detection) -> list[Step]:
    return [
        Step("discovery", "Virülans faktörü / hedef keşfi (küratörlü DB)",
             ["diamond", "VFDB", "Victors", "BV-BRC"]),
        Step("funnel", "Antijen eleme hunisi (melez)",
             ["DeepLoc", "SignalP", "DeepTMHMM", "human-BLAST", "VaxiJen"]),
        Step("epitope", "B/T-hücre epitop tahmini (sliding-window)",
             ["BepiPred", "MHCflurry", "NetMHCIIpan"]),
        Step("mhc_panel", "Çok-organizmalı MHC bağlanma haritası",
             ["NetMHCpan", "IPD-MHC paneli"]),
        Step("survival", "Sağ kalım elemesi (alerjenite/toksisite/kapsam)",
             ["AllerTOP", "ToxinPred", "IEDB-coverage"]),
        # NOT: yapısal doğrulama (AlphaFold peptit-MHC) + docking/MD (HADDOCK/
        # GROMACS) adımları ŞİMDİLİK ÇIKARILDI — odak aday-belirlemede. GPU
        # gelince/gerektiğinde geri eklenir.
        Step("scoring", "Adaylık puanı + sıralama", ["VaxForge scorer"]),
        Step("report", "Rapor + veri paketi + HTML panosu",
             ["PDF", "CSV/FASTA/GenBank/PDB", "HTML"]),
    ]


def build_plan(det: Detection, has_gpu: bool = False) -> list[Step]:
    steps: list[Step] = list(_PREP.get(det.molecule, []))
    steps += _core(det)
    if has_gpu:
        for s in steps:
            if s.deferred:
                s.deferred = False
                s.note = s.note.replace("-> ayrı makine/bulut", "(yerel GPU)")
    for s in steps:
        if s.deferred:
            s.status = "deferred"
    return steps


def plan_table(steps: list[Step]) -> list[dict]:
    return [
        {
            "#": i + 1,
            "adım": s.title,
            "araçlar": ", ".join(s.tools),
            "GPU": "evet" if s.gpu else "",
            "durum": s.status,
            "not": s.note,
        }
        for i, s in enumerate(steps)
    ]
