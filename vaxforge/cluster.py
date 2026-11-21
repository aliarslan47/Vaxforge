"""Örtüşen aday kümeleme.

Epitop tahmini kayan pencereyle çalışır: aynı proteinin 1 aa kaymış komşuları
ayrı peptit olarak çıkar. Bunlar biyolojik olarak TEK bir epitop bölgesidir, ama
huniden ayrı ayrı geçtikleri için aday sayısı şişer. (BoHV-1 koşusunda 703
adayın 484'ü böyle komşuydu; tam kopya yoktu, o yüzden kopya elemesi bunları
yakalamıyordu.)

Bu modül aynı protein + aynı peptit türü içinde, konumları `min_overlap` aa ya
da daha fazla örtüşen peptitleri tek kümede toplar ve kümeyi en yüksek adaylık
puanlı temsilcisiyle gösterir. Eleme DEĞİL, temsil: elenen üyeler temsilcinin
`metrics["cluster_members"]` alanında dizi olarak durur, yani hiçbir bilgi
kaybolmaz ve karar geri alınabilir.

Neden en yüksek puanlı temsilci: küme içindeki peptitler büyük ölçüde aynı
rezidüleri paylaşır, dolayısıyla puan farkları pencere kenarlarından gelir;
en yüksek puanlı pencere o bölgenin en iyi temsilcisidir.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any

from vaxforge.models import Peptide


def _overlap(a: Peptide, b: Peptide) -> int:
    """İki peptidin protein üzerindeki örtüşen rezidü sayısı (negatifse 0)."""
    a1, a2 = a.start, a.start + len(a.seq)
    b1, b2 = b.start, b.start + len(b.seq)
    return max(0, min(a2, b2) - max(a1, b1))


def cluster(
    peptides: list[Peptide],
    min_overlap: int = 4,
    enabled: bool = True,
) -> tuple[list[Peptide], dict[str, Any]]:
    """Örtüşen adayları kümele; (temsilciler, özet) döndür.

    Girdi adaylık puanına göre sıralı gelmek ZORUNDA DEĞİL; burada puana göre
    azalan sırada gezilir, böylece her kümenin temsilcisi kümenin en iyisidir.
    Çıktı, girdideki sıralamayı (puan sırası) korur.

    `enabled=False` ise hiçbir şey yapılmaz ve özet `enabled: False` ile döner —
    eski şişkin davranışı karşılaştırma için geri almak mümkün kalsın.
    """
    summary: dict[str, Any] = {
        "enabled": bool(enabled),
        "min_overlap": int(min_overlap),
        "n_before": len(peptides),
        "n_after": len(peptides),
        "n_absorbed": 0,
        "n_clusters_multi": 0,
        "largest_cluster": 1 if peptides else 0,
    }
    if not enabled or len(peptides) < 2:
        return peptides, summary

    # Puana göre azalan: her kümenin ilk görülen üyesi temsilci olur.
    order = sorted(peptides, key=lambda p: -p.candidacy)
    reps_by_group: dict[tuple[str, str], list[Peptide]] = defaultdict(list)
    absorbed: set[int] = set()
    members: dict[int, list[Peptide]] = defaultdict(list)

    for pep in order:
        group = (pep.parent, pep.kind)
        host = None
        for rep in reps_by_group[group]:
            if _overlap(rep, pep) >= min_overlap:
                host = rep
                break
        if host is None:
            reps_by_group[group].append(pep)
        else:
            absorbed.add(id(pep))
            members[id(host)].append(pep)

    # Temsilcilere küme bilgisini yaz.
    for rep in (p for g in reps_by_group.values() for p in g):
        mem = members.get(id(rep), [])
        rep.metrics["cluster_size"] = len(mem) + 1
        rep.metrics["cluster_members"] = [
            {"peptide": m.seq, "start": m.start, "candidacy": round(m.candidacy, 4)}
            for m in sorted(mem, key=lambda m: m.start)
        ]
        if mem:
            rep.notes.append(
                f"örtüşen {len(mem)} komşu adayı temsil ediyor "
                f"(≥{min_overlap} aa örtüşme, aynı protein)"
            )

    kept = [p for p in peptides if id(p) not in absorbed]
    sizes = [int(p.metrics.get("cluster_size", 1)) for p in kept]
    summary.update(
        n_after=len(kept),
        n_absorbed=len(absorbed),
        n_clusters_multi=sum(1 for s in sizes if s > 1),
        largest_cluster=max(sizes) if sizes else 0,
    )
    return kept, summary
