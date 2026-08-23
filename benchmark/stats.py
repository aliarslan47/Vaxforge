"""VaxForge benchmark istatistik modülü — yeniden kullanılabilir.

Reverse vaccinology benchmark'ında "bilinen koruyucu antijeni geri bulma" kıyasının
istatistiksel yükünü taşır. Küçük-n gerçeğine dürüst: her fonksiyon n'i raporlar,
n çok küçükse (ör. tek antijen) testi yine hesaplar ama caveat döndürür.

Testler:
  - permutation_rank_test : koruyucu antijenlerin rank'ı rastgeleden iyi mi? (dağılımsız)
  - fisher_enrichment     : 2x2 aday-üretme zenginleşmesi (NERVE2'nin p-değeriyle aynı dil)
  - roc_auc_ci            : candidacy skoru sınıflandırıcı → AUC + bootstrap %95 CI
  - delong_auc_compare    : iki aracın AUC farkı anlamlı mı (eşleşmiş, aynı girdi)
  - jaccard               : iki epitop/host setinin örtüşmesi (çok-konak betimsel)

Atıflar (yayında): Fisher 1922; Mann-Whitney 1947; permütasyon testi (Good 2005);
DeLong 1988 (AUC karşılaştırma); bootstrap CI (Efron 1979).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score


@dataclass
class StatResult:
    name: str
    statistic: float | None
    p_value: float | None
    ci: tuple[float, float] | None = None
    n: int = 0
    detail: dict = field(default_factory=dict)
    caveat: str | None = None

    def as_dict(self) -> dict:
        return {
            "test": self.name, "statistic": self.statistic, "p_value": self.p_value,
            "ci95": list(self.ci) if self.ci else None, "n": self.n,
            "detail": self.detail, "caveat": self.caveat,
        }


def permutation_rank_test(scores: Sequence[float], positive_idx: Sequence[int],
                          n_perm: int = 10000, seed: int = 42,
                          statistic: str = "mean_rank") -> StatResult:
    """Koruyucu antijenler rastgeleden yüksek mi sıralanıyor? Dağılımsız permütasyon.

    scores: tüm proteinlerin candidacy skoru (yüksek=iyi).
    positive_idx: koruyucu antijenlerin scores içindeki indeksleri.
    statistic: 'mean_rank' (ort. rank, düşük=iyi) veya 'best_rank' (en iyi rank).
    Null: pozitifler skorlardan rastgele n_pos protein.
    """
    scores = np.asarray(scores, dtype=float)
    N = len(scores)
    pos = np.asarray(sorted(set(positive_idx)), dtype=int)
    n_pos = len(pos)
    # rank 1 = en yüksek skor
    order = np.argsort(-scores, kind="mergesort")
    rank = np.empty(N, dtype=float)
    rank[order] = np.arange(1, N + 1)

    def stat(idx):
        r = rank[idx]
        return float(np.mean(r)) if statistic == "mean_rank" else float(np.min(r))

    obs = stat(pos)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = stat(rng.choice(N, size=n_pos, replace=False))
    # düşük rank = iyi → gözlenen kadar VEYA daha iyi (küçük) permütasyon oranı
    p = float((np.sum(null <= obs) + 1) / (n_perm + 1))
    caveat = None
    if n_pos < 3 or N < 10:
        caveat = (f"küçük-n (n_pos={n_pos}, N={N}) → test hesaplandı ama güç düşük; "
                  "sonuç betimsel yorumlanmalı")
    return StatResult("permutation_rank_test", obs, p, n=n_pos,
                      detail={"statistic_type": statistic, "N_total": N,
                              "null_mean": float(np.mean(null)),
                              "positive_ranks": rank[pos].tolist()},
                      caveat=caveat)


def fisher_enrichment(prot_candidate: int, prot_total: int,
                      bg_candidate: int, bg_total: int) -> StatResult:
    """2x2 Fisher exact: koruyucular aday-üretmede zenginleşmiş mi? + fold-enrichment."""
    table = [[prot_candidate, prot_total - prot_candidate],
             [bg_candidate, bg_total - bg_candidate]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    p_rate = prot_candidate / max(1, prot_total)
    b_rate = bg_candidate / max(1, bg_total)
    fold = (p_rate / b_rate) if b_rate > 0 else None
    caveat = None
    if prot_total < 3:
        caveat = f"küçük-n (koruyucu={prot_total}) → Fisher hesaplandı ama güç düşük"
    return StatResult("fisher_enrichment", fold, float(p), n=prot_total,
                      detail={"table": table, "odds_ratio": float(odds),
                              "protective_rate": p_rate, "background_rate": b_rate,
                              "fold_enrichment": fold}, caveat=caveat)


def roc_auc_ci(y_true: Sequence[int], y_score: Sequence[float],
               n_boot: int = 2000, seed: int = 42) -> StatResult:
    """candidacy skoru sınıflandırıcı → ROC-AUC + bootstrap %95 CI."""
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    n_pos, n_neg = int(y_true.sum()), int((1 - y_true).sum())
    if n_pos == 0 or n_neg == 0:
        return StatResult("roc_auc_ci", None, None, n=len(y_true),
                          caveat="AUC tanımsız (tek sınıf)")
    auc = float(roc_auc_score(y_true, y_score))
    rng = np.random.default_rng(seed)
    boots = []
    idx = np.arange(len(y_true))
    for _ in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        if y_true[s].sum() in (0, len(s)):
            continue
        boots.append(roc_auc_score(y_true[s], y_score[s]))
    ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))) if boots else None
    caveat = None
    if n_pos < 3:
        caveat = f"küçük-n (pozitif={n_pos}) → AUC/CI gürültülü, betimsel yorumla"
    return StatResult("roc_auc_ci", auc, None, ci=ci, n=len(y_true),
                      detail={"n_pos": n_pos, "n_neg": n_neg}, caveat=caveat)


def delong_auc_compare(y_true: Sequence[int], score_a: Sequence[float],
                       score_b: Sequence[float], n_boot: int = 2000,
                       seed: int = 42) -> StatResult:
    """İki aracın AUC farkı (aynı girdi, eşleşmiş). Bootstrap-tabanlı iki-yanlı p."""
    y_true = np.asarray(y_true, dtype=int)
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    if y_true.sum() in (0, len(y_true)):
        return StatResult("delong_auc_compare", None, None, caveat="tek sınıf")
    auc_a, auc_b = roc_auc_score(y_true, a), roc_auc_score(y_true, b)
    rng = np.random.default_rng(seed)
    idx = np.arange(len(y_true))
    diffs = []
    for _ in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        if y_true[s].sum() in (0, len(s)):
            continue
        diffs.append(roc_auc_score(y_true[s], a[s]) - roc_auc_score(y_true[s], b[s]))
    diffs = np.asarray(diffs)
    # iki-yanlı ampirik p (fark 0'ı içeriyor mu)
    p = float(2 * min((np.mean(diffs <= 0)), (np.mean(diffs >= 0)))) if len(diffs) else None
    ci = (float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))) if len(diffs) else None
    return StatResult("delong_auc_compare", float(auc_a - auc_b), p, ci=ci, n=len(y_true),
                      detail={"auc_a": float(auc_a), "auc_b": float(auc_b)})


def jaccard(set_a: set, set_b: set) -> StatResult:
    """İki setin (ör. host epitop repertuarları) Jaccard örtüşmesi."""
    a, b = set(set_a), set(set_b)
    union = a | b
    j = len(a & b) / len(union) if union else None
    return StatResult("jaccard", j, None, n=len(union),
                      detail={"intersection": len(a & b), "union": len(union),
                              "only_a": len(a - b), "only_b": len(b - a)})


if __name__ == "__main__":
    # kendi kendine test
    sc = [0.9, 0.85, 0.8, 0.75] + list(np.linspace(0.1, 0.7, 60))
    print(permutation_rank_test(sc, [0, 1, 2, 3]).as_dict())
    print(fisher_enrichment(4, 4, 18, 60).as_dict())
    print(roc_auc_ci([1, 1, 1, 1] + [0] * 60, sc).as_dict())
    print(jaccard({"A", "B", "C"}, {"B", "C", "D"}).as_dict())
