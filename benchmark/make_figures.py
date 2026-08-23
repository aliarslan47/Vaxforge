"""Benchmark figürleri — kaydedilen JSON sonuçlarından üretir (tekrarlanabilir).

Figürler → benchmark/figures/:
  fig1_antigen_recovery.png : araç × 4CMenB antijeni kurtarma matrisi (NERVE 3'ü eliyor)
  fig2_recall_bar.png       : recall@antijen + AUC (204-set, adil ortak girdi)
  fig3_fold_vs_recall.png   : fold-enrichment vs recall ödünleşimi (NERVE yüksek-fold/düşük-recall)
  fig4_multihost_jaccard.png: rabies G-epitop host×host Jaccard (çok-konak, VaxForge'a özgü)
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BENCH = Path(__file__).resolve().parent
FIG = BENCH / "figures"
FIG.mkdir(exist_ok=True)
ANT = ["fHbp", "NadA", "NHBA(VFDB-)", "PorA"]
TOOLS = ["VaxForge", "Vaxign-ML", "NERVE 2.0"]
COL = {"VaxForge": "#1f77b4", "Vaxign-ML": "#2ca02c", "NERVE 2.0": "#d62728"}


def _load():
    h204 = json.load(open(BENCH / "results" / "menb204_headtohead.json"))
    full = json.load(open(BENCH / "results" / "menb_full_competitors.json"))
    rab = json.load(open(BENCH / "results" / "rabies" / "rabies_summary.json"))
    return {t["tool"]: t for t in h204["tools"]}, full, rab


def fig1_recovery(h204):
    """Antijen kurtarma matrisi: yeşil=kurtardı(rank), kırmızı=elendi."""
    M = np.zeros((3, 4))
    labels = [["" for _ in ANT] for _ in TOOLS]
    for i, t in enumerate(TOOLS):
        r = h204[t]["ranks"]
        for j, a in enumerate(ANT):
            v = r[a]
            if isinstance(v, int):
                M[i, j] = 1; labels[i][j] = f"#{v}"
            else:
                M[i, j] = 0; labels[i][j] = "elendi"
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.imshow(M, cmap=plt.cm.RdYlGn, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(ANT)
    ax.set_yticks(range(3)); ax.set_yticklabels(TOOLS)
    for i in range(3):
        for j in range(4):
            ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=10, fontweight="bold")
    ax.set_title("4CMenB koruyucu antijenlerinin kurtarılması (204-set, default)\n"
                 "yeşil=aday listesinde (rank) · kırmızı=elendi", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "fig1_antigen_recovery.png", dpi=150); plt.close(fig)


def fig2_recall(h204):
    fig, ax = plt.subplots(figsize=(6, 4))
    recalls = [int(h204[t]["recall"].split("/")[0]) for t in TOOLS]
    aucs = [h204[t]["auc"] for t in TOOLS]
    x = np.arange(3)
    bars = ax.bar(x, recalls, color=[COL[t] for t in TOOLS], width=0.6)
    ax.set_ylim(0, 4.5); ax.set_ylabel("recall@antijen (4 üzerinden)")
    ax.set_xticks(x); ax.set_xticklabels(TOOLS)
    for b, r, a in zip(bars, recalls, aucs):
        txt = f"{r}/4" + (f"\nAUC {a}" if a is not None else "")
        ax.text(b.get_x() + b.get_width() / 2, r + 0.08, txt, ha="center", va="bottom", fontsize=9)
    ax.set_title("MenB: bilinen 4 koruyucu antijenin geri bulunması (204-set, özdeş girdi)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "fig2_recall_bar.png", dpi=150); plt.close(fig)


def fig3_fold_recall(h204, full):
    """fold-enrichment vs recall — ödünleşim. 204 (dolu) + tam-proteom (içi boş) noktalar.

    recall=%100'de üç nokta (VaxForge-204, Vaxign-ML-204, Vaxign-ML-tam) fold ~3.1-3.6'da
    üst üste biner → etiketler leader-line (ok) ile ayrı köşelere taşınır ki örtüşme olmasın.
    """
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    # her nokta için etiket ofsetini elle ayarla (örtüşmeyi önle); dar kümedekiler leader-line ile
    off204 = {"VaxForge": (10, 12), "Vaxign-ML": (10, -20), "NERVE 2.0": (8, 6)}
    for t in TOOLS:
        rec = int(h204[t]["recall"].split("/")[0]) / 4 * 100
        fold = h204[t]["fold_enrichment"]
        ax.scatter(rec, fold, s=130, color=COL[t], label=f"{t} (204)", zorder=3)
        dx, dy = off204[t]
        arrow = dict(arrowstyle="-", color=COL[t], lw=0.7) if abs(dy) > 8 else None
        ax.annotate(f"{t} (204)", (rec, fold), textcoords="offset points", xytext=(dx, dy),
                    fontsize=8, color=COL[t], arrowprops=arrow,
                    ha="left" if dx >= 0 else "right")
    # tam proteom rakipler (içi boş nokta)
    offfull = {"Vaxign-ML": (10, -38), "NERVE 2.0": (10, -6)}
    for tool, key in [("Vaxign-ML", "vaxignml"), ("NERVE 2.0", "nerve")]:
        rec = int(full[key]["recall"].split("/")[0]) / 4 * 100
        fold = full[key]["fold"]
        if fold:
            ax.scatter(rec, fold, s=110, facecolors="none", edgecolors=COL[tool], linewidths=2, zorder=3)
            dx, dy = offfull[tool]
            ax.annotate(f"{tool} (tam 2003)", (rec, fold), textcoords="offset points", xytext=(dx, dy),
                        fontsize=7, color=COL[tool], ha="left",
                        arrowprops=dict(arrowstyle="-", color=COL[tool], lw=0.7))
    ax.axhline(16.57, ls="--", c="gray", lw=1); ax.text(2, 18.5, "NERVE2 yayınlanmış fold 16.57", fontsize=7, color="gray")
    ax.set_xlim(0, 118); ax.set_ylim(0, 70)
    ax.set_xlabel("recall@antijen (%)"); ax.set_ylabel("fold-enrichment")
    ax.set_title("Ödünleşim: yüksek fold ≠ yüksek recall\n(NERVE aşırı-seçici → yüksek fold ama 3 antijen kaçar)", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig3_fold_vs_recall.png", dpi=150); plt.close(fig)


def fig4_jaccard(rab):
    hosts = rab["hosts"]
    jg = rab["multihost"]["jaccard_overlap_G"]
    n = len(hosts); M = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            key = f"{hosts[i]}~{hosts[j]}"
            v = jg.get(key, {}).get("statistic")
            if v is None:
                key = f"{hosts[j]}~{hosts[i]}"; v = jg.get(key, {}).get("statistic")
            M[i, j] = M[j, i] = v if v is not None else np.nan
    fig, ax = plt.subplots(figsize=(5, 4.2))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(hosts, rotation=30, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(hosts)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center",
                    color="white" if M[i, j] < 0.6 else "black", fontsize=9)
    fig.colorbar(im, label="Jaccard (G-epitop repertuar örtüşmesi)")
    ax.set_title("Rabies glikoprotein G: host×host epitop örtüşmesi\n"
                 "(<1.0 → host seçimi epitopları değiştirir; NERVE/Vaxign yapamaz)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "fig4_multihost_jaccard.png", dpi=150); plt.close(fig)


def main():
    h204, full, rab = _load()
    fig1_recovery(h204); fig2_recall(h204); fig3_fold_recall(h204, full); fig4_jaccard(rab)
    print("Figürler:", *(p.name for p in sorted(FIG.glob("*.png"))))


if __name__ == "__main__":
    main()
