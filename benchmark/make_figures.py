"""Benchmark figürleri — kaydedilen JSON sonuçlarından üretir (tekrarlanabilir).

Biçim kuralları (kullanıcı): figürlerde TEPE BAŞLIĞI YOK; her figürde yalnız x ve y ekseni
açıklaması bulunur; etiketler İngilizce (yayın standardı). Açıklama/anlam figür ALTYAZISINDA
(makale metni) verilir, figürün üstünde değil.

Figürler → benchmark/figures/:
  fig1_antigen_recovery.png : araç × 4CMenB antijeni kurtarma matrisi
  fig2_recall_bar.png       : recall@antijen + AUC (204-set)
  fig3_fold_vs_recall.png   : fold-enrichment vs recall ödünleşimi
  fig4_multihost_jaccard.png: rabies G-epitop host×host Jaccard
  fig5_bacteria_panel.png   : 3 bakteri × 3 araç recall@antijen + AUC
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.lines as mlines  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
import numpy as np  # noqa: E402

BENCH = Path(__file__).resolve().parent
FIG = BENCH / "figures"
FIG.mkdir(exist_ok=True)
ANT = ["fHbp", "NadA", "NHBA(VFDB-)", "PorA"]     # results JSON anahtarları
ANT_LABEL = ["fHbp", "NadA", "NHBA", "PorA"]       # figürde gösterilen temiz etiket
TOOLS = ["VaxForge", "Vaxign-ML", "NERVE 2.0"]
COL = {"VaxForge": "#1f77b4", "Vaxign-ML": "#2ca02c", "NERVE 2.0": "#d62728"}
plt.rcParams.update({"font.size": 11, "axes.titlesize": 11})


def _load():
    h204 = json.load(open(BENCH / "results" / "menb204_headtohead.json"))
    full = json.load(open(BENCH / "results" / "menb_full_competitors.json"))
    rab = json.load(open(BENCH / "results" / "rabies" / "rabies_summary.json"))
    return {t["tool"]: t for t in h204["tools"]}, full, rab


def fig1_recovery(h204):
    """Antijen kurtarma matrisi: yeşil hücre = kurtarıldı (rank yazılı), kırmızı = elendi (✗)."""
    M = np.zeros((3, 4))
    labels = [["" for _ in ANT] for _ in TOOLS]
    for i, t in enumerate(TOOLS):
        r = h204[t]["ranks"]
        for j, a in enumerate(ANT):
            v = r[a]
            if isinstance(v, int):
                M[i, j] = 1; labels[i][j] = f"#{v}"
            else:
                M[i, j] = 0; labels[i][j] = "✗"
    cmap = ListedColormap(["#ef9a9a", "#a5d6a7"])
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(ANT_LABEL)
    ax.set_yticks(range(3)); ax.set_yticklabels(TOOLS)
    for i in range(3):
        for j in range(4):
            ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=13,
                    fontweight="bold", color="#1b5e20" if M[i, j] else "#b71c1c")
    ax.set_xlabel("4CMenB protective antigen")
    ax.set_ylabel("Reverse vaccinology tool")
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=3)
    ax.tick_params(which="minor", length=0)
    rec = mpatches.Patch(color="#a5d6a7", label="Recovered (rank shown)")
    elim = mpatches.Patch(color="#ef9a9a", label="Eliminated (✗)")
    ax.legend(handles=[rec, elim], loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "fig1_antigen_recovery.png", dpi=150); plt.close(fig)


def fig2_recall(h204):
    fig, ax = plt.subplots(figsize=(6, 4.2))
    recalls = [int(h204[t]["recall"].split("/")[0]) for t in TOOLS]
    aucs = [h204[t]["auc"] for t in TOOLS]
    x = np.arange(3)
    bars = ax.bar(x, recalls, color=[COL[t] for t in TOOLS], width=0.6)
    ax.set_ylim(0, 4.6)
    ax.set_ylabel("Recall over known antigens (of 4)")
    ax.set_xlabel("Reverse vaccinology tool")
    ax.set_xticks(x); ax.set_xticklabels(TOOLS)
    for b, r, a in zip(bars, recalls, aucs):
        txt = f"{r}/4" + (f"\nAUC {a}" if a is not None else "")
        ax.text(b.get_x() + b.get_width() / 2, r + 0.08, txt, ha="center", va="bottom", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig2_recall_bar.png", dpi=150); plt.close(fig)


def fig3_fold_recall(h204, full):
    """fold-enrichment vs recall ödünleşimi. Renk=araç, dolu daire=204-set, içi boş kare=tam proteom.
    Yüksek recall bandı (sağ) hafif yeşille işaretli; yayınlanmış NERVE fold referans çizgisi."""
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.axvspan(90, 105, color="#a5d6a7", alpha=0.15, zorder=0)
    pts = []
    for t in TOOLS:
        rec = int(h204[t]["recall"].split("/")[0]) / 4 * 100
        pts.append((rec, h204[t]["fold_enrichment"], t, True))
    for tool, key in [("Vaxign-ML", "vaxignml"), ("NERVE 2.0", "nerve")]:
        f = full[key]["fold"]
        if f:
            pts.append((int(full[key]["recall"].split("/")[0]) / 4 * 100, f, tool, False))

    def draw(a, s=170):
        for rec, fold, t, filled in pts:
            a.scatter(rec, fold, s=s, marker=("o" if filled else "s"),
                      color=(COL[t] if filled else "white"),
                      edgecolors=COL[t], linewidths=1.8, zorder=3)

    draw(ax)
    ax.axhline(16.57, ls="--", c="gray", lw=1)
    ax.text(2, 15.0, "NERVE 2.0 published fold = 16.6", fontsize=8, color="gray", va="top")
    ax.set_xlim(0, 108); ax.set_ylim(0, 70)
    ax.set_xlabel("Recall over known antigens (%)")
    ax.set_ylabel("Fold enrichment")
    # tam-recall kümesi üst üste biner → yakınlaştırma-inset ile ayır
    axin = ax.inset_axes([0.40, 0.34, 0.34, 0.40])
    draw(axin, s=110)
    axin.set_xlim(99.0, 101.0); axin.set_ylim(2.7, 4.0)
    axin.set_xticks([100]); axin.set_yticks([3.0, 3.5, 4.0])
    axin.tick_params(labelsize=7)
    axin.grid(alpha=0.3)
    for rec, fold, t, filled in pts:
        if rec >= 90:
            axin.annotate(t + (" (204)" if filled else " (full)"), (rec, fold),
                          textcoords="offset points", xytext=(7, 0), fontsize=6.5,
                          va="center", color=COL[t])
    ax.indicate_inset_zoom(axin, edgecolor="gray")
    tool_h = [mlines.Line2D([], [], marker="o", linestyle="None", markerfacecolor=COL[t],
                            markeredgecolor=COL[t], markersize=10, label=t) for t in TOOLS]
    style_h = [mlines.Line2D([], [], marker="o", linestyle="None", markerfacecolor="gray",
                             markeredgecolor="gray", markersize=9, label="204-protein set"),
               mlines.Line2D([], [], marker="s", linestyle="None", markerfacecolor="white",
                             markeredgecolor="gray", markersize=9, label="full proteome")]
    ax.legend(handles=tool_h + style_h, fontsize=8, loc="upper right", framealpha=0.95)
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
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(hosts, rotation=30, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(hosts)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center",
                    color="white" if M[i, j] < 0.6 else "black", fontsize=10)
    ax.set_xlabel("Host")
    ax.set_ylabel("Host")
    fig.colorbar(im, ax=ax, label="Jaccard overlap of glycoprotein G epitopes")
    fig.tight_layout(); fig.savefig(FIG / "fig4_multihost_jaccard.png", dpi=150); plt.close(fig)


def fig5_bacteria_panel():
    bugs = [("saureus", "S. aureus"), ("listeria", "L. monocytogenes"),
            ("salmonella", "S. Typhimurium")]
    data = {}
    for key, _ in bugs:
        h = json.load(open(BENCH / "results" / f"{key}_headtohead.json"))
        data[key] = {t["tool"]: t for t in h["tools"]}
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    x = np.arange(len(bugs)); w = 0.26
    for k, t in enumerate(TOOLS):
        vals = [int(data[key][t]["recall"].split("/")[0]) for key, _ in bugs]
        bars = ax.bar(x + (k - 1) * w, vals, width=w, color=COL[t], label=t)
        for b, (key, _) in zip(bars, bugs):
            tool = data[key][t]
            rec = tool["recall"]; auc = tool.get("auc")
            lbl = rec + (f"\nAUC {auc:.2f}" if isinstance(auc, (int, float)) else "")
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.06, lbl,
                    ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x); ax.set_xticklabels([n for _, n in bugs], fontsize=10, style="italic")
    ax.set_ylim(0, 6)
    ax.set_ylabel("Recall over known antigens (of 5)")
    ax.set_xlabel("Bacterial species (full proteome)")
    ax.legend(loc="upper right", fontsize=8, ncol=3)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig5_bacteria_panel.png", dpi=150); plt.close(fig)


def main():
    h204, full, rab = _load()
    fig1_recovery(h204); fig2_recall(h204); fig3_fold_recall(h204, full); fig4_jaccard(rab)
    fig5_bacteria_panel()
    print("Figürler:", *(p.name for p in sorted(FIG.glob("*.png"))))


if __name__ == "__main__":
    main()
