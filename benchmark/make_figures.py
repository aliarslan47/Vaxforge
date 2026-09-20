"""Benchmark figürleri — kaydedilen JSON sonuçlarından üretir (tekrarlanabilir).

Biçim kuralları (kullanıcı): figürlerde TEPE BAŞLIĞI YOK; her figürde yalnız x ve y ekseni
açıklaması; etiketler İngilizce (yayın standardı). Anlam figür ALTYAZISINDA (makale metni).

Figürler → benchmark/figures/ (numara = makalede ilk-anılma sırası):
  fig1_antigen_recovery.png : antijen × araç RANK nokta-grafiği (✗ = geri getirilmedi)
  fig2_recall_bar.png       : recall@antijen + AUC (204-set)
  fig3_fold_vs_recall.png   : araç başına fold-enrichment (bar) + üstünde recall etiketi
  fig4_bacteria_panel.png   : 3 bakteri × 3 araç recall@antijen + AUC
  fig5_multihost_jaccard.png: rabies G-epitop host-ÇİFTİ Jaccard örtüşmesi (yatay bar)
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
ANT = ["fHbp", "NadA", "NHBA(VFDB-)", "PorA"]     # results JSON anahtarları
ANT_LABEL = ["fHbp", "NadA", "NHBA", "PorA"]       # temiz etiket
TOOLS = ["VaxForge", "Vaxign-ML", "NERVE 2.0"]
COL = {"VaxForge": "#1f77b4", "Vaxign-ML": "#2ca02c", "NERVE 2.0": "#d62728"}
plt.rcParams.update({"font.size": 11})


def _load():
    h204 = json.load(open(BENCH / "results" / "menb204_headtohead.json"))
    full = json.load(open(BENCH / "results" / "menb_full_competitors.json"))
    rab = json.load(open(BENCH / "results" / "rabies" / "rabies_summary.json"))
    return {t["tool"]: t for t in h204["tools"]}, full, rab


def fig1_recovery(h204):
    """Her antijenin (satır) her araçtaki aday-listesi RANK'ı (nokta). ✗ = geri getirilmedi.
    Düşük rank = daha iyi. Araçlar küçük dikey ofsetle ayrılır (üst üste binmesin)."""
    max_rank = max(v for t in TOOLS for v in h204[t]["ranks"].values() if isinstance(v, int))
    elim_x = max_rank + 5
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    offs = {"VaxForge": 0.22, "Vaxign-ML": 0.0, "NERVE 2.0": -0.22}
    for t in TOOLS:
        r = h204[t]["ranks"]
        xr, yr, ye = [], [], []
        for j, a in enumerate(ANT):
            v = r[a]
            if isinstance(v, int):
                xr.append(v); yr.append(j + offs[t])
            else:
                ye.append(j + offs[t])
        ax.scatter(xr, yr, s=95, color=COL[t], label=t, zorder=3)
        for xx, yy in zip(xr, yr):
            ax.annotate(str(int(xx)), (xx, yy), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=7.5, color=COL[t])
        ax.scatter([elim_x] * len(ye), ye, s=95, marker="x", color=COL[t], zorder=3)
    ax.axvline(elim_x - 2.5, ls=":", color="gray", lw=1)
    ax.set_yticks(range(len(ANT_LABEL))); ax.set_yticklabels(ANT_LABEL)
    ax.set_ylim(-0.6, len(ANT_LABEL) - 0.4)
    ax.set_xlim(0, elim_x + 2)
    ticks = list(range(0, max_rank + 1, 5))
    ax.set_xticks(ticks + [elim_x]); ax.set_xticklabels([str(v) for v in ticks] + ["not\nrecovered"])
    ax.set_xlabel("Rank in candidate list (lower is better)")
    ax.set_ylabel("4CMenB protective antigen")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", alpha=0.3)
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
    """Araç başına fold-enrichment (bar). Her barın üstünde o koşunun recall'ı yazılı → yüksek
    fold'un yüksek recall demek OLMADIĞI doğrudan görünür (NERVE: dev fold, recall 1/4)."""
    labels, folds, recalls, colors = [], [], [], []
    for t in TOOLS:
        labels.append(t + "\n(204-set)")
        folds.append(h204[t]["fold_enrichment"])
        recalls.append(h204[t]["recall"]); colors.append(COL[t])
    for tool, key in [("Vaxign-ML", "vaxignml"), ("NERVE 2.0", "nerve")]:
        f = full[key]["fold"]
        if f:
            labels.append(tool + "\n(full proteome)")
            folds.append(f); recalls.append(full[key]["recall"]); colors.append(COL[tool])
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    bars = ax.bar(x, folds, color=colors, width=0.62)
    for b, rec in zip(bars, recalls):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.0, f"recall {rec}",
                ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.axhline(16.57, ls="--", c="gray", lw=1)
    ax.text(len(labels) - 0.4, 18.0, "NERVE 2.0 published fold 16.6", fontsize=7.5,
            color="gray", ha="right")
    ax.set_ylim(0, 72)
    ax.set_ylabel("Fold enrichment of known antigens among candidates")
    ax.set_xlabel("Tool (input set)")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig3_fold_vs_recall.png", dpi=150); plt.close(fig)


def fig4_bacteria_panel():
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
    fig.tight_layout(); fig.savefig(FIG / "fig4_bacteria_panel.png", dpi=150); plt.close(fig)


def fig5_jaccard(rab):
    """Rabies glikoprotein G: her HOST ÇİFTİ için epitop örtüşmesi (Jaccard), yatay bar.
    Simetrik host×host matris yerine okunur çift-bazlı gösterim."""
    hosts = rab["hosts"]
    jg = rab["multihost"]["jaccard_overlap_G"]
    pairs, vals = [], []
    for i in range(len(hosts)):
        for j in range(i + 1, len(hosts)):
            v = jg.get(f"{hosts[i]}~{hosts[j]}", {}).get("statistic")
            if v is None:
                v = jg.get(f"{hosts[j]}~{hosts[i]}", {}).get("statistic")
            if v is not None:
                pairs.append(f"{hosts[i]}–{hosts[j]}"); vals.append(v)
    order = np.argsort(vals)
    pairs = [pairs[k] for k in order]; vals = [vals[k] for k in order]
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    ax.barh(range(len(pairs)), vals, color="#4c72b0", height=0.6)
    for i, v in enumerate(vals):
        ax.text(v + 0.015, i, f"{v:.2f}", va="center", fontsize=9)
    ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs)
    ax.set_xlim(0, 1.08)
    ax.axvline(1.0, ls=":", color="gray", lw=1)
    ax.set_xlabel("Jaccard overlap of glycoprotein G epitopes (1 = identical repertoire)")
    ax.set_ylabel("Host pair")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig5_multihost_jaccard.png", dpi=150); plt.close(fig)


def main():
    h204, full, rab = _load()
    fig1_recovery(h204); fig2_recall(h204); fig3_fold_recall(h204, full)
    fig4_bacteria_panel(); fig5_jaccard(rab)
    # eski (yeniden numaralanmış) dosyaları temizle
    for old in ["fig4_multihost_jaccard.png", "fig5_bacteria_panel.png"]:
        p = FIG / old
        if p.exists():
            p.unlink()
    print("Figürler:", *(p.name for p in sorted(FIG.glob("*.png"))))


if __name__ == "__main__":
    main()
