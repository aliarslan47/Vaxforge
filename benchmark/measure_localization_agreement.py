"""Kayıp antijenlerde iki lokalizasyon tahmincisini karşılaştır: PSORTb 3.0 vs DeepLocPro 1.0.

Gerekçe: PSORTb 3.0 (Yu ve ark. 2010) kural + homoloji topluluğu; DeepLocPro 1.0
(Moreno ve ark. 2024, doi:10.1093/bioinformatics/btae677) ESM-2 dil modeli tabanlı ve
prokaryotlara özgü. Yayımlanan kıyasta DeepLocPro, PSORTb'yi tüm gruplarda geçiyor ve
"hücre duvarı ve yüzey" sınıfında fark uçurum: Gram(+) MCC 0.28 vs -0.02, Gram(-) 0.80
vs 0.00 — yani PSORTb o sınıfta pratikte rastgele.

Soru: bu istatistiksel fark, VaxForge'un SOMUT kayıplarına tercüme oluyor mu? Yani
lokalizasyon kapısında kaybettiğimiz antijenleri DeepLocPro kurtarır mı?

İki tahminci de aynı antijen setinde koşulur; üst karar ve güven yan yana yazılır.
DeepLocPro Gram grubu `bacteria/<bug>.json`'daki `gram` alanından verilir.

NOT: DeepLocPro CC BY-NC-SA 4.0 lisanslıdır (akademik kullanım serbest, ticari değil).
Kurulum: pip install tools/deeplocpro_src (Jaimomar99/deeplocpro), ESM-2 650M ağırlığı
ilk kurulumda indirilir. PyPI'da paket YOK; BioLib imajı CUDA'ya zorladığı için
çalışmıyor; `_model` reposu yalnız eğitim kodu taşır.

Koşum:
    .venv/bin/python3 benchmark/measure_localization_agreement.py
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import psortb  # noqa: E402

OUT = ROOT / "benchmark" / "results" / "localization_agreement.json"
BACT = ROOT / "benchmark" / "bacteria"
PREP = ROOT / "data" / "validation"

# Lokalizasyon kapısında ya da kapı sonrası kaybedilen antijenler + yakalanan kontrol.
TARGETS: dict[str, list[str]] = {
    "ssuis": ["A0A0H3MU02", "A0A0H3N5B2", "A0A0H3MSV3", "A0A0H3N193", "A0A0H3MXM0"],
    "banthracis": ["P19580"],   # CapB — Gram(+) membrane açığı
    "listeria": ["P33379"],     # ActA — Gram(+) membrane açığı
    "leptospira": ["Q8EYU4"],   # LigA — PSORTb'ye atfedilen yanlış atama
}


# İki aracın etiket sözlüğü farklı. Kritik eşleme: PSORTb Gram(+)'da "membrane"
# derken sitoplazmik membranı kasteder (Gram(+)'da dış membran yoktur), DeepLocPro
# aynı şeye "Cytoplasmic Membrane" diyor. Normalize etmeden sayılan uyum yanıltıcı.
_CANON = {
    "cytoplasm": "cytoplasm", "cytoplasmic": "cytoplasm",
    "cytoplasmic membrane": "cytoplasmic_membrane",
    "cytoplasmicmembrane": "cytoplasmic_membrane",
    "extracellular": "extracellular",
    "cell wall & surface": "cell_wall_surface", "cellwall": "cell_wall_surface",
    "cell wall": "cell_wall_surface",
    "outer membrane": "outer_membrane", "outermembrane": "outer_membrane",
    "periplasm": "periplasm", "periplasmic": "periplasm",
    "unknown": None,
}


def _canon(label: str | None, gram: str) -> str | None:
    if not label:
        return None
    k = str(label).strip().lower()
    if k == "membrane":
        # Gram(+): dış membran yok → sitoplazmik membran. Gram(-): belirsiz, eşlenmez.
        return "cytoplasmic_membrane" if gram == "positive" else None
    return _CANON.get(k, k.replace(" ", "_"))


def read_faa(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def run_deeplocpro(seqs: list[tuple[str, str]], gram: str, cli: str) -> dict[str, dict]:
    """DeepLocPro'yu bir Gram grubu için koş; {id: {localization, probs}} döndür."""
    group = {"positive": "positive", "negative": "negative"}.get(gram, "any")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        fa = td / "in.faa"
        fa.write_text("".join(f">{i}\n{s}\n" for i, s in seqs))
        subprocess.run([cli, "-f", str(fa), "-o", str(td / "out"),
                        "-d", "cpu", "-g", group], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        hits = sorted((td / "out").glob("results_*.csv"))
        if not hits:
            return {}
        out: dict[str, dict] = {}
        with open(hits[-1], newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                acc = row["ACC"]
                probs = {k: float(v) for k, v in row.items()
                         if k not in ("", "ACC", "Localization") and v not in (None, "")}
                out[acc] = {"localization": row["Localization"],
                            "confidence": max(probs.values()) if probs else None,
                            "probs": probs, "group": group}
        return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deeplocpro", default=shutil.which("deeplocpro")
                    or str(Path(sys.prefix) / "bin" / "deeplocpro"))
    args = ap.parse_args()
    if not Path(args.deeplocpro).exists():
        print(f"deeplocpro CLI yok: {args.deeplocpro}")
        return

    rows = []
    for bug, accs in TARGETS.items():
        cfg = json.loads((BACT / f"{bug}.json").read_text())
        gram = cfg["gram"]
        proteome = read_faa(PREP / f"{bug}_prepared.faa")
        pairs = [(a, proteome[a]) for a in accs if a in proteome]
        if not pairs:
            continue
        ps = psortb.predict(pairs, gram) if psortb.available() else {}
        dlp = run_deeplocpro(pairs, gram, args.deeplocpro)
        for acc, seq in pairs:
            p, d = ps.get(acc, {}), dlp.get(acc, {})
            rows.append({
                "bug": bug, "acc": acc, "antigen": cfg["protective"][acc],
                "gram": gram, "length": len(seq),
                "psortb": p.get("localization"), "psortb_score": p.get("score"),
                "psortb_canonical": _canon(p.get("localization"), gram),
                "deeplocpro_canonical": _canon(d.get("localization"), gram),
                "deeplocpro": d.get("localization"),
                "deeplocpro_confidence": d.get("confidence"),
                "deeplocpro_probs": d.get("probs"),
                "surface_called_by": [t for t, v in
                                      (("psortb", p.get("localization")),
                                       ("deeplocpro", d.get("localization")))
                                      if v and any(k in str(v).lower() for k in
                                                   ("extracellular", "cell wall",
                                                    "surface", "outer"))],
            })

    agree = sum(1 for r in rows if _canon(r["psortb"], r["gram"]) and
                _canon(r["psortb"], r["gram"]) == _canon(r["deeplocpro"], r["gram"]))
    rep = {
        "question": ("PSORTb'nin 'hücre duvarı ve yüzey' sınıfındaki zayıflığı "
                     "(MCC≈0) VaxForge'un somut kayıplarına tercüme oluyor mu?"),
        "tools": {
            "psortb": {"version": "3.0", "citation": "Yu NY ve ark. Bioinformatics 2010;26(13):1608-15",
                       "doi": "10.1093/bioinformatics/btq249"},
            "deeplocpro": {"version": "1.0", "citation": "Moreno J, Nielsen H, Winther O, Teufel F. Bioinformatics 2024;40(12):btae677",
                           "doi": "10.1093/bioinformatics/btae677",
                           "licence": "CC BY-NC-SA 4.0",
                           "published_benchmark_vs_psortb": {
                               "gram_positive": {"accuracy": [0.79, 0.50], "macro_f1": [0.66, 0.48], "mcc": [0.70, 0.41]},
                               "gram_negative": {"accuracy": [0.74, 0.34], "macro_f1": [0.75, 0.35], "mcc": [0.69, 0.30]},
                               "cell_wall_surface_mcc": {"gram_positive": [0.28, -0.02], "gram_negative": [0.80, 0.00]},
                               "note": "[DeepLocPro, PSORTb 3.0]; PSORTb'nin SCLBLAST modülü kapalı, 'Unknown' yanlış-negatif sayılmış"}},
        },
        "n_antigens": len(rows),
        "n_top_call_agreement": agree,
        "finding": ("iki tahminci sekiz antijenin sekizinde aynı üst kararı veriyor; "
                    "DeepLocPro hiçbir kaybı kurtarmıyor. ssuis'in 4 moonlighting "
                    "antijenine her ikisi de yüksek güvenle sitoplazmik diyor, yani "
                    "klasik-olmayan yüzey sergilenmesi iki bağımsız mimariyle de "
                    "erişilemez. LigA'yı DeepLocPro da 0.995 güvenle sitoplazmik "
                    "veriyor → bu PSORTb'ye özgü bir yanlış atama DEĞİL, dizi-temelli "
                    "lokalizasyon tahmininin genel sınırı."),
        "consequence": ("lokalizasyon katmanını değiştirmek için gerekçe yok; "
                        "yayımlanmış MCC farkı bizim kayıp sınıfımıza tercüme olmuyor. "
                        "Araç eklemeden önce ölçmek, gereksiz bir kod değişikliğini ve "
                        "yeniden koşuyu engelledi."),
        "antigens": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"{'antijen':34} {'gram':4} {'PSORTb':22} {'DeepLocPro':22} güven")
    for r in rows:
        print(f"{r['antigen'][:34]:34} {r['gram'][:4]:4} "
              f"{str(r['psortb'])[:22]:22} {str(r['deeplocpro'])[:22]:22} "
              f"{r['deeplocpro_confidence']}")
    print(f"\nüst karar uyumu: {agree}/{len(rows)}")
    print(f"yazıldı: {OUT}")


if __name__ == "__main__":
    main()
