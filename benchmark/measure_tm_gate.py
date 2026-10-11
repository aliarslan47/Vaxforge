"""TM kapısı (max_helices) ne kadarını eliyor ve bedeli ne? İddia değil, sayı.

Gerekçe: ileriye dönük ground-truth'ta SspA'nın (subtilizin-benzeri S8 peptidaz)
kaybının bağlayıcı kısıtı ÖLÇÜLDÜ ve lokalizasyon değil TM kapısı çıktı
(benchmark/results/funnel_dropout.json: `cell_wall ∈ izinli` ama `2 TM > 1`).
Bu, şimdiye kadarki tek kayıp mekanizmasından (lokalizasyon) farklı, ikinci bir
mekanizma. Soru: tek vaka mı, yoksa kapı Gram(+) yüzey proteinleri sınıfında
sistematik olarak mı yanlış yerde duruyor?

Hipotez (sınanıyor, varsayılmıyor): Gram(+) hücre-duvarına çapalı yüzey proteinleri
N-ucunda sinyal peptidi, C-ucunda LPXTG/membran çapası taşır. Hidrofobik sayıcı bu
ikisini iki ayrı TM sarmalı sayarsa, "çok-TM = zor ifade edilir" gerekçesiyle yazılmış
kapı tam olarak hedeflediğimiz sınıfı eler.

Ölçüm (üç sayı, hiçbiri PSORTb gerektirmiyor → ucuz):
  1) tüm proteomda TM dağılımı;
  2) sinyal-peptidi(+) havuzunda TM>max oranı — kapının yüzey havuzundaki bedeli;
  3) deneysel koruyucu/prospektif antijenlerde TM — kapının bilinen antijenlerdeki
     yanlış-pozitif oranı. Asıl karar sayısı bu.

Koşum: .venv/bin/python3 benchmark/measure_tm_gate.py --bug ssuis
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import signalp, tmhmm_local  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.sequtils import count_tm_helices  # noqa: E402

OUT = ROOT / "benchmark" / "results" / "tm_gate_scope.json"


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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bug", default="ssuis")
    args = ap.parse_args()

    bug_cfg = json.loads((ROOT / "benchmark" / "bacteria" / f"{args.bug}.json")
                         .read_text(encoding="utf-8"))
    proteome = ROOT / "data" / "validation" / f"{args.bug}_prepared.faa"
    if not proteome.exists():
        print(f"proteom yok: {proteome}")
        return 1
    seqs = read_faa(proteome)
    pairs = list(seqs.items())

    cfg = ThresholdConfig.load()
    max_tm = int(cfg.resolve("bacteria")["transmembrane"].params["max_helices"].value)
    sp_min = float(cfg.resolve("bacteria")["signal_peptide"].params["min_probability"].value)

    tm_real = tmhmm_local.predict(pairs) if tmhmm_local.available() else {}
    tm_method = "GERÇEK (tmhmm_local)" if tm_real else "Kyte-Doolittle hidropati (yedek)"
    tm = {a: tm_real.get(a, count_tm_helices(s)) for a, s in pairs}

    sp_real = signalp.predict(pairs, "bacteria") if signalp.available() else {}
    sp_method = "GERÇEK (SignalP)" if sp_real else "heuristik (funnel yedeği yok → atlandı)"
    # signalp.predict -> {id: {sp_prob, has_sp}} (funnel.py:197 ile aynı anahtar)
    sp = {a: float(sp_real.get(a, {}).get("sp_prob", 0.0)) for a, _ in pairs}

    # 1) proteom geneli
    dist = Counter(tm.values())
    n = len(pairs)
    n_over = sum(1 for a in tm if tm[a] > max_tm)

    # 2) sinyal-peptidi(+) havuzu (yüzeye yönelen sınıf)
    sp_pos = [a for a in sp if sp[a] >= sp_min] if sp_real else []
    n_sp_over = sum(1 for a in sp_pos if tm[a] > max_tm)

    # 3) deneysel antijenler: koruyucu liste + prospektif küratörlük
    antigens: dict[str, str] = dict(bug_cfg.get("protective", {}))
    pros = ROOT / "benchmark" / "results" / f"prospective_groundtruth_{args.bug}.json"
    if pros.exists():
        for a in json.loads(pros.read_text(encoding="utf-8"))["antigens"]:
            if a.get("target_accession"):
                antigens.setdefault(a["target_accession"], a["name"])
    ag_rows = []
    for acc, name in sorted(antigens.items()):
        if acc not in tm:
            continue
        ag_rows.append({"accession": acc, "name": name, "length": len(seqs[acc]),
                        "tm_helices": tm[acc], "signalp": round(sp.get(acc, 0.0), 3),
                        "dropped_by_tm_gate": tm[acc] > max_tm})
    n_ag_dropped = sum(1 for r in ag_rows if r["dropped_by_tm_gate"])

    rep = {
        "bug": args.bug, "gram": bug_cfg["gram"],
        "question": "TM kapısı (max_helices) ne kadarını eliyor, bedeli ne?",
        "max_helices": max_tm, "signalp_min_probability": sp_min,
        "tm_method": tm_method, "signalp_method": sp_method,
        "proteome": {
            "n_proteins": n,
            "tm_distribution": {str(k): v for k, v in sorted(dist.items())},
            "n_over_gate": n_over,
            "pct_over_gate": round(100.0 * n_over / n, 1) if n else None,
        },
        "signal_peptide_pool": {
            "n_signalp_positive": len(sp_pos),
            "n_over_gate": n_sp_over,
            "pct_over_gate": round(100.0 * n_sp_over / len(sp_pos), 1) if sp_pos else None,
            "note": ("SignalP yoksa bu bölüm boştur; funnel'ın kendi yedeği heuristik "
                     "olduğu için burada taklit edilmiyor (dürüst boş)"),
        },
        "experimental_antigens": {
            "n": len(ag_rows), "n_dropped_by_tm_gate": n_ag_dropped,
            "pct_dropped": round(100.0 * n_ag_dropped / len(ag_rows), 1) if ag_rows else None,
            "rows": ag_rows,
        },
        "interpretation": (
            "Karar sayısı experimental_antigens.pct_dropped. Tek antijen düşüyorsa kapı "
            "yanlış kalibre DEĞİL, o antijen sınıf-dışı bir vaka olabilir; birden fazlası "
            "düşüyorsa kapı Gram(+) yüzey sınıfında yeniden kalibre edilmeli. Kalibrasyon "
            "yapılmadan kapı DEĞİŞTİRİLMEZ (FILTER_CALIBRATION.md kuralı)."),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"TM yöntemi: {tm_method} · SignalP: {sp_method} · max_helices={max_tm}")
    print(f"proteom: {n} protein · TM dağılımı: "
          f"{dict(sorted(dist.items()))} · kapı üstü: {n_over} (%{rep['proteome']['pct_over_gate']})")
    if sp_pos:
        print(f"sinyal-peptidi(+) havuzu: {len(sp_pos)} · kapı üstü: {n_sp_over} "
              f"(%{rep['signal_peptide_pool']['pct_over_gate']})")
    print(f"\ndeneysel antijenler ({len(ag_rows)}):")
    for r in ag_rows:
        flag = "  <<< TM KAPISI ELER" if r["dropped_by_tm_gate"] else ""
        print(f"  {r['accession']:12} TM={r['tm_helices']} n={r['length']:5} "
              f"{r['name'][:46]:46}{flag}")
    print(f"\nTM kapısının elediği deneysel antijen: {n_ag_dropped}/{len(ag_rows)} "
          f"(%{rep['experimental_antigens']['pct_dropped']})")
    print(f"yazıldı: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
