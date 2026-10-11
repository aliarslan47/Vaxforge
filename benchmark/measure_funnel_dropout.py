"""Belirli proteinler huninin HANGİ kapısında düşüyor? Tahmin etme, kapı kapı ölç.

Gerekçe: ileriye dönük ground-truth'ta (PROSPECTIVE_GROUNDTRUTH.md) SspA hedef
proteomda %100 kimlikte var, PSORTb ona `cell_wall 9.98` diyor — yani lokalizasyon
kapısı için İZİNLİ bir yer — ama koşuda aday üretmemiş. Şimdiye kadar bulduğumuz
kayıp mekanizmasının (lokalizasyon) dışında bir şey oluyor. Hangi kapı olduğunu
iddia etmek yerine huniyi o protein için koşup trace'i okuyoruz.

Alt-küme güvenliği: funnel'ın tüm sert kapıları protein-başına MUTLAK eşiklerle
çalışıyor (lokalizasyon ∈ izinli-liste, TM ≤ max_helices, self-homoloji < hard eşik);
popülasyona göre göreli bir kapı yok. Bu yüzden iki proteinle koşmak tam proteomla
koşmakla aynı kararı verir. Kontrol proteini (koşuda sıralanmış bir antijen) bu
varsayımı her koşuda sınar: o GEÇMELİ.

Koşum:
    .venv/bin/python3 benchmark/measure_funnel_dropout.py --bug ssuis \
        --acc A0A0H3N6J4 --control A0A0H3MU02
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import funnel, ingest  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402
from vaxforge.detect import detect  # noqa: E402
from vaxforge.hosts import HostRegistry  # noqa: E402

OUT = ROOT / "benchmark" / "results" / "funnel_dropout.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bug", default="ssuis")
    ap.add_argument("--acc", action="append", required=True,
                    help="ölçülecek aksesyon (birden çok kez verilebilir)")
    ap.add_argument("--control", action="append", default=[],
                    help="koşuda aday üretmiş aksesyon — alt-küme varsayımını sınar")
    ap.add_argument("--host", default="human")
    args = ap.parse_args()

    cfg_path = ROOT / "benchmark" / "bacteria" / f"{args.bug}.json"
    bug_cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    gram = bug_cfg["gram"]

    proteome = ROOT / "data" / "validation" / f"{args.bug}_prepared.faa"
    if not proteome.exists():
        print(f"proteom yok: {proteome}")
        return 1

    det = detect(str(proteome))
    wanted = list(dict.fromkeys(args.acc + args.control))
    records = [p for p in ingest.load_proteins(str(proteome), det) if p.id in wanted]
    found = {p.id for p in records}
    for a in wanted:
        if a not in found:
            print(f"proteomda yok: {a}")

    cfg = ThresholdConfig.load()
    tools = cfg.resolve("bacteria")
    hosts = HostRegistry.load()
    host_objs = [h for h in hosts.all() if h.name == args.host] if hasattr(hosts, "all") \
        else [hosts.get(args.host)]

    survivors, summary = funnel.run(records, tools, profile="bacteria",
                                    hosts=host_objs, gram=gram)
    surv_ids = {p.id for p in survivors}

    rows = []
    for p in records:
        gates = [{"step": t.get("step"), "ok": bool(t.get("ok")),
                  "reason": t.get("reason")} for t in p.trace
                 if str(t.get("step", "")).startswith("funnel:")]
        blocking = [g for g in gates if not g["ok"]]
        rows.append({
            "accession": p.id,
            "role": "kontrol" if p.id in (args.control or []) else "ölçülen",
            "passed_funnel": p.id in surv_ids,
            "blocking_gates": [g["step"] for g in blocking],
            "gates": gates,
            "localization": p.annotations.get("localization"),
            "self_homology": p.annotations.get("self_homology"),
            "antigenicity": p.annotations.get("antigenicity"),
            "funnel_score": p.annotations.get("funnel_score"),
        })

    rep = {
        "bug": args.bug, "gram": gram, "host": args.host,
        "question": "verilen proteinler huninin hangi kapısında düşüyor?",
        "subset_safety": ("funnel'ın sert kapıları protein-başına mutlak eşikler; "
                          "göreli kapı yok. Kontrol proteini varsayımı sınar."),
        "allowed_localizations_key": ("allowed_gram_positive" if gram == "positive"
                                      else "allowed"),
        "n_input": len(records), "n_passed": len(surv_ids),
        "funnel_summary_keys": sorted(summary.keys()),
        "proteins": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    for r in rows:
        print(f"\n=== {r['accession']} ({r['role']}) — huni: "
              f"{'GEÇTİ' if r['passed_funnel'] else 'DÜŞTÜ'}")
        print(f"    lokalizasyon={r['localization']} self-homoloji={r['self_homology']} "
              f"antijenite={r['antigenicity']} skor={r['funnel_score']}")
        for g in r["gates"]:
            print(f"    [{'ok ' if g['ok'] else 'DÜŞ'}] {str(g['step']):28} {g['reason']}")
    print(f"\nyazıldı: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
