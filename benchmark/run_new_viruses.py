"""Otonom orkestratör: 7 yeni virüsü SIRAYLA koşar (küçük→büyük, pox sona).

Kullanıcı kuralı: her virüs bitişinde kayıt + her 30 dk heartbeat (checkpoint.py → STATUS.md/PROGRESS.jsonl).
Bir virüs hata verirse ATLANIR (log) → temiz olanlar her hâlükârda teslim edilir.
Her virüs = run_rabies.py (env-parametrize) subprocess; sonuç {name}_summary.json.

Sıra (protein sayısı): getah3 · rvfv4 · veev10 · reovirus11 · jev15 · vaccinia223 · cowpox233.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
sys.path.insert(0, str(BENCH))
import checkpoint as cp  # noqa: E402

GT = json.loads((BENCH / "new_virus_groundtruth.json").read_text())

# 3 temiz virüs (elle) + 4 prepli (json'dan). (slug, taxon, gt_id, gt_name, gt_acc, virus_name)
CLEAN = {
    "rvfv":     ("NCBITaxon:11588", "YP_003848705", "glycoprotein (Gn/Gc)", "YP_003848705.1", "RVFV"),
    "reovirus": ("NCBITaxon:351073", "YP_003199424", "sigma-1", "YP_003199424.1", "Reovirus"),
    "jev":      ("NCBITaxon:11072", "NP_775666", "envelope protein E", "NP_775666.1", "JEV"),
}
# koşu sırası (küçük→büyük). VXF_ONLY env verilirse yalnız o alt-küme koşulur (ör. re-run).
ORDER = ["getah", "rvfv", "veev", "reovirus", "jev", "vaccinia", "cowpox"]
_only = os.environ.get("VXF_ONLY", "").strip()
if _only:
    keep = [s.strip() for s in _only.split(",") if s.strip()]
    ORDER = [s for s in ORDER if s in keep]
NAMES = {"getah": "Getah", "veev": "VEEV", "vaccinia": "Vaccinia", "cowpox": "Cowpox"}


def cfg_for(slug):
    prot = ROOT / "data" / "validation" / f"{slug}_proteome.faa"
    if slug in CLEAN:
        taxon, gt_id, gt_name, gt_acc, vname = CLEAN[slug]
    else:
        g = GT[slug]
        taxon, gt_id, gt_name, gt_acc = g["taxon"], (g.get("gt_id") or "___none___"), g["gt_name"], g.get("gt_acc", "")
        vname = NAMES.get(slug, slug)
    return prot, taxon, gt_id, gt_name, gt_acc, vname


_current = {"slug": "-"}


def heartbeat_loop(stop):
    while not stop.wait(1800):  # her 30 dk
        cp.heartbeat("run", f"{_current['slug']} sürüyor")


def main():
    stop = threading.Event()
    hb = threading.Thread(target=heartbeat_loop, args=(stop,), daemon=True)
    hb.start()
    cp.log_event("orkestratör başladı", f"sıra: {', '.join(ORDER)}", phase="run")

    done = []
    for slug in ORDER:
        _current["slug"] = slug
        prot, taxon, gt_id, gt_name, gt_acc, vname = cfg_for(slug)
        if not prot.exists():
            cp.log_event("ATLANDI", f"{slug}: proteom yok", phase="run")
            continue
        n_prot = sum(1 for l in prot.read_text().splitlines() if l.startswith(">"))
        outdir = BENCH / "results" / slug
        env = dict(os.environ,
                   VIRUS_PROTEOME=str(prot), VIRUS_TAXON=taxon,
                   VIRUS_GT=gt_id, VIRUS_GT_NAME=gt_name, VIRUS_GT_ACC=gt_acc,
                   VIRUS_LABEL=f"{slug}_proteome.faa ({n_prot} protein)",
                   VIRUS_NAME=vname, VIRUS_OUTDIR=str(outdir))
        cp.log_event("virüs başladı", f"{slug} ({n_prot} protein) gt={gt_name}", phase="run")
        t0 = time.time()
        try:
            r = subprocess.run([str(ROOT / ".venv" / "bin" / "python"), str(BENCH / "run_rabies.py")],
                               env=env, cwd=str(ROOT), capture_output=True, text=True, timeout=14 * 3600)
            dt = round(time.time() - t0)
            log = outdir / "run.log"
            outdir.mkdir(parents=True, exist_ok=True)
            log.write_text((r.stdout or "") + "\n---STDERR---\n" + (r.stderr or ""))
            sumf = outdir / f"{vname.lower()}_summary.json"
            if r.returncode == 0 and sumf.exists():
                s = json.loads(sumf.read_text())
                g = s.get("ground_truth", {})
                msg = (f"{slug} BİTTİ {dt}s · gt {gt_name} rank #{g.get('protein_rank')}/"
                       f"{s.get('n_candidate_proteins')} · aday {s.get('n_candidate_proteins')}")
                cp.log_event("virüs bitti", msg, phase="run")
                done.append((slug, dt, g.get("protein_rank"), s.get("n_candidate_proteins")))
            else:
                cp.log_event("virüs HATA", f"{slug} rc={r.returncode} (bkz {log})", phase="run")
        except subprocess.TimeoutExpired:
            cp.log_event("virüs TIMEOUT", f"{slug} 14h aşıldı", phase="run")
        except Exception as e:
            cp.log_event("virüs HATA", f"{slug}: {e}", phase="run")

    stop.set()
    summary = " | ".join(f"{s}:#{rk}/{nc}({dt}s)" for s, dt, rk, nc in done)
    cp.log_event("orkestratör bitti", f"{len(done)}/{len(ORDER)} tamam · {summary}", phase="done")
    print("TÜM KOŞULAR BİTTİ:", summary)


if __name__ == "__main__":
    raise SystemExit(main())
