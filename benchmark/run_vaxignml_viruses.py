"""Vaxign-ML'i 9 virüste koş (rabies zaten var) → virüs kolu head-to-head tamamlansın.

VaxForge çok-konak; Vaxign-ML host'suz tek-skor (protegenicity). NERVE virüste koşamaz (bacterial-only).
Her virüs: rootless podman (run_competitors.sh deseni) → result.tsv → protegenicity desc STABLE sırala →
ground-truth accession rank'ı. Çıktı results/vaxignml_<slug>_summary.json (rabies summary formatı).
Kayıt: checkpoint.py. Bir virüs hata verirse ATLANIR.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
sys.path.insert(0, str(BENCH))
import checkpoint as cp  # noqa: E402

IMG = "docker.io/e4ong1031/vaxign-ml:latest"
# slug -> (ground-truth token [sample içinde substring], gt adı)
GT = {
    "vsv":      ("NP_041715", "glycoprotein G"),
    "emcv":     ("_matp_4 ", "capsid VP1 (1D)"),
    "getah":    ("YP_164439", "structural polyprotein (E2)"),
    "rvfv":     ("YP_003848705", "glycoprotein (Gn/Gc)"),
    "veev":     ("NP_741966", "E2 envelope glycoprotein"),
    "reovirus": ("YP_003199424", "sigma-1"),
    "jev":      ("NP_775666", "envelope protein E"),
    "vaccinia": ("YP_232970", "L1R"),
    "cowpox":   ("NP_619885", "L1R ortholog (CPXV099)"),
}
ORDER = ["vsv", "emcv", "getah", "rvfv", "veev", "reovirus", "jev", "vaccinia", "cowpox"]


def run_container(prep: Path, out: Path) -> Path | None:
    (out / "_FEATURE" / "PSORTB").mkdir(parents=True, exist_ok=True)
    tsv = out / f"{prep.stem}.result.tsv"
    if tsv.exists() and tsv.stat().st_size > 0:
        return tsv  # zaten koşulmuş (vsv smoke)
    cmd = ["podman", "run", "--rm",
           "-v", f"{prep}:{prep}", "-v", f"{out}:{out}",
           "-v", f"{out}/_FEATURE/PSORTB:/tmp/results",
           IMG, "python3.6", "VaxignML.py", "-i", str(prep), "-o", str(out), "-t", "virus"]
    with open(out / "run.log", "w") as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=3 * 3600)
    return tsv if (r.returncode == 0 and tsv.exists()) else None


def parse_rank(tsv: Path, token: str):
    rows = []
    for line in tsv.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 3:
            sample, pred, prot = parts[0], parts[1], parts[2]
            try:
                rows.append((sample, float(pred), float(prot)))
            except ValueError:
                continue
    if not rows:
        return None, None, None, 0
    # protegenicity desc, STABLE (bağlarda girdi sırası korunur — rabies summary'siyle tutarlı)
    order = sorted(range(len(rows)), key=lambda i: -rows[i][2])
    for rank, idx in enumerate(order, 1):
        if token in rows[idx][0]:
            return rank, rows[idx][2], rows[idx][1], len(rows)
    return None, None, None, len(rows)


def main():
    cp.log_event("vaxignml başladı", f"9 virüs: {', '.join(ORDER)}", phase="vaxignml")
    done = []
    for slug in ORDER:
        token, gt_name = GT[slug]
        prep = ROOT / "data" / "validation" / f"{slug}_proteome.faa"
        out = BENCH / "results" / f"vaxignml_{slug}"
        if not prep.exists():
            cp.log_event("ATLANDI", f"{slug}: proteom yok", phase="vaxignml"); continue
        cp.log_event("vaxignml virüs", f"{slug} başladı", phase="vaxignml")
        t0 = time.time()
        try:
            tsv = run_container(prep, out)
            if not tsv:
                cp.log_event("vaxignml HATA", f"{slug}: container/tsv yok", phase="vaxignml"); continue
            rank, prot, pred, n = parse_rank(tsv, token)
            summary = {"tool": "Vaxign-ML", "pathogen": slug, "n_proteins": n,
                       "gt_name": gt_name, "gt_token": token,
                       "gt_rank": rank, "gt_protegenicity": prot, "gt_prediction": pred,
                       "runtime_s": round(time.time() - t0, 1),
                       "note": "virus tipi; host seçimi YOK (VaxForge çok-konak avantajı burada oynamaz). "
                               "protegenicity desc stable sıralama."}
            (BENCH / "results" / f"vaxignml_{slug}_summary.json").write_text(
                json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
            cp.log_event("vaxignml bitti", f"{slug}: {gt_name} rank #{rank}/{n} ({round(time.time()-t0)}s)", phase="vaxignml")
            done.append((slug, rank, n))
        except Exception as e:
            cp.log_event("vaxignml HATA", f"{slug}: {e}", phase="vaxignml")
    summary = " | ".join(f"{s}:#{r}/{n}" for s, r, n in done)
    cp.log_event("vaxignml TAMAM", f"{len(done)}/9 · {summary}", phase="done")
    print("VAXIGN-ML BİTTİ:", summary)


if __name__ == "__main__":
    raise SystemExit(main())
