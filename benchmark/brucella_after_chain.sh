#!/bin/bash
# Ana zincir (mbovis+leptospira+banthracis) bitene dek bekle, SONRA brucella VaxForge'u
# DÜZELTİLMİŞ kodla tek başına koştur (netMHCpan deadlock fix sonrası). Çakışma yok.
set -uo pipefail
ROOT=/home/msi-nb/vaxforge; cd "$ROOT"
ST="$ROOT/benchmark/results/BRUCELLA_RERUN.status"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$ST"; }
log "bekliyor: ana zincir bitsin (run_new6_chain)..."
while pgrep -f run_new6_chain >/dev/null 2>&1; do sleep 60; done
log "ana zincir bitti → brucella VaxForge baslat (fixed kod)"
if BUG=brucella .venv/bin/python benchmark/run_bacterium.py >>"$ROOT/benchmark/results/brucella/vaxforge_run.log" 2>&1; then
  r=$(python3 -c "import json;print(json.load(open('benchmark/results/brucella/vaxforge_summary.json'))['recall_at_antigen'])" 2>/dev/null)
  log "brucella VaxForge BITTI (recall=$r)"
else
  log "brucella VaxForge HATA (log: results/brucella/vaxforge_run.log)"
fi
log "=== BRUCELLA RERUN TAMAM ==="
