#!/bin/bash
set -uo pipefail
ROOT=/home/msi-nb/vaxforge; cd "$ROOT"
ST="$ROOT/benchmark/results/VACCINE_CHAIN.status"
BUGS=(brucella_rb51 brucella_s19 brucella_s2)
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$ST"; }
log "=== ASI SUSU ZINCIRI BASLADI (rb51/s19/s2) ==="
for b in "${BUGS[@]}"; do
  log "RAKIP baslat: $b"
  if SKIP_PREP=1 bash benchmark/run_competitors.sh "$b" >>"benchmark/results/${b}_competitors.log" 2>&1; then
    log "RAKIP bitti: $b"
  else log "RAKIP HATA: $b"; fi
done
log "=== FAZ 1 (rakipler) BITTI ==="
for b in "${BUGS[@]}"; do
  mkdir -p "benchmark/results/${b}"
  log "VAXFORGE baslat: $b"
  if BUG="$b" .venv/bin/python benchmark/run_bacterium.py >>"benchmark/results/${b}/vaxforge_run.log" 2>&1; then
    r=$(python3 -c "import json;print(json.load(open('benchmark/results/${b}/vaxforge_summary.json'))['recall_at_antigen'])" 2>/dev/null)
    log "VAXFORGE bitti: $b (recall=$r)"
  else log "VAXFORGE HATA: $b"; fi
done
log "=== ASI SUSU ZINCIRI TAMAMEN BITTI ==="
