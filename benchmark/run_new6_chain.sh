#!/bin/bash
# Yeni 6 bakteri gece-boyu zinciri: ÖNCE rakipler (Vaxign-ML+NERVE2), SONRA VaxForge.
# Genom sırası: coxiella ssuis brucella mbovis leptospira banthracis.
# İlerleme: benchmark/results/NEW6_CHAIN.status (append-only) — kopmaya dayanıklı.
set -uo pipefail
ROOT=/home/msi-nb/vaxforge
cd "$ROOT"
ST="$ROOT/benchmark/results/NEW6_CHAIN.status"
ORDER=(coxiella ssuis brucella mbovis leptospira banthracis)
COMP_TODO=(ssuis brucella mbovis leptospira banthracis)   # coxiella rakip zaten bitti
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$ST"; }

log "=== ZINCIR BASLADI ==="
# FAZ 1 — rakipler
for b in "${COMP_TODO[@]}"; do
  log "RAKIP baslat: $b"
  if SKIP_PREP=1 bash benchmark/run_competitors.sh "$b" >>"$ROOT/benchmark/results/${b}/competitors_run.log" 2>&1; then
    n=$(wc -l < "$ROOT/benchmark/results/${b}/vaxignml/${b}_prepared.result.tsv" 2>/dev/null || echo NA)
    c=$(wc -l < "$ROOT/benchmark/results/${b}/nerve/vaccine_candidates.csv" 2>/dev/null || echo NA)
    log "RAKIP bitti: $b (vaxign tsv=$n, nerve aday=$c)"
  else
    log "RAKIP HATA: $b (log: results/${b}/competitors_run.log) — devam"
  fi
done
log "=== FAZ 1 (rakipler) BITTI ==="

# FAZ 2 — VaxForge
for b in "${ORDER[@]}"; do
  log "VAXFORGE baslat: $b"
  if BUG="$b" .venv/bin/python benchmark/run_bacterium.py >>"$ROOT/benchmark/results/${b}/vaxforge_run.log" 2>&1; then
    log "VAXFORGE bitti: $b"
  else
    log "VAXFORGE HATA: $b (log: results/${b}/vaxforge_run.log) — devam"
  fi
done
log "=== ZINCIR TAMAMEN BITTI ==="
