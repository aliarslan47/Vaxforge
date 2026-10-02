#!/bin/bash
set -uo pipefail
ROOT=/home/msi-nb/vaxforge; cd "$ROOT"
ST="$ROOT/benchmark/results/VAXIGN_RERUN.status"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$ST"; }
log "=== VAXIGN RE-RUN (bare) BASLADI: brucella mbovis ==="
for b in brucella mbovis; do
  gf=$(python3 -c "import json;print(json.load(open('benchmark/bacteria/${b}.json'))['gram_flag'])")
  [ "$gf" = "p" ] && t="gram+" || t="gram-"
  PREP="$ROOT/data/validation/${b}_vaxign.faa"
  OUT="$ROOT/benchmark/results/${b}/vaxignml"
  mkdir -p "$OUT/_FEATURE/PSORTB"
  log "VAXIGN baslat: $b (t=$t)"
  if podman run --rm -v "$PREP:$PREP" -v "$OUT:$OUT" -v "$OUT/_FEATURE/PSORTB:/tmp/results" \
      docker.io/e4ong1031/vaxign-ml:latest \
      python3.6 VaxignML.py -i "$PREP" -o "$OUT" -t "$t" >>"$OUT/rerun.log" 2>&1; then
    # <bug>_vaxign.result.tsv -> <bug>_prepared.result.tsv (downstream beklediği ad)
    if [ -f "$OUT/${b}_vaxign.result.tsv" ]; then
      cp "$OUT/${b}_vaxign.result.tsv" "$OUT/${b}_prepared.result.tsv"
      log "VAXIGN bitti: $b (tsv satir=$(wc -l < "$OUT/${b}_prepared.result.tsv"))"
    else
      log "VAXIGN HATA: $b — result.tsv olusmadi (bkz $OUT/rerun.log)"
    fi
  else
    log "VAXIGN HATA: $b — konteyner cikis!=0 (bkz $OUT/rerun.log)"
  fi
done
log "=== VAXIGN RE-RUN BITTI ==="
