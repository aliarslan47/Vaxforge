#!/bin/bash
# Bakteri başına rakipler: NERVE 2.0 (-g p|n) + Vaxign-ML (-t bacteria), TAM PROTEOM.
# Kullanım:  bash benchmark/run_competitors.sh <saureus|listeria|salmonella>
# Çıktı: benchmark/results/<bug>/{nerve,vaxignml}/  + competitors.log
set -uo pipefail

BUG="${1:?kullanım: run_competitors.sh <bug>}"
ROOT="/home/msi-nb/vaxforge"
BENCH="$ROOT/benchmark"
CFG="$BENCH/bacteria/${BUG}.json"
RES="$BENCH/results/${BUG}"
mkdir -p "$RES/nerve" "$RES/vaxignml"
LOG="$RES/competitors.log"

GRAM_FLAG=$(python3 -c "import json;print(json.load(open('$CFG'))['gram_flag'])")
# Vaxign-ML -t değeri: gram+ / gram- / virus  (NERVE -g ise p/n)
if [ "$GRAM_FLAG" = "p" ]; then VXG_TYPE="gram+"; else VXG_TYPE="gram-"; fi
echo "=== [$BUG] rakipler başladı $(date) · gram=$GRAM_FLAG · vaxign-t=$VXG_TYPE ===" | tee "$LOG"

# 0) çıplak-aksesyon başlıklı tam proteom hazırla (üç araç da aynı dosyayı kullanır)
PREP="$ROOT/data/validation/${BUG}_prepared.faa"
if [ -n "${SKIP_PREP:-}" ] && [ -f "$PREP" ]; then
  echo "hazırlanmış proteom kullanılıyor (SKIP_PREP): $PREP" | tee -a "$LOG"
else
  python3 "$BENCH/prep_proteome.py" "$BUG" 2>&1 | tee -a "$LOG"
fi

# 1) Vaxign-ML (XGBoost, hızlı). Rootless podman doğrudan (VaxignML.sh sudo istiyordu → bypass).
echo "--- Vaxign-ML $(date) ---" | tee -a "$LOG"
mkdir -p "$RES/vaxignml/_FEATURE/PSORTB"
podman run --rm -v "$PREP:$PREP" -v "$RES/vaxignml:$RES/vaxignml" \
  -v "$RES/vaxignml/_FEATURE/PSORTB:/tmp/results" \
  docker.io/e4ong1031/vaxign-ml:latest \
  python3.6 VaxignML.py -i "$PREP" -o "$RES/vaxignml" -t "$VXG_TYPE" 2>&1 | tee -a "$LOG"
ls -la "$RES/vaxignml"/*.result.tsv 2>&1 | tee -a "$LOG"

# 2) NERVE 2.0 (docker→podman, psortb servisi). Proteom NERVE dizinine kopyalanır (/workdir mount).
echo "--- NERVE 2.0 $(date) ---" | tee -a "$LOG"
# psortb'u kendimiz taze başlat + hazır bekle (eski durmuş konteyner ad çakışmasını önle).
podman rm -f psortb 2>/dev/null | tee -a "$LOG"
[ ! "$(podman network ls | grep nerve-network)" ] && podman network create nerve-network --attachable 2>&1 | tee -a "$LOG"
podman run --rm -p 8080:8080 --network nerve-network --name psortb -d francecosta/psortb_http_api:v0.0.1 2>&1 | tee -a "$LOG"
echo "psortb hazır bekleniyor..." | tee -a "$LOG"
for i in $(seq 1 30); do
  curl -sf http://localhost:8080/ >/dev/null 2>&1 && { echo "psortb hazır ($((i*3))s)" | tee -a "$LOG"; break; }
  sleep 3
done
cp "$PREP" "$BENCH/tools/NERVE/${BUG}_prepared.faa"
( cd "$BENCH/tools/NERVE" && bash NERVE_auto.sh -g "$GRAM_FLAG" -p1 "${BUG}_prepared.faa" -wd "./out_${BUG}/" ) 2>&1 | tee -a "$LOG"
# aday + elenen çıktılarını results altına kopyala
cp "$BENCH/tools/NERVE/out_${BUG}/vaccine_candidates.csv" "$RES/nerve/" 2>&1 | tee -a "$LOG"
cp "$BENCH/tools/NERVE/out_${BUG}/discarded_proteins.csv" "$RES/nerve/" 2>&1 | tee -a "$LOG"

echo "=== [$BUG] rakipler bitti $(date) ===" | tee -a "$LOG"
echo "VaxignML tsv satır: $(wc -l < "$RES/vaxignml/${BUG}_prepared.result.tsv" 2>/dev/null || echo NA)" | tee -a "$LOG"
echo "NERVE aday satır: $(wc -l < "$RES/nerve/vaccine_candidates.csv" 2>/dev/null || echo NA)" | tee -a "$LOG"
