#!/bin/bash
# Bir bakteri için TAM kol: prep (bir kez) → rakipler + VaxForge EŞZAMANLI → compare.
# Kullanım:  bash benchmark/run_bug_full.sh <saureus|listeria|salmonella>
set -uo pipefail
BUG="${1:?kullanım: run_bug_full.sh <bug>}"
ROOT="/home/msi-nb/vaxforge"
BENCH="$ROOT/benchmark"
RES="$BENCH/results/${BUG}"
mkdir -p "$RES"
source "$ROOT/.venv/bin/activate" 2>/dev/null

echo "=== [$BUG] KOL BAŞLADI $(date) ===" | tee "$RES/orchestrate.log"

# 0) prep bir kez (yarışı önle) → sonra iki iş SKIP_PREP ile paylaşır
python3 "$BENCH/prep_proteome.py" "$BUG" 2>&1 | tee -a "$RES/orchestrate.log"
export SKIP_PREP=1

# 1) rakipler + VaxForge EŞZAMANLI (bağımsızlar; compare ikisini de bekler)
echo "--- rakipler + VaxForge eşzamanlı başlatılıyor ---" | tee -a "$RES/orchestrate.log"
bash "$BENCH/run_competitors.sh" "$BUG" > "$RES/../${BUG}_competitors_run.log" 2>&1 &
CPID=$!
BUG="$BUG" python3 -u "$BENCH/run_bacterium.py" > "$RES/../${BUG}_vaxforge_full.log" 2>&1 &
VPID=$!
echo "competitors PID=$CPID · vaxforge PID=$VPID" | tee -a "$RES/orchestrate.log"
wait $CPID; CRC=$?
echo "rakipler bitti rc=$CRC $(date)" | tee -a "$RES/orchestrate.log"
wait $VPID; VRC=$?
echo "vaxforge bitti rc=$VRC $(date)" | tee -a "$RES/orchestrate.log"

# 2) head-to-head (ikisi de bitti)
if [ -f "$RES/vaxforge_summary.json" ] && ls "$RES/vaxignml/"*.result.tsv >/dev/null 2>&1 && [ -f "$RES/nerve/vaccine_candidates.csv" ]; then
  echo "--- compare_bacterium ---" | tee -a "$RES/orchestrate.log"
  python3 "$BENCH/compare_bacterium.py" "$BUG" 2>&1 | tee -a "$RES/orchestrate.log"
  echo "=== [$BUG] KOL BİTTİ $(date) → results/${BUG}_headtohead.json ===" | tee -a "$RES/orchestrate.log"
else
  echo "=== [$BUG] EKSİK ÇIKTI — compare atlandı (vaxforge_summary/vaxignml tsv/nerve csv kontrol et) ===" | tee -a "$RES/orchestrate.log"
fi
