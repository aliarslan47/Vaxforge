#!/bin/bash
# İki-fix'li (PSORTb chunking + netMHCpan paralel-chunk) kodla ÜÇ bakteriyi de
# yeniden koş. Rakipler zaten geçerli → yalnız VaxForge + compare. SIRALI (her koşu
# ~12 çekirdek kullanır → oversubscribe olmasın). Küçük→büyük proteom sırası.
set -uo pipefail
ROOT="/home/msi-nb/vaxforge"; BENCH="$ROOT/benchmark"
source "$ROOT/.venv/bin/activate" 2>/dev/null
echo "=== RE-RUN ALL BAŞLADI $(date) ===" | tee "$BENCH/results/rerun_all.log"
for BUG in saureus listeria salmonella; do
  echo "--- [$BUG] VaxForge re-run (2-fix) $(date) ---" | tee -a "$BENCH/results/rerun_all.log"
  SKIP_PREP=1 BUG="$BUG" python3 -u "$BENCH/run_bacterium.py" > "$BENCH/results/${BUG}_vaxforge_full.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then echo "!! [$BUG] VaxForge rc=$rc — durdu" | tee -a "$BENCH/results/rerun_all.log"; fi
  python3 "$BENCH/compare_bacterium.py" "$BUG" > "$BENCH/results/${BUG}_compare_final.log" 2>&1
  echo "--- [$BUG] BİTTİ $(date): $(grep -m1 Recall@antijen "$BENCH/results/${BUG}_vaxforge_full.log" 2>/dev/null) ---" | tee -a "$BENCH/results/rerun_all.log"
done
echo "=== TÜM RE-RUN BİTTİ $(date) ===" | tee -a "$BENCH/results/rerun_all.log"
