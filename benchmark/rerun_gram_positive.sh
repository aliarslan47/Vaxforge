#!/bin/bash
# Gram(+) lokalizasyon düzeltmesinden sonra YALNIZ Gram-pozitif dalı yeniden koş.
#
# Gerekçe: config/thresholds.yaml'a localization.allowed_gram_positive eklendi
# ([extracellular, cell_wall, membrane]); Gram(+)'de dış zar yoktur, sitoplazmik zar
# hücrenin dış sınırıdır ve yüzeye çıkan proteinler PSORTb'de 'Cytoplasmic Membrane'
# etiketi alıyordu → Gram-negatif mantığıyla yazılmış liste onları eliyordu.
# Bkz. benchmark/MISS_ATTRIBUTION.md §2.
#
# KAPSAM: 5 Gram(+) tür. Gram(−) 5 tür (coxiella, brucella*, salmonella, leptospira)
# etkilenmiyor — izinli listeleri değişmedi, bit-bit aynı sonucu verirler.
#
# BEKLENTİ: CapB (banthracis) + ActA (listeria) geri gelir → genel recall 37/48 → 39/48.
# Gram(+) türlerde aday sayısı artacağı için fold-zenginleşme bir miktar DÜŞER
# (ödünleşmenin ölçülmesi gereken tarafı; makalede açıkça çerçevelendi).
#
# Rakipler (Vaxign-ML, NERVE 2.0) etkilenmez → yalnız VaxForge + compare koşulur.
# SIRALI: her koşu ~12 çekirdek kullanır, oversubscribe olmasın. Küçük→büyük proteom.
set -uo pipefail
ROOT="/home/msi-nb/vaxforge"; BENCH="$ROOT/benchmark"
source "$ROOT/.venv/bin/activate" 2>/dev/null
LOG="$BENCH/results/rerun_grampos.log"
echo "=== GRAM(+) RE-RUN BAŞLADI $(date) ===" | tee "$LOG"
python3 "$BENCH/checkpoint.py" "Gram(+) re-run başladı" "5 tür sıralı: ssuis listeria saureus mbovis banthracis; allowed_gram_positive aktif" >/dev/null 2>&1

for BUG in ssuis listeria saureus mbovis banthracis; do
  echo "--- [$BUG] VaxForge re-run (gram-aware loc) $(date) ---" | tee -a "$LOG"
  SKIP_PREP=1 BUG="$BUG" python3 -u "$BENCH/run_bacterium.py" \
      > "$BENCH/results/${BUG}_vaxforge_grampos.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "!! [$BUG] VaxForge rc=$rc — bu tür durdu, sıradakine geçiliyor" | tee -a "$LOG"
  fi
  python3 "$BENCH/compare_bacterium.py" "$BUG" \
      > "$BENCH/results/${BUG}_compare_grampos.log" 2>&1
  REC="$(grep -m1 'Recall@antijen' "$BENCH/results/${BUG}_vaxforge_grampos.log" 2>/dev/null)"
  echo "--- [$BUG] BİTTİ $(date): ${REC:-recall okunamadı} ---" | tee -a "$LOG"
  python3 "$BENCH/checkpoint.py" "Gram(+) re-run: $BUG bitti" "${REC:-recall okunamadı} (rc=$rc)" >/dev/null 2>&1
done

echo "=== GRAM(+) RE-RUN BİTTİ $(date) ===" | tee -a "$LOG"
python3 "$BENCH/checkpoint.py" "Gram(+) re-run TAMAM" "5 tür koşuldu; CapB/ActA durumu için compare loglarına bak" >/dev/null 2>&1
