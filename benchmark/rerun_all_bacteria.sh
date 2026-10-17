#!/bin/bash
# FAZ 2 — 10 bakteri türünün TEK KOD SÜRÜMÜNDE yeniden koşumu.
#
# GEREKÇE: depodaki bakteri kıyas verisi tümüyle kalibrasyon-ÖNCESİ. Üç kapı
# kalibrasyonu arada girdi:
#   1) toksisite kapısı uzunluk sınıfına ayrıldı (2026-10-10)
#   2) otoimmünite epitop granülaritesine taşındı, protein düzeyi yumuşadı (2026-10-10)
#   3) TM kapısı sinyal-duyarlı sayıma geçti (2026-10-18)
# (1) ve (2) tüm türleri, (3) de tüm bakterileri etkiliyor (virüste TM yumuşak →
# virüs kolu ETKİLENMİYOR, yeniden koşulmaz). Yarım kalan Gram(+) koşumu bilerek
# geri alınmıştı; tablonun iç tutarlı olması için 10 türün hepsi aynı kodla koşar.
#
# RAKİPLER YENİDEN KOŞULMAZ: Vaxign-ML ve NERVE 2.0 dış araçlar, bizim kapılarımızdan
# etkilenmiyorlar; mevcut sayıları geçerli kalır (benchmark/results/vaxignml_*).
#
# SIRALI, küçük→büyük proteom. Her koşu ~12 çekirdek ister; makinede 16 çekirdek
# ve 12 GB RAM var → paralel koşum RAM'de boğulur, oversubscribe edilmez.
# Beklenen süre: tür başına ~2-3 saat, toplam ~25-30 saat.
#
# Kopma dayanıklılığı: her tür bitiminde checkpoint.py'ye kayıt + ayrı log.
# Bir tür çökerse sıradakine geçilir (rc loga yazılır), koşum durmaz.
set -uo pipefail
ROOT="/home/msi-nb/vaxforge"; BENCH="$ROOT/benchmark"
source "$ROOT/.venv/bin/activate" 2>/dev/null
LOG="$BENCH/results/rerun_all.log"
STAMP="$(date '+%Y%m%d_%H%M')"

BUGS="coxiella ssuis saureus listeria brucella leptospira mbovis salmonella banthracis"

{
echo "=== 10 TÜR YENİDEN KOŞUM BAŞLADI $(date) ==="
echo "sıra: $BUGS + menb (ayrı betik)"

# --- kalibrasyon-öncesi sonuçları arşivle (SİLME, yeniden adlandırma) ---
for BUG in $BUGS menb; do
  d="$BENCH/results/$BUG"
  if [ -d "$d" ] && [ ! -d "${d}_pre_calib_${STAMP}" ]; then
    cp -a "$d" "${d}_pre_calib_${STAMP}" && echo "arşiv: ${BUG} -> ${BUG}_pre_calib_${STAMP}"
  fi
done

python3 "$BENCH/checkpoint.py" "Faz 2 başladı" \
  "10 tür sıralı, tek kod sürümü; TM sinyal-duyarlı + toksisite + otoimmünite kalibre" >/dev/null 2>&1

for BUG in $BUGS; do
  echo "--- [$BUG] VaxForge $(date) ---"
  SKIP_PREP=1 BUG="$BUG" python3 -u "$BENCH/run_bacterium.py" \
      > "$BENCH/results/${BUG}_vaxforge_postcal_all.log" 2>&1
  rc=$?
  [ $rc -ne 0 ] && echo "!! [$BUG] rc=$rc — bu tür durdu, sıradakine geçiliyor"
  python3 "$BENCH/compare_bacterium.py" "$BUG" \
      > "$BENCH/results/${BUG}_compare_postcal_all.log" 2>&1
  REC="$(grep -m1 'Recall@antijen' "$BENCH/results/${BUG}_vaxforge_postcal_all.log" 2>/dev/null)"
  echo "--- [$BUG] BİTTİ $(date): ${REC:-recall okunamadı} (rc=$rc) ---"
  python3 "$BENCH/checkpoint.py" "Faz 2: $BUG bitti" "${REC:-recall okunamadı} (rc=$rc)" >/dev/null 2>&1
done

# --- 10. tür: MenB kendi betiğiyle ---
echo "--- [menb] VaxForge $(date) ---"
python3 -u "$ROOT/scripts/validate_menb.py" > "$BENCH/results/menb_vaxforge_postcal_all.log" 2>&1
rc=$?
[ $rc -ne 0 ] && echo "!! [menb] rc=$rc"
python3 -u "$BENCH/compare_menb.py" > "$BENCH/results/menb_compare_postcal_all.log" 2>&1
echo "--- [menb] BİTTİ $(date) (rc=$rc) ---"
python3 "$BENCH/checkpoint.py" "Faz 2: menb bitti" "rc=$rc" >/dev/null 2>&1

echo "=== 10 TÜR YENİDEN KOŞUM BİTTİ $(date) ==="
python3 "$BENCH/checkpoint.py" "Faz 2 TAMAM" \
  "10 tür tek kod sürümünde koşuldu; stats + figür + makale güncellemesi sırada" >/dev/null 2>&1
} >> "$LOG" 2>&1
