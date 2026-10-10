# Filtre kalibrasyonu — kapıların bağımsız kontrol setlerindeki davranışı

Üretim: `benchmark/calibrate_filters.py` · veri: `benchmark/results/filter_calibration.json` · seed: 0

Her hücre, o kontrol setinde **işaretlenen** (yani filtrenin eleyeceği) peptit/protein yüzdesidir. `tox_protein` dışındaki her satırda yüksek değer = **yanlış-pozitif**: filtrenin elememesi gereken diziyi eliyor.

Set boyutları: `self9`=1000, `self15`=1000, `iedb9`=1000, `iedb15`=1000, `tox_protein`=300, `nontox_protein`=300

## Kapılar

| kapı | insan self 9-mer | insan self 15-mer | IEDB doğrulanmış epitop (9-mer) | IEDB doğrulanmış epitop (13-17-mer) | toksin proteini (tam uzunluk) — **duyarlılık** | non-toksin proteini (tam uzunluk) |
|---|---|---|---|---|---|---|
| `toxicity_ToxinPred2_hybrid_t0.6` | %75.2 | %66.3 | %63.1 | %57.7 | %94.3 | %4.7 |
| `toxicity_ToxinPred3_ET_t0.38` | %31.7 | %27.0 | %23.8 | %23.0 | %30.0 | %3.0 |
| `toxicity_ToxinPred3_hybrid_t0.38` | %22.2 | %14.0 | %14.6 | %14.0 | %9.0 | %1.0 |
| `allergenicity_FAO_WHO_6mer` | %3.8 | %7.9 | %2.0 | %6.1 | %65.3 | %83.0 |
| `mimicry_human_exact_kmer(soft)` | %100.0 | %100.0 | %1.3 | %1.7 | %51.3 | %66.3 |

## Eşik taraması

### ToxinPred2_hybrid

| eşik | insan self 9-mer | insan self 15-mer | IEDB doğrulanmış epitop (9-mer) | IEDB doğrulanmış epitop (13-17-mer) | toksin proteini (tam uzunluk) — **duyarlılık** | non-toksin proteini (tam uzunluk) |
|---|---|---|---|---|---|---|
| 0.05 | %100.0 | %100.0 | %100.0 | %99.9 | %100.0 | %48.0 |
| 0.1 | %100.0 | %100.0 | %100.0 | %99.7 | %100.0 | %43.7 |
| 0.15 | %100.0 | %100.0 | %100.0 | %99.4 | %100.0 | %38.3 |
| 0.2 | %100.0 | %99.9 | %100.0 | %99.2 | %99.0 | %31.0 |
| 0.25 | %100.0 | %99.9 | %100.0 | %99.0 | %97.0 | %23.7 |
| 0.3 | %100.0 | %99.9 | %100.0 | %98.9 | %95.0 | %17.0 |
| 0.35 | %100.0 | %99.9 | %100.0 | %98.8 | %94.7 | %11.7 |
| 0.4 | %99.9 | %99.5 | %99.9 | %98.3 | %94.3 | %9.3 |
| 0.45 | %99.5 | %97.9 | %98.7 | %95.3 | %94.3 | %7.0 |
| 0.5 | %97.8 | %92.6 | %94.7 | %89.7 | %94.3 | %6.0 |
| 0.55 | %90.2 | %81.9 | %84.8 | %78.8 | %94.3 | %5.0 |
| 0.6 | %75.9 | %66.9 | %64.0 | %58.0 | %94.3 | %4.7 |
| 0.65 | %51.1 | %48.3 | %36.9 | %40.0 | %94.0 | %3.7 |
| 0.7 | %29.8 | %32.1 | %20.4 | %28.5 | %93.7 | %2.3 |
| 0.75 | %18.1 | %20.9 | %10.9 | %17.8 | %92.3 | %1.7 |
| 0.8 | %14.3 | %13.2 | %8.5 | %11.4 | %92.3 | %1.3 |
| 0.85 | %8.6 | %6.3 | %5.0 | %5.5 | %92.3 | %1.0 |
| 0.9 | %3.5 | %2.5 | %1.3 | %1.5 | %92.3 | %1.0 |
| 0.95 | %0.3 | %0.3 | %0.0 | %0.1 | %92.3 | %0.7 |

### ToxinPred3_ET

| eşik | insan self 9-mer | insan self 15-mer | IEDB doğrulanmış epitop (9-mer) | IEDB doğrulanmış epitop (13-17-mer) | toksin proteini (tam uzunluk) — **duyarlılık** | non-toksin proteini (tam uzunluk) |
|---|---|---|---|---|---|---|
| 0.05 | %99.9 | %99.6 | %100.0 | %100.0 | %100.0 | %100.0 |
| 0.1 | %97.7 | %94.4 | %97.5 | %93.3 | %98.7 | %97.7 |
| 0.15 | %89.4 | %81.0 | %87.4 | %76.6 | %89.7 | %88.3 |
| 0.2 | %75.4 | %63.8 | %71.4 | %59.1 | %68.3 | %64.7 |
| 0.25 | %59.3 | %49.4 | %54.4 | %43.4 | %46.3 | %32.3 |
| 0.3 | %45.5 | %40.4 | %39.4 | %33.6 | %35.0 | %10.7 |
| 0.35 | %36.6 | %31.4 | %28.6 | %26.7 | %31.0 | %4.7 |
| 0.4 | %29.0 | %23.9 | %21.5 | %20.0 | %29.3 | %2.3 |
| 0.45 | %22.8 | %17.1 | %16.3 | %13.9 | %27.7 | %1.3 |
| 0.5 | %17.5 | %12.3 | %11.7 | %8.8 | %24.0 | %0.3 |
| 0.55 | %13.8 | %8.9 | %7.5 | %6.9 | %19.0 | %0.3 |
| 0.6 | %10.1 | %6.4 | %4.5 | %4.9 | %13.0 | %0.3 |
| 0.65 | %7.1 | %5.0 | %3.0 | %3.5 | %6.3 | %0.3 |
| 0.7 | %4.9 | %2.6 | %1.8 | %2.2 | %3.0 | %0.3 |
| 0.75 | %4.0 | %1.9 | %0.7 | %1.5 | %1.7 | %0.0 |
| 0.8 | %2.7 | %0.8 | %0.6 | %0.9 | %1.0 | %0.0 |
| 0.85 | %1.8 | %0.7 | %0.2 | %0.6 | %0.3 | %0.0 |
| 0.9 | %0.7 | %0.4 | %0.1 | %0.5 | %0.0 | %0.0 |
| 0.95 | %0.0 | %0.1 | %0.1 | %0.0 | %0.0 | %0.0 |

## Yorum ve karar

1. **Toksisite, epitop düzeyinde sert filtre olamaz.** ToxinPred2 (Model 2 hibrit, eşik 0.6) insan self 9-mer'lerinin ve deneysel doğrulanmış IEDB epitoplarının çoğunluğunu "toksik" işaretliyor; influenza M1(58-66) GILGFVFTL dahil. Araç tam uzunluk proteinde ise yüksek duyarlılık ve düşük yanlış-pozitifle çalışıyor. Sorun aracın kalitesi değil, **uzunluk sınıfı**: ToxinPred2 'toxicity of proteins' için eğitilmiştir.

2. **Epitop uzunluğunda altın standart yok.** Bir toksinden kesilen 9-mer'in kendisi toksik değildir (toksisite katlanmış proteinin işlevidir), dolayısıyla bu uzunlukta duyarlılık ölçülemez. Ölçülebilen tek şey yanlış-pozitif oranıdır; bu da mimari kararı dayatıyor: epitopta **yumuşak işaretleme**, proteinde/MEV'de **sert eleme**.

3. **ToxinPred3 peptit düzeyi için doğru araç** (Rathore ve ark. 2024, Comput Biol Med 179:108926) ama protein düzeyinde ToxinPred2'nin yerini almaz — tablodaki `tox_protein` satırı bunu gösteriyor.

4. **Alerjenite ölçütü uzunluğa göre seçilmeli.** FAO/WHO 6-mer kuralı 9-mer'de özgül; tam proteinde non-alerjen proteinleri de işaretliyor. MEV konstruktu bir protein olduğu için orada Codex/FAO-WHO birincil ölçütü (80 aa pencerede ≥%35 kimlik, diamond) kullanılır.

5. **Mimikri kapısı sağlam.** İnsan self k-mer'lerinde ~%100, IEDB epitoplarında ~%0 işaretleme — tasarlandığı gibi ayırıyor ve zaten yumuşak.


> Kural: bir sert kapı, bu setlerde kalibre edilmeden üretime girmez.
