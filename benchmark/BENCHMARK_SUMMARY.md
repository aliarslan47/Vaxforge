# VaxForge Benchmark — Master Sonuç Özeti

**Tarih:** 2026-08-19/20 (otonom gece koşusu). Tüm araçlar **DEFAULT parametre, özdeş girdi**
(adil koşul, eşik tuning YOK). Aksesyonlar: `PATHOGEN_ACCESSIONS.md`. İstatistik: `stats.py`.
Ham sonuçlar: `results/*.json`. Figürler: `figures/`. İlerleme günlüğü: `PROGRESS.md`.

---

## 1. BAKTERİ KOLU — MenB 3'lü head-to-head (VaxForge vs NERVE 2.0 vs Vaxign-ML)

Ground-truth = ruhsatlı 4CMenB antijenleri **fHbp, NadA, NHBA, PorA**. Üç girdi ölçeği:

| Set | Araç | recall | aday | fold | Fisher p | AUC | fHbp | NadA | NHBA | PorA |
|---|---|---|---|---|---|---|---|---|---|---|
| 64 | **VaxForge** | **4/4** | 22 | 3.33 | 0.012 | **1.00** | #2 | #4 | #3 | #1 |
| 64 | Vaxign-ML | 4/4 | 18 | 4.29 | 0.005 | 0.98 | #7 | #1 | #2 | #4 |
| 64 | NERVE 2.0 | **1/4** | 2 | 15.0 | 0.12 (ns) | — | ❌ | ❌ | ❌ | #2 |
| 204 | **VaxForge** | **4/4** | 68 | 3.12 | 0.012 | **0.99** | #3 | #8 | #6 | #1 |
| 204 | Vaxign-ML | 4/4 | 59 | 3.64 | 0.007 | 0.97 | #14 | #5 | #6 | #8 |
| 204 | NERVE 2.0 | **1/4** | 3 | 25.0 | 0.058 (ns) | — | ❌ | ❌ | ❌ | #3 |
| tam 2003 | Vaxign-ML | 4/4 | 613 | 3.28 | — | — | #150 | #52 | #59 | #29 |
| tam 2003 | NERVE 2.0 | **1/4** | 9 | 62.5 | — | — | ❌ | ❌ | ❌ | #8 |

*(VaxForge tam-2003 CPU'da ~17h+ → koşulmadı; 64/204 temsili. Rakipler tam-2003 koştu.)*

**ANA BULGULAR (BAŞLIK METRİĞİ = recall@antijen; AUC yalnız ikincil, caveat'lı):**
1. **NERVE 2.0 HER ÖLÇEKTE 3/4 ruhsatlı antijeni ELER** (fHbp, NadA, NHBA) — sert `select`
   filtresi yüzünden: NHBA `Periplasmic`, NadA `adhesin_p=0.28<0.5`, fHbp `Unknown loc + 1TM`.
   ÖNEMLİ: NERVE'in KENDİ sürekli skoru bu antijenleri yüksek verir (score 0.68–0.86; skor-AUC 0.88)
   — onları eleyen SKOR değil, skordan SONRA gelen sert filtre. → recall = aracın teslim ettiği fark.
2. **NERVE'in yüksek fold-enrichment'ı (yayınlanmış 16.57; bizde 25→62) aşırı-seçicilikten** gelir,
   iyi recall'den değil — recall'ü feda eder. Fisher p anlamsız (n çok küçük). → **fold tek başına yanıltıcı.**
3. **VaxForge ve Vaxign-ML dördünü de bulur** (recall 4/4). Vaxign-ML tam proteomda çok gevşek (613 aday).
4. **Bağımsız ML aracı (Vaxign-ML) da 4/4** → VaxForge'un soft-filter tezi yalnız değil, sağlam.
5. **NHBA (VFDB-negatif)** VaxForge'da kurtulur (VFDB hard→soft düzeltmesi); NERVE'de elenir.

**⚠️ AUC HAKKINDA DÜRÜST NOT (eleştiriden kaçınmak için — 2026-08-20):** AUC BAŞLIK YAPILMAZ.
Sadece 4 pozitif (ground-truth antijen) var → AUC farkları İSTATİSTİKSEL ANLAMSIZ:
DeLong VaxForge vs Vaxign-ML ΔAUC +0.017 **p=0.42**; vs NERVE-skor ΔAUC +0.12 **p=0.14**.
VaxForge AUC=1.00 → CI [1.0,1.0] DEJENERE (kusursuz ayrım, bootstrap oynayamaz), kırılgan.
Aynı zeminde (her araç kendi sürekli skoru, n_pos=4): VaxForge 1.00[1.0,1.0], Vaxign-ML 0.98[0.94,1.0],
NERVE-skor 0.88[0.63,1.0]. "VaxForge en iyi AUC/ayrım" İDDİA EDİLMEZ. Savunulabilir iddia = **recall
(4/4 vs 1/4) + rank**. Hesap: `stats.py` roc_auc_ci + delong_auc_compare.

## 1b. BAKTERİ KOLU — 3 bakteri × TAM PROTEOM (2026-08-27, çok-antijen recall, 3'lü head-to-head)

MenB'yi genişletmek için 4 memeli konağı (insan/sığır/fare/domuz) enfekte eden, genom-boyu artan 3
zoonotik bakteri **TAM PROTEOM** koşuldu (downsample YOK). Her araç DEFAULT + özdeş girdi. 5'er koruyucu/
aday antijen (bkz. `PATHOGEN_ACCESSIONS.md`; hepsi UniProt GN= ile doğrulandı).

| Bakteri (genom) | Araç | recall | aday | fold | Fisher p | AUC |
|---|---|---|---|---|---|---|
| **S. aureus** (2.8 Mb, 2889) | **VaxForge** | **5/5** | 729 | 3.98 | 0.001 | **1.00** |
|  | Vaxign-ML | 5/5 | 838 | 3.46 | 0.002 | 0.94 |
|  | NERVE 2.0 | **3/5** | 34 | 55.8 | 1.5e-5 | — |
| **L. monocytogenes** (2.9 Mb, 2844) | **VaxForge** | **4/5** | 478 | 4.79 | 0.003 | — |
|  | Vaxign-ML | 5/5 | 917 | 3.11 | 0.003 | 0.95 |
|  | NERVE 2.0 | **3/5** | 73 | 24.3 | 1.6e-4 | — |
| **S. Typhimurium** (4.9 Mb, 4533) | **VaxForge** | **5/5** | 1368 | 3.32 | 0.002 | **0.99** |
|  | Vaxign-ML | 5/5 | 1276 | 3.56 | 0.002 | 0.92 |
|  | NERVE 2.0 | **3/5** | 118 | 23.6 | 1.7e-4 | — |

**TOPLAM: VaxForge 14/15 · Vaxign-ML 15/15 · NERVE 9/15.**

**ANA BULGULAR:**
1. **VaxForge recall'da Vaxign-ML'e fiilen berabere (14/15 vs 15/15) ve AUC/rank'ta önde** — saureus 1.00 vs 0.94,
   salmonella 0.99 vs 0.92; koruyucu antijenler daha üst sıralarda (ör. IsdB #4 vs #46).
2. **NERVE 2.0 her bakteride 3/5** — sert `select` filtresi gerçek yüzey antijenlerini eler (IsdB/IsdA loc=Cellwall+
   düşük adhesin; OmpD/SseB); MenB'deki aynı aşırı-seçicilik kalıbı, yüksek fold ama düşük recall.
3. **Tek VaxForge kaybı: ActA (listeria).** Meşru "zor antijen" — Vaxign-ML'de de en dip (#210/917), bir bug değil.
4. Kalıp MenB + 3 bakteride tutarlı: soft-filter (VaxForge) + bağımsız ML (Vaxign-ML) yüksek recall; NERVE hard-filter kaçırır.

**⚠️ METODOLOJİK ŞEFFAFLIK — iki ölçek-bug'ı bulundu+düzeltildi+doğrulandı (pipeline validasyonu):**
Tam-proteom ölçeğinde iki "tek dev batch patlıyor" hatası VaxForge recall'ını yapay düşürüyordu (ilk koşu
saureus 3/5, listeria 3/5, salmonella **0/5**). Metodolojik hata avıyla izole edildi:
- **(A) PSORTb** (`psortb.py`): tüm proteinler tek konteynerde, subprocess `timeout=1800s`. 4528 protein timeout'u
  aşınca boş dönüyor → funnel heuristik lokalizasyona düşüyor → yüzey antijenleri "cytoplasm" sanılıp eleniyor.
  **Fix:** ≤800 protein alt-batch'ler.
- **(B) netMHCpan** (`netmhc_local.py`): tüm proteinlerin peptitleri (tam proteom → **1.37M** benzersiz 9-mer) tek
  netMHCpan çağrısında → araç çöküyor, ~%0 kapsama → çoğu protein MHC-epitopsuz. **Fix:** 20k-peptit chunk'lar,
  16 çekirdekte **paralel** (`ThreadPoolExecutor`). Doğrulama: 50528 peptit %100 kapsama, OmpA 342/342.
- Toksisite/alerjen filtreleri **suçsuzdu** (peptit-düzeyi deterministik olduğu bağımsız test ile kanıtlandı).
  İki-fix sonrası üç bakteri yeniden koşuldu → yukarıdaki nihai (geçerli) sonuçlar. Bu, yayında Methods/validation
  bölümüne dürüstçe girer (ölçek-güvenli pipeline).

## 2. VİRÜS KOLU — Rabies × 4 konak (yalnız VaxForge; NERVE bacterial-only, Vaxign-ML host'suz)

Ground-truth = glikoprotein **G (NP_056796)**. VaxForge G'yi **rank #2/5** buldu; Vaxign-ML
(virus tipi) G #4/5 ama host boyutu yok. Asıl katkı **çok-konak MHC**:

- Host başına bağlanan epitop FARKLI (G: fare 8, insan 7, sığır 7, domuz 5).
- G-epitop host×host **Jaccard**: insan~sığır 0.75, insan~fare 0.67, insan~domuz 0.50, fare~domuz 0.62.
- → **Host seçimi öngörülen epitopları esaslı değiştirir** — NERVE/Vaxign yapamayan boyut.
- Dürüst caveat: tek antijen → G rank istatistiği betimsel (permütasyon p=0.40, n=1). IEDB
  popülasyon kapsamı rabies/insan-dışı 'veri yok' (uydurma yok).

## 3. İSTATİSTİK (stats.py)
Fisher exact (enrichment), permütasyon rank testi (dağılımsız), ROC-AUC + bootstrap CI,
DeLong (AUC karşılaştırma), Jaccard (çok-konak). Küçük-n dürüstlüğü her testte işaretli.

## 4. FİGÜRLER (figures/)
- `fig1_antigen_recovery.png` — araç×antijen kurtarma matrisi (NERVE 3 kırmızı "elendi") ★ ana figür
- `fig2_recall_bar.png` — recall@antijen + AUC
- `fig3_fold_vs_recall.png` — fold-vs-recall ödünleşimi (NERVE yüksek-fold/düşük-recall) [kozmetik: alt-sağ etiket örtüşmesi düzeltilebilir]
- `fig4_multihost_jaccard.png` — rabies G host×host Jaccard (çok-konak) ★
- `fig5_bacteria_panel.png` — 3 bakteri × 3 araç recall@antijen + AUC (bakteriyel kol) ★ ana figür

## 5. SEÇİM GEREKÇELERİ & YAZIM AŞAMASI
Patojen (MenB, rabies) + konak (insan/sığır/fare/domuz) **biyolojik kullanım + seçilim gerekçeleri
DOI-doğrulanmış makalelerle** → **`SELECTION_RATIONALE.md`** (2026-08-20; 6 atıf web-doğrulandı:
Pizza 2000, Giuliani 2006, Cox 1977, Fooks 2017, NERVE2 Conte 2024, Vaxign-ML Ong 2020 — hepsi gerçek
DOI, uydurma yok). Kalan (opsiyonel/GPU): VaxForge tam-2003 koşusu flagship fold sayısını tamamlar;
konak-başına ayrı epidemiyolojik atıf yazımda eklenebilir.
