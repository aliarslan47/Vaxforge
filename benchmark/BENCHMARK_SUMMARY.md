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

## 5. YAZIM AŞAMASINDA EKLENECEK (kullanıcı kuralı)
Patojen (MenB, rabies) + konak (insan/sığır/fare/domuz) **biyolojik kullanım + seçilim gerekçeleri
DOI-doğrulanmış makalelerle** desteklenecek. Uydurma DOI YASAK. Ayrıca VaxForge tam-2003 koşusu
(GPU gelince) flagship fold sayısını tamamlar.
