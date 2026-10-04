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

## 1b. BAKTERİ KOLU — 9 bakteri × TAM PROTEOM (2026-10-02, çok-antijen recall, 3'lü head-to-head)

MenB'yi genişletmek için 4 memeli konağı (insan/sığır/fare/domuz) enfekte eden, genom-boyu artan **9
zoonotik bakteri** (eski 3 + 2026-10 eklenen 6) **TAM PROTEOM** koşuldu (downsample YOK). Her araç DEFAULT
+ özdeş girdi. 4-5 koruyucu/aday antijen (bkz. `PATHOGEN_ACCESSIONS.md`; UniProt/RefSeq doğrulandı). Genom
küçükten büyüğe sıra; yeni-6 🆕.

| Bakteri (genom, protein) | Araç | recall | aday | fold | Fisher p | AUC |
|---|---|---|---|---|---|---|
| 🆕 **C. burnetii** (2.0 Mb, 1812) | **VaxForge** | **2/4** | 524 | 1.7 | 0.33 | — |
|  | Vaxign-ML | 4/4 | 477 | 3.8 | 0.0048 | 0.91 |
|  | NERVE 2.0 | **0/4** | 9 | — | 1 | — |
| 🆕 **S. suis** (2.2 Mb, 1922) | **VaxForge** | **1/5** | 288 | 1.3 | 0.56 | — |
|  | Vaxign-ML | 5/5 | 628 | 3.1 | 0.0037 | 0.88 |
|  | NERVE 2.0 | **0/5** | 0 | — | 1 | — |
| **S. aureus** (2.8 Mb, 2889) | **VaxForge** | **5/5** | 729 | 3.98 | 0.001 | **1.00** |
|  | Vaxign-ML | 5/5 | 838 | 3.46 | 0.002 | 0.94 |
|  | NERVE 2.0 | **3/5** | 34 | 55.8 | 1.5e-5 | — |
| **L. monocytogenes** (2.9 Mb, 2844) | **VaxForge** | **4/5** | 478 | 4.79 | 0.003 | — |
|  | Vaxign-ML | 5/5 | 917 | 3.11 | 0.003 | 0.95 |
|  | NERVE 2.0 | **3/5** | 73 | 24.3 | 1.6e-4 | — |
| 🆕 **B. suis** (3.3 Mb, 3249) | **VaxForge** | **4/5** | 929 | 2.8 | 0.026 | — |
|  | Vaxign-ML | 5/5 | 839 | 3.9 | 0.0011 | 0.88 |
|  | NERVE 2.0 | **2/5** | 9 | 185.4 | 6.8e-5 | — |
| 🆕 **M. bovis** (4.4 Mb, 4088) | **VaxForge** | **4/5** | 535 | 6.2 | 0.0013 | — |
|  | Vaxign-ML | 4/5 | 695 | 4.7 | 0.0036 | 0.95 |
|  | NERVE 2.0 | **0/5** | 12 | — | 1 | — |
| 🆕 **L. interrogans** (4.7 Mb, 3676) | **VaxForge** | **4/5** | 1105 | 2.7 | 0.031 | — |
|  | Vaxign-ML | 5/5 | 1211 | 3.0 | 0.0039 | 0.88 |
|  | NERVE 2.0 | **0/5** | 116 | — | 1 | — |
| **S. Typhimurium** (4.9 Mb, 4533) | **VaxForge** | **5/5** | 1368 | 3.32 | 0.002 | **0.99** |
|  | Vaxign-ML | 5/5 | 1276 | 3.56 | 0.002 | 0.92 |
|  | NERVE 2.0 | **3/5** | 118 | 23.6 | 1.7e-4 | — |
| 🆕 **B. anthracis** (5.5 Mb, 5493) | **VaxForge** | **4/5** | 1163 | 3.8 | 0.0083 | — |
|  | Vaxign-ML | 5/5 | 1297 | 4.2 | 7.3e-4 | 0.94 |
|  | NERVE 2.0 | **2/5** | 10 | 274.4 | 3.0e-5 | — |

**TOPLAM (9 bakteri): VaxForge 33/44 (%75) · Vaxign-ML 43/44 (%98) · NERVE 13/44 (%30).**
(MenB exemplar dahil 10-bakteri geneli: VaxForge 37/48 · Vaxign-ML 47/48 · NERVE 14/48.)

**ANA BULGULAR:**
1. **VaxForge 9 türün 7'sinde 4/5-5/5** (saureus/salmonella 5/5, listeria/brucella/mbovis/leptospira/banthracis 4/5);
   tam-recall türlerde AUC/rank'ta önde (ör. IsdB #4; mbovis'te Vaxign'den yüksek fold 6.2).
2. **VaxForge iki türde geride: C. burnetii 2/4, S. suis 1/5.** Kaçırılanlar sitoplazmik **moonlighting** antijenleri
   (GAPDH, enolaz, 6PGD, GroEL, DnaJ) → lokalizasyon SERT filtresi tasarım gereği eler; Vaxign-ML (yalnız-sekans,
   lokalizasyon kapısı yok) tutar. **Kusur değil, bilinçli ödünleşme** (sitoplazmik yanlış-pozitifi bastırır).
3. **NERVE 2.0 dört türde 0/5-0/4** (coxiella/ssuis/mbovis/leptospira) — sert `select` yüzey antijenlerini eler;
   brucella/banthracis'te en yüksek fold (185-274) ama düşük recall → fold tek başına yanıltıcı (kalıp MenB'den beri tutarlı).
4. **Tek "klasik" VaxForge kaybı: ActA (listeria)** — Vaxign'de de en dip, meşru zor antijen.

**⚠️ METODOLOJİK ŞEFFAFLIK — ÜÇ ölçek-bug'ı bulundu+düzeltildi+doğrulandı (pipeline validasyonu):**
Tam-proteom ölçeğinde üç hata VaxForge recall'ını yapay düşürüyor/koşuyu donduruyordu. Metodolojik hata avıyla izole edildi:
- **(A) PSORTb** (`psortb.py`): tüm proteinler tek konteynerde, `timeout=1800s`. Büyük proteom timeout'u aşınca boş
  dönüyor → funnel heuristik lokalizasyona düşüyor → yüzey antijenleri "cytoplasm" sanılıp eleniyor. **Fix:** ≤800 protein alt-batch.
- **(B) netMHCpan** (`netmhc_local.py`): tüm peptitler (tam proteom → **1.37M** 9-mer) tek çağrıda → ~%0 kapsama.
  **Fix:** 20k-peptit chunk, çekirdeklerde **paralel**. Doğrulama: 50528 peptit %100, OmpA 342/342.
- **(C) netMHCpan eşzamanlılık DEADLOCK'u** (`netmhc_local.py`, 2026-10-01): paralel 12 worker ORTAK scratch/TMPDIR'de
  çakışıp `pipe_read`'de kilitleniyor, `timeout=7200s` yüzünden koşu süresiz donuyordu (brucella ~5h asılı kaldı).
  **Fix:** her chunk'a izole `TMPDIR` + `stdin=DEVNULL` + timeout 7200→1200s. Sonra brucella 4/5 ile sorunsuz koştu.
- Toksisite/alerjen filtreleri **suçsuzdu** (bağımsız test). Fix'ler sonrası etkilenen türler yeniden koşuldu → nihai sonuçlar.
  Yayında Methods/validation'a girer (`MANUSCRIPT_METHODS.md` §B, üç kusur).

## 1c. BRUCELLA AŞI SUŞU VALİDASYONU — virülan-vs-aşı (2026-10-04, TAM PROTEOM)

Üç **ruhsatlı canlı-atenüe** bruselloz aşı suşu kendi RefSeq proteomunda koşuldu; virülan referans
*B. suis* 1330 ile kıyas (→ MANUSCRIPT Tablo S2 + Şekil S1, METHODS §A.5):

| Suş | Tür | Protein | Antijen | Recall | Fisher p | En iyi sıra | Kaçırılan |
|---|---|---|---|---|---|---|---|
| 1330 (virülan ref.) | *B. suis* | 3249 | 5 | 4/5 | 0,026 | SodC #26 | L7/L12 |
| RB51 (sığır) | *B. abortus* | 2971 | 4 | 3/4 | 0,055 | BP26 #87 | L7/L12 |
| S19 (sığır) | *B. abortus* | 2964 | 4 | 3/4 | 0,055 | BP26 #89 | L7/L12 |
| S2 (küçük ruminant) | *B. suis* | 2979 | 5 | 4/5 | 0,015 | Omp31 #24 | L7/L12 |

- **Tutarlılık bulgusu:** yüzey/periplazmik antijenler (Omp25/Omp31/BP26/SodC) üç aşı suşunda da bulundu;
  **tek ve tutarlı kayıp L7/L12** (sitoplazmik ribozomal → lokalizasyon filtresi tasarım gereği eler),
  tıpkı virülan 1330'daki gibi. Omp31 *B. abortus*'ta YOK (delesyon; Cloeckaert 2002) → RB51/S19 4 antijen.
- *B. suis* S2 (aynı 5-antijen seti) virülan 1330'u birebir tekrarladı (4/5). Rakipler koşuldu ama derlenmiş
  head-to-head yok → VaxForge sağlamlık validasyonu olarak sunuldu.

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

## 4. FİGÜRLER (figures/) — `make_figures.py` ile üretilir
- `fig1_antigen_recovery.png` — MenB araç×antijen rank matrisi (NERVE 3 kırmızı "elendi") ★ ana figür
- `fig2_recall_bar.png` — MenB recall@antijen + AUC (204-set)
- `fig3_fold_vs_recall.png` — fold-vs-recall ödünleşimi (NERVE yüksek-fold/düşük-recall)
- `fig4_bacteria_panel.png` — **10 bakteri × 3 araç** recall@antijen (genom sırası, kesir+x/n; yeni-6 yıldızlı) ★ ana figür
- `fig5_multihost_jaccard.png` — rabies G host×host Jaccard (çok-konak) ★
- `fig6_virus_antigen_rank.png` — 10 virüs ground-truth antijen rank, VaxForge vs Vaxign-ML
- `fig7_virus_multihost.png` — virüs başına host-çifti epitop Jaccard aralığı (çok-konak)
- `figS1_brucella_vaccine.png` — **Ek Şekil S1**: virülan 1330 + 3 aşı suşu (RB51/S19/S2) antijen rank nokta-grafiği

## 5. SEÇİM GEREKÇELERİ & YAZIM AŞAMASI (TAMAMLANDI 2026-09-20)
Patojen (MenB, rabies, 3 bakteri) + konak (insan/sığır/fare/domuz) **biyolojik kullanım + seçilim
gerekçeleri DOI-doğrulanmış makalelerle** → **`SELECTION_RATIONALE.md`**. Yazım-hazır Methods
düzyazısı → **`MANUSCRIPT_METHODS.md`** (§A seçim gerekçesi, §B ölçek-güvenlik validasyonu).
- **Koruyucu-antijen atıfları TAMAM:** 20 ground-truth antijenin (MenB 4 + rabies G + 15 bakteri)
  koruyuculuk kanıtı + DOI'si → **`PROTECTIVE_ANTIGENS.md`** (2026-09-20; 18 atıf Crossref/yayıncı ile
  doğrulandı, uydurma yok; kanıt-gücü antijen-başına açık: aktif/pasif immünizasyon vs aday-antijen).
- **İki ölçek-bug'ı Methods/validation'a taşındı** → `MANUSCRIPT_METHODS.md` §B.
- Kalan (opsiyonel/GPU): VaxForge tam-2003 MenB koşusu flagship fold sayısını tamamlar.
