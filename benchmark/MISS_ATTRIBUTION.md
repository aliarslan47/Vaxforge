# Kaçan antijenlerin mekanistik atfı — 10 bakteri, 48 antijen

Son güncelleme: 2026-10-04 · Veri: `results/miss_attribution.json` · Analiz: diamond blastp
(`tools/db/human.dmnd`, UniProt Swiss-Prot insan proteomu, e-value 1e-4) + PSORTb + TMHMM.

## Soru

10-bakteri head-to-head'de VaxForge 37/48 (%77), Vaxign-ML 47/48 (%98), NERVE2 14/48 (%29).
VaxForge'un 11 kaybı önceki iki vakadaki gibi bir ölçek-bug'ı mı, yoksa yöntemsel mi?

## Yanıt: 11 kaybın tamamı atfedildi, ölçek-bug'ı YOK

| Sebep | n | Antijenler |
|---|---|---|
| Konak-homoloji (otoimmünite) sert filtresi | 6 | gnd, GroEL/HtpB, enolaz, GAPDH, Mip, DnaJ |
| Lokalizasyon filtresi — Gram(+) `membrane` açığı | 2 | CapB (*B. anthracis*), ActA (*L. monocytogenes*) |
| PSORTb yanlış atama (üst-araç limiti) | 1 | LigA (*Leptospira*) |
| Epitop/bağlanma aşaması (tahmin düzeyi) | 2 | ESAT-6 (*M. bovis*), L7/L12 (*Brucella*) |

## 1. Konak-homoloji filtresi — tasarım gereği, 6/11

`config/thresholds.yaml` → `human_homology.min_identity = 35`, `hard_filter: true`.

Kontenjans (insan %kimlik eşiği %35):

|  | insan ≥%35 | insan <%35 | toplam |
|---|---|---|---|
| **kaçan** | 6 | 5 | 11 |
| **yakalanan** | 0 | 37 | 37 |

Fisher exact p = 3.76e-05, OR = ∞.

İlişki tek yönde **mutlak**: yakalanan 37 antijenin hiçbirinde insan hit'i yok (0/37, çoğunda
e-value bile üretmiyor). Kaçan 6'sı %37.0–52.3 aralığında. Bu 6'sının hepsi korunmuş housekeeping
proteini (şaperon / glikolitik enzim / izomeraz) — insan ortoloğu bulunması kaçınılmaz sınıf.

**Yorum.** VaxForge bu antijenleri beceremediği için değil, otoimmünite gerekçesiyle **bilerek**
atıyor. Vaxign-ML konak-homoloji elemesi yapmıyor (bilinen antijenler üzerinde eğitilmiş
sınıflandırıcı), dolayısıyla onları üst sıraya koyabiliyor. %77 vs %98 farkının yarısı
başarısızlık değil, ölçülebilir bir **güvenlik–duyarlılık ödünleşmesi**.

**Dürüstlük kaydı.** (a) %35 bir parametre tercihidir, evrensel bir eşik değil. (b) Homoloji ≠
patoloji: bu proteinlerin gerçekten otoimmünite yapacağı ayrı ve kanıtlanmamış bir iddiadır;
doğru ifade "ihtiyatlı filtre uygulanıyor"dur. (c) HSP60/enolaz otoimmün çapraz-reaktivite
literatürü makaleye girecekse DOI ile doğrulanmalı — bu analizde atıf doğrulaması YAPILMADI.

## 2. Lokalizasyon — Gram-pozitif `membrane` açığı, 2/11 (DÜZELTİLEBİLİR)

`config/thresholds.yaml` izinli lokalizasyonlar:

```yaml
bacteria: [extracellular, outer_membrane, cell_wall, periplasm]
parasite: [extracellular, plasma_membrane, secreted, gpi_anchored]
```

Parazit listesinde `plasma_membrane` var, bakteri listesinde yok. Gram-pozitif bakteride dış zar
yoktur — sitoplazmik zar hücrenin dış sınırıdır ve yüzeye çıkan proteinler PSORTb'de `membrane`
etiketi alır. CapB ve ActA tam bu yüzden eleniyor (ikisi de Gram(+), ikisi de `membrane` skor 10.0,
TM=1). Bakteri listesi Gram-negatif mantığıyla yazılmış; Gram-pozitife uygulandığında meşru yüzey
antijenlerini kesiyor.

Bu, 1. maddeden farklı olarak savunulabilir bir tasarım tercihi değil, **tutarsızlık**.

**DÜZELTİLDİ (2026-10-10).** `config/thresholds.yaml`'a `localization.allowed_gram_positive`
eklendi (`[extracellular, cell_wall, membrane]`); `funnel.run` artık `profile == "bacteria" and
gram == "positive"` olduğunda bu listeyi seçiyor (`evaluate.candidate_rows` rapor tarafında da
gram-duyarlı). Gömülü çok-sarmallı proteinleri `transmembrane.max_helices` ayıklamaya devam
ediyor, yani özgüllük güvencesi korunuyor. Regresyon testi:
`tests/test_modules.py::test_localization_allowed_is_gram_aware`.

**Yeniden koşum kapsamı — düzeltme:** burada başlangıçta "10 bakterinin tamamı" yazılmıştı, bu
yanlıştı. İzinli liste yalnız Gram-pozitif dalda değişti, dolayısıyla etkilenen **5 tür**:
`saureus`, `listeria`, `ssuis`, `mbovis`, `banthracis`. Beş Gram-negatif tür (`coxiella`,
`brucella`, `salmonella`, `leptospira` + Brucella aşı suşları) bit-bit aynı sonucu verir.
Beklenen etki: CapB ve ActA geri gelirse recall 37/48 → 39/48; Gram(+) türlerde aday sayısı
artacağı için fold-zenginleşme bir miktar düşecektir (ödünleşmenin ölçülmesi gereken tarafı).

## 3. PSORTb yanlış ataması, 1/11

LigA (*Leptospira*) → PSORTb `cytoplasm` skor 9.97. Literatürde bilinen yüzey lipoproteini.
VaxForge'un değil, üst-aracın hatası; limitasyon olarak raporlanır.

## 4. Epitop/bağlanma aşaması, 2/11 — **YENİDEN SINIFLANDIRILDI (2026-10-10): KALİBRASYON HATASI**

**Bu bölümün ilk hâli yanlıştı.** "Tahmin düzeyinde meşru ıska" denen şey ölçüldü ve
mühendislik/kalibrasyon hatası olduğu bulundu. Ölçüm: `benchmark/calibrate_filters.py`
(sonuç: `results/filter_calibration.json`).

Toksisite kapısı ToxinPred2 Model 2 (hibrit) + eşik 0.6 ile çalışıyordu. Bağımsız
kontrol setlerinde (n=1000):

| kontrol seti | "toksik" işaretlenen |
|---|---|
| insan Swiss-Prot 9-mer (self, toksik değil) | **%75.2** |
| insan Swiss-Prot 15-mer | %66.3 |
| IEDB'de deneysel **doğrulanmış** epitoplar (9-mer) | **%63.1** |
| IEDB doğrulanmış epitoplar (13-17-mer) | %57.7 |
| tam uzunluk toksin proteinleri (duyarlılık) | %94.3 |
| tam uzunluk non-toksin proteinleri | %4.7 |

Influenza M1(58-66) **GILGFVFTL** — literatürün en çok doğrulanmış zararsız insan CD8
epitobu — `ml_score 0.663` ile "toksik" işaretlenip eleniyordu.

Yani kapı epitop uzunluğunda bilgi taşımıyor; rastgele bir 4'te-3 kırımıydı. Az sayıda
epitop üreten bir proteinin tümüyle kaybolma olasılığı bu yüzden yüksek: **ESAT-6 ve
L7/L12 buradan kayboldu.** Toksisite hipotezinin "çürütülmesi" (bkz. aşağıdaki eski
paragraf) de bu yüzden yanıltıcıydı: ölçülen oran proteinler arasında ayırt edici
değildi, çünkü skor dağılımı eşiğin tam üstüne yığılıyor.

**Kök neden araç değil, uzunluk sınıfı.** ToxinPred2'nin başlığı "toxicity of
*proteins*"tir; AAC özellikleri 9-mer'de anlamsızdır. Aynı araç tam proteinde kusursuz
(%94.3 duyarlılık / %4.7 yanlış-pozitif) — bu yüzden MEV konstruktu düzeyinde korunuyor.

**Düzeltme (uygulandı).** Epitop düzeyinde ToxinPred3 (peptit için eğitilmiş, Rathore
2024, doi:10.1016/j.compbiomed.2024.108926) **yumuşak** olarak kullanılıyor: işaretler,
adaylık puanını ×0.9 yapar, eleme yalnız kalibre edilmiş yüksek-güven eşiğinin üstünde.
Sert toksisite elemesi protein/MEV düzeyinde ToxinPred2 hibrit ile sürüyor.
Regresyon testi: `tests/test_modules.py::test_toxicity_gate_is_soft_and_length_appropriate`.

**Yan bulgu — alerjenite.** Aynı tarama FAO/WHO 6-mer kuralının 9-mer'de özgül
(yanlış-pozitif ~%2-3) ama tam proteinde kullanılamaz olduğunu gösterdi (rastgele
non-alerjen proteinlerin neredeyse tamamı işaretlendi). MEV konstruktu bir protein
olduğu için orada Codex/FAO-WHO'nun birincil ölçütüne (80 aa pencerede ≥%35 kimlik,
diamond) geçildi; peptit düzeyinde 6-mer kuralı kalıyor.

### Eski (2026-10-04) metin — kayıt için

ESAT-6 ve L7/L12 üç sert filtreyi de geçti (funnel'da hayatta), 0 aday peptit üretti.

Toksisite hipotezi **test edildi ve çürütüldü**: ESAT-6'nın 9-mer toksik işaretleme oranı %60 —
dört proteinin en düşüğü — buna rağmen 0 peptit; kontrol CFP-10 %82 ile daha yüksekken 2 peptit
üretti. Uzunluk da tek başına açıklamıyor (CFP-10 100 aa, ESAT-6 95 aa). Dolayısıyla kayıp
MHC-bağlanma/skorlama eşiğinde gerçekleşiyor: sistematik kusur değil, tahmin düzeyinde meşru ıska.
L7/L12'de %90 toksik işaretleme katkıda bulunuyor olabilir, ancak tek başına kanıtlanamadı.

## Sonuç — **2026-10-10'da düzeltildi**

İlk hâlinde burada "**mühendislik hatası yok**" yazıyordu. Bu **yanlıştı**: 4. maddedeki
2 kayıp, ölçüldüğünde kalibre edilmemiş bir kapıdan (tam-protein toksisite modelinin
9-mer'lere sert filtre olarak uygulanması) çıktı. Ayrıca 2. maddedeki Gram(+) açığının
düzeltilmesi ActA'yı geri getirmedi; ActA huniyi geçmeye başladı ama aynı toksisite
kapısında kayboldu — yani iki kayıp daha aynı adrese taşındı.

%77 recall'ın güncel bileşimi:

| sebep | n | durum |
|---|---|---|
| Konak-homoloji sert filtresi (otoimmünite) | 6 | **tasarım tercihi** — savunulabilir, ödünleşme olarak raporlanır |
| Gram(+) `membrane` yapılandırma açığı | 2 | düzeltildi (CapB/ActA huniyi geçiyor) |
| PSORTb yanlış ataması (LigA) | 1 | üst-araç limiti |
| Toksisite kapısı kalibrasyon hatası | 2 (+ActA, +LLO) | **düzeltildi** — mühendislik hatası |

Yani iki özdeş ders: hem Gram listesi hem toksisite eşiği, bir aracın **kendi geçerlilik
alanı dışında** kullanılmasıydı. Yeni kural: her sert kapı, bağımsız kontrol setlerinde
(insan self k-mer + doğrulanmış IEDB epitopları + uzunluk-eşli protein setleri)
kalibre edilmeden üretimde kullanılmaz — `benchmark/calibrate_filters.py`.

Yeniden üretim: `results/miss_attribution.json` (her antijen için insan %kimlik + yakalanma durumu
+ atfedilen sebep).
