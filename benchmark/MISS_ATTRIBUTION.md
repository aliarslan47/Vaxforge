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

## 1. Konak-homoloji filtresi — 6/11, **YENİDEN SINIFLANDIRILDI (2026-10-10): GRANÜLARİTE HATASI**

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

### Yeniden sınıflandırma (2026-10-10): "tasarım tercihi" kısmen yanlıştı

Yukarıdaki "bilerek atıyor" yorumu kapının **varlığını** haklı çıkarıyor ama
**granülaritesini** sorgulamıyordu. Ölçüldüğünde ikisi ayrıştı
(`benchmark/measure_homology_granularity.py`, S. suis, k=9, insan Swiss-Prot):

| Antijen | protein-bütünü kimlik | 9-mer | self ile TAM eşleşen | oran | güvenli k-mer |
|---|---|---|---|---|---|
| gnd | %52.3 | 467 | 14 | %3.00 | 453 |
| Enolaz | %49.8 | 427 | 8 | %1.87 | 419 |
| GAPDH | %43.9 | 328 | 7 | %2.13 | 321 |
| DnaJ | %37.0 | 370 | 5 | %1.35 | 365 |
| Suilizin (yakalanan) | %0 | 489 | 0 | %0 | 489 |

Protein-bütünü kimlik, epitop-düzeyi otoimmünite riskinin vekili değil: enolazda
8 riskli 9-mer'i elemek için 419 self-olmayan epitop da atılıyordu. Yani kayıp
"güvenlik için ödenen bedel" değil, **gereğinden pahalı ödenmiş bir bedeldi**.

Düzeltme: protein düzeyi homoloji YUMUŞAK (işaret + skor ×0.75; eleme yalnız
≥%90, yani konak proteininin kendisi), epitop düzeyi mimikri SERT (konak self ile
tam eşleşen epitop elenir; kalibrasyon: doğrulanmış IEDB epitoplarında %1.3 YP).
Otoimmünite güvenliği kaybolmuyor, riskin gerçekte bulunduğu granülariteye taşınıyor.
Korunan ihtiyat: BLOSUM-yakın eşleşmeler elenmiyor ama skor cezası alıyor (TCR
dejenerasyonu), ve protein-düzeyi işaret raporda görünür kalıyor.

**Bu 6 antijenin geri döneceği iddiası DEĞİLDİR.** Artık huniden geçebilirler;
lokalizasyon/TM/antijenite kapıları ve ×0.75 ceza hâlâ önlerinde. Gerçek etki
yeniden koşumla ölçülecek.

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

## 3. Lokalizasyon tahmininin genel sınırı, 1/11 — **YENİDEN SINIFLANDIRILDI (2026-10-10)**

LigA (*Leptospira*) → PSORTb `cytoplasm` skor 9.97. Literatürde bilinen yüzey lipoproteini.
İlk hâlinde buna "**PSORTb yanlış ataması**" dendi, yani tek bir üst-aracın kusuru.
**Bu nitelendirme fazla iyimserdi.**

İkinci, mimari olarak tamamen farklı bir tahminciyle sınandı: DeepLocPro 1.0 (Moreno ve
ark. 2024, doi:10.1093/bioinformatics/btae677) — ESM-2 dil modeli tabanlı, prokaryotlara
özgü, PSORTdb 4.0'ın deneysel-doğrulanmış verisiyle eğitilmiş. Yayımlanan kıyasta PSORTb
3.0'ı tüm gruplarda geçiyor ve "hücre duvarı ve yüzey" sınıfındaki fark uçurum:
Gram(+) MCC 0.28 vs −0.02, Gram(−) 0.80 vs 0.00 — yani PSORTb o sınıfta pratikte rastgele.
Bu istatistiğin bizim kayıplarımıza tercüme olup olmadığını ölçtük
(`benchmark/measure_localization_agreement.py`, veri: `results/localization_agreement.json`):

| Antijen | Gram | PSORTb | DeepLocPro | güven |
|---|---|---|---|---|
| Suilizin *(yakalanan kontrol)* | + | extracellular 9.73 | Extracellular | 0.999 |
| Enolaz | + | cytoplasm 9.97 | Cytoplasmic | 0.988 |
| GAPDH | + | cytoplasm 9.97 | Cytoplasmic | 0.875 |
| gnd | + | cytoplasm 7.50 | Cytoplasmic | 0.845 |
| DnaJ | + | cytoplasm 9.97 | Cytoplasmic | 0.987 |
| CapB | + | membrane 10.0 | Cytoplasmic Membrane | 0.939 |
| ActA | + | membrane 10.0 | Cytoplasmic Membrane | 0.656 |
| **LigA** | − | **cytoplasm 9.97** | **Cytoplasmic** | **0.995** |

**Üst karar uyumu 8/8.** (PSORTb'nin Gram(+) `membrane` etiketi ile DeepLocPro'nun
`Cytoplasmic Membrane` etiketi aynı sınıftır — Gram(+)'da dış membran yoktur; betik
etiketleri normalize ederek sayar.)

Sonuç: LigA'yı DeepLocPro da 0.995 güvenle sitoplazmik veriyor. Dolayısıyla bu **tek bir
aracın yanlış ataması değil, dizi-temelli lokalizasyon tahmininin genel sınırı** — ve bu
kayıp aracı değiştirerek kapatılamaz. Aynı ölçüm, 1. maddedeki moonlighting sonucunu da
tek araçlı bir kabulden iki araçlı bir bulguya yükseltiyor.

**Yan sonuç (mühendislik kararı).** Lokalizasyon katmanını DeepLocPro'ya taşımak
planlanmıştı; bu ölçüm gerekçeyi ortadan kaldırdı. Yayımlanmış MCC farkı gerçek ama bizim
kayıp sınıfımıza tercüme olmuyor. Aracı eklemeden önce ölçmek, gereksiz bir kod
değişikliğini ve 10 türün yeniden koşulmasını engelledi. DeepLocPro yine de yerel olarak
kuruldu (`tools/deeplocpro_src`, CC BY-NC-SA 4.0) ve bu ölçüm yeniden üretilebilir.

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
| Konak-homoloji kapısı — yanlış granülarite | 6 | **düzeltildi** — protein düzeyi yumuşak, epitop düzeyi sert; kapının varlığı tasarım tercihi, protein-bütünü eleme mühendislik hatasıydı |
| Gram(+) `membrane` yapılandırma açığı | 2 | düzeltildi (CapB/ActA huniyi geçiyor) |
| Lokalizasyon tahmininin genel sınırı (LigA) | 1 | **yeniden sınıflandırıldı** — PSORTb'ye özgü değil; DeepLocPro 1.0 da %99.5 güvenle sitoplazmik diyor (8/8 uyum) |
| Toksisite kapısı kalibrasyon hatası | 2 (+ActA, +LLO) | **düzeltildi** — mühendislik hatası |

Yani **üç** özdeş ders: Gram listesi, toksisite eşiği ve konak-homoloji granülaritesi —
üçü de bir ölçümün **kendi geçerlilik alanı dışında** kullanılmasıydı (yanlış organizma
sınıfı, yanlış uzunluk sınıfı, yanlış granülarite). Bu sayımla 11 kaybın 10'u mühendislik
kaynaklı çıkıyor; ilk raporun "mühendislik hatası yok" sonucu tümüyle geri alınmıştır. Yeni kural: her sert kapı, bağımsız kontrol setlerinde
(insan self k-mer + doğrulanmış IEDB epitopları + uzunluk-eşli protein setleri)
kalibre edilmeden üretimde kullanılmaz — `benchmark/calibrate_filters.py`.

Yeniden üretim: `results/miss_attribution.json` (her antijen için insan %kimlik + yakalanma durumu
+ atfedilen sebep).

---

## 2026-10-12 — Bağımsız sette dördüncü mekanizma: TM kapısı

Bu belgenin atfı, bizim seçtiğimiz 48 antijen üzerindeydi. İleriye dönük (post-2019-07-31)
bağımsız set genişletilince (PROSPECTIVE_GROUNDTRUTH.md, 9 → 15 antijen) **bu listede
olmayan** bir mekanizma ölçüldü.

**SspA (subtilizin-benzeri S8 peptidaz, *S. suis*).** Hedef proteomda %100 kimlikte var.
PSORTb `cell_wall 9.98` — Gram(+) izinli listesinde geçerli bir yer, yani §2 ve §3'teki
lokalizasyon mekanizmalarının **ikisi de değil**. Self-homoloji %0.0, konak-homoloji
kapısı da değil (§1 değil). Toksisite/epitop aşamasına hiç gelmiyor (§4 değil). Huni o
protein için koşulup trace okundu: kapı **TM sayısı**, `2 TM > 1`.

| mekanizma | bu belgedeki sayım | bağımsız sette |
|---|---|---|
| konak-homoloji granülaritesi | 6 | 0 |
| lokalizasyon (Gram açığı + genel sınır) | 3 | 3 (6155, 575, 3685 — hepsi sitoplazmik) |
| toksisite kalibrasyonu | 2 | 0 |
| **TM kapısı** | **0 — bu belgede yok** | **1 (SspA)** |

İki ders. Birincisi doğrulama: lokalizasyonun bağlayıcı kısıt olduğu bulgusu, seçim
yanlılığına kapalı bir sette (üç kayıptan üçü) kendiliğinden tekrarlıyor. İkincisi yeni:
**48 antijenlik ilk set, var olan bir mekanizmayı göstermeye yetmemiş.** TM kapısı o sette
hiç bağlayıcı olmamıştı, çünkü o antijenlerin hiçbiri TM>1 değildi (tm_gate_scope.json'da
16 antijenin 15'i TM≤1). Bağımsız set, atfın kendi kapsamını da sınıyor.

Kapı **değiştirilmedi**: tek antijen, gösterilmiş bir yanlış kalibrasyon değil. Kapsam
ölçülüp kalibrasyon kuyruğuna alındı (FILTER_CALIBRATION.md, 2026-10-12), çünkü kural
aynı kalıyor — ölçmeden kapıya dokunulmaz.

Yeniden üretim: `benchmark/measure_funnel_dropout.py` → `results/funnel_dropout.json`,
`benchmark/measure_tm_gate.py` → `results/tm_gate_scope.json`.
