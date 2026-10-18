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

---

## 2026-10-12 — TM kapısı (`max_helices`) kalibrasyon kuyruğuna alındı, DEĞİŞTİRİLMEDİ

Ölçüm: `benchmark/measure_tm_gate.py` → `results/tm_gate_scope.json`; tek-protein teşhisi
`benchmark/measure_funnel_dropout.py` → `results/funnel_dropout.json`.

**Nasıl gündeme geldi.** İleriye dönük ground-truth'ta SspA (subtilizin-benzeri S8
peptidaz, *S. suis* A7'den V5 karışımının bileşeni) hedef proteomda %100 kimlikte var,
PSORTb ona `cell_wall 9.98` diyor — Gram(+) izinli listesinde geçerli bir yer — ama koşuda
aday üretmedi. Huni o protein için koşuldu: lokalizasyon ve self-homoloji kapılarını
geçiyor, **TM kapısında düşüyor** (`2 TM > 1`, tmhmm_local). Kontrol proteini (suilizin,
koşuda sıra 9) aynı alt-kümede geçti, yani ölçüm tam proteomla aynı kararı veriyor.
Bu, şimdiye kadar atfedilen tüm kayıpların mekanizmasından (lokalizasyon) **farklı**.

**Hipotez (sınandı, varsayılmadı).** Gram(+) hücre-duvarına çapalı yüzey proteinleri
N-ucunda sinyal peptidi, C-ucunda LPXTG/membran çapası taşır. Hidrofobik sayıcı bu ikisini
iki ayrı sarmal sayarsa, "çok-TM = zor ifade edilir" gerekçesiyle yazılmış kapı tam olarak
hedeflediğimiz sınıfı eler. SspA'nın profili (n=1692, sinyal peptidi(+), hücre duvarı,
TM=2) bu desene uyuyor.

**Kapsam ölçüldü (*S. suis* BM407, 1922 protein, tmhmm_local + SignalP ikisi de GERÇEK):**

| havuz | n | TM > 1 | oran |
|---|---|---|---|
| tüm proteom | 1922 | 367 | %19.1 |
| sinyal-peptidi(+) (yüzeye yönelen sınıf) | 72 | 16 | **%22.2** |
| deneysel koruyucu/prospektif antijenler | 16 | 1 (SspA) | **%6.2** |

**Karar: kapı değiştirilmedi.** Karar sayısı son satır: deneysel antijenlerin yalnız biri
düşüyor, yani kapının *gösterilmiş* bir yanlış kalibrasyonu yok — SspA sınıf-dışı tek bir
vaka olabilir. Buna karşılık kapı, yüzeye yönelen havuzun %22'sini kesiyor; bu, yanlış
kalibrasyon kanıtı değil ama görmezden gelinecek bir bedel de değil. Bu yüzden kapı
FILTER_CALIBRATION kuyruğuna alındı: değişiklikten önce koşulması gereken deney, bu
belgenin üstündeki desenin aynısı — **uzunluk- ve sınıf-eşli kontrol setleri**: pozitif =
deneysel doğrulanmış Gram(+) hücre-duvarı/yüzey antijenleri, negatif = gerçek çok-geçişli
integral membran proteinleri (taşıyıcılar, permeazlar). Ölçülecek şey, `max_helices`'in
2'ye çıkarılmasının (ya da "sinyal peptidi(+) ise ilk sarmalı saymama" kuralının) pozitif
sette kazandırdığı ile negatif sette kaybettirdiği. Tek bir türün 16 antijeni bu kararı
vermeye yetmez.

Bu, DeepLocPro çapraz denetiminde verilen kararın aynısı: **ölç, sonra değiştirip
değiştirmeyeceğine karar ver.** Ölçmeden kapıya dokunmak, 10 türün yeniden koşulmasını
gerektiren ve gerekçesi olmayan bir değişiklik olurdu.

## 2026-10-18 — TM kapısı kalibre EDİLDİ, sinyal-duyarlı sayıma geçildi

Kuyruktaki deney koşuldu. Kapı eşiği (`max_helices: 1`) **değişmedi**; değişen, neyin
sarmal sayıldığı: sinyal peptidi(+) proteinde ilk sarmal artık salgı sinyali sayılıyor,
membran çapası değil.

**Kontrol setleri.** Pozitif = 5 Gram(+) türün (ssuis, listeria, saureus, mbovis,
banthracis) deneysel doğrulanmış koruyucu antijenleri, GERÇEK PSORTb ile lokalizasyon
kapısını geçenler (n=25). Negatif = aynı proteomlardan **yalnız anotasyonla** seçilmiş
çok-geçişli integral membran proteinleri (permeaz, ABC/MFS taşıyıcı, simport/antiport,
PTS, efluks), uzunluk çeyreklerine eşlenmiş (n=102). Negatifleri TM sayısıyla seçmek
ölçümü döngüsel yapardı; seçim yapılmadı. Seed 1307, araçların üçü de gerçek
(tmhmm_local + SignalP + PSORTb Gram+).

**Ön kayıtlı karar ölçütü** (sonuçlara bakılmadan yazıldı): mevcut kural yerine bir
kural ancak duyarlılığı artırır VE özgüllükten en çok 5 yüzde puanı kaybederse önerilir;
birden fazlası uyarsa özgüllüğü en yüksek olanı.

| kural | duyarlılık | özgüllük |
|---|---|---|
| R1 mevcut, `TM ≤ 1` | 24/25 (%96.0) | 89/102 (%87.3) |
| **R3 sinyal-duyarlı, `etkin TM ≤ 1`** | **25/25 (%100)** | **89/102 (%87.3)** |
| R2 `TM ≤ 2` | 25/25 (%100) | 88/102 (%86.3) |
| R4 sinyal-duyarlı, `etkin TM ≤ 2` | 25/25 (%100) | 88/102 (%86.3) |

**Seçilen: R3.** Mevcut kuralı kesin domine ediyor — kazandığı pozitifi özgüllükten
hiçbir şey vermeden kazanıyor. Eşiği gevşeten R2/R4 aynı duyarlılığa çıkıyor ama bir
puan özgüllük veriyor, yani gereksiz.

**Mekanizma öngörülmüştü, ölçüm doğruladı.** Gram(+) duvar-çapalı yüzey proteini
N-ucunda sinyal peptidi, C-ucunda LPXTG/membran çapası taşır; sayıcı ikisini iki ayrı
sarmal sayınca "çok-TM = zor ifade edilir" gerekçeli kapı tam hedeflenen sınıfı eliyordu.
Kapının elediği tek pozitif, hipotezi doğuran proteinin kendisi: *S. suis* **SspA**
(`A0A0H3N6J4`, TM=2, SignalP 0.77, PSORTb `cell_wall`).

**Dürüst sınır.** Duyarlılık kazancı **tek antijene** dayanıyor (n=1). Buna karşılık
özgüllük kaybı **sıfır**, yani R3'ün ölçülen hiçbir eksende bedeli yok. Benimseme
gerekçesi "muhtemelen daha iyi" değil, "hiçbir şey kaybettirmiyor".

**Ölçüt sapması, kayda geçiyor.** Pozitif sete `mixture=true` bileşenleri **dahil
edildi** — `PROSPECTIVE_GROUNDTRUTH.md`'deki kural ise onları dışlar. O kural *recall
atfı* içindir: koruma beş bileşenli karışımdan geldiğinde tek proteine atfedilemez.
Buradaki soru atıf değil **yapı**: sayıcı bu yapısal sınıfı yanlış mı sayıyor?
Deneysel üretilmiş, PSORTb'nin duvar/yüzey dediği bir karışım bileşeni o sınıfın meşru
üyesidir. İlk koşum bileşenleri dışlamıştı ve sonuç şuydu: pozitif setteki hiçbir
proteinin TM'si 1'den büyük değil, dört kural da 23/23 — **ölçüm kendi sorusunu
göremiyordu**. O koşum silinmedi: `results/tm_gate_calibration.pre_mixture.json`.
Her iki tabaka raporlanıyor; karışım hariç tabakada R1 ile R3 aynı (%100), yani
değişiklik o tabakada hiçbir şeyi bozmuyor.

**Güç sınırı.** PSORTb 11 aday pozitifi izinli listenin dışında buldu: 6'sı kanonik
sitoplazmik (enolaz, GAPDH, DnaJ, gnd, ribozomal S13, LdcB — moonlighting yüzey
iddiaları dizi-temelli tahminciyle doğrulanamıyor), 5'i `unknown` (Ide-Ssuis,
HIP 2085/8455/10530, HIP 6155). `unknown` olanlar "sitoplazmik" değil, PSORTb'nin
**kararsız** kaldığı vakalar; pozitif sete alınmadılar ve bu, ölçümün gücünü düşürüyor.
Dizi-temelli lokalizasyonun genel sınırı (LigA vakasının aynısı), bu kapıya özgü değil.

Uygulama: `config/thresholds.yaml` → `transmembrane.signal_peptide_aware: true`
(false mevcut-öncesi davranışa döner), `vaxforge/funnel.py` etkin TM sayımı +
trace'te `tm_helices_effective`. Regresyon: `test_tm_gate_is_signal_peptide_aware`.
Koşum: `benchmark/calibrate_tm_gate.py`, veri: `results/tm_gate_calibration.json`.

## 2026-10-18 — Phobius ile bağımsız çapraz-doğrulama: yama DEĞİŞMEDİ

**Soru.** Yukarıdaki düzeltme elle yazılmış bir yamadır: "sinyal peptidi(+) ise ilk
sarmalı düş". Phobius ise aynı problemi **modelleyerek** çözer — sinyal peptidini ve
TM topolojisini tek HMM'de ortak çözümler, yani bir N-ucu bölgesinin salgı sinyali mi
membran çapası mı olduğuna kendisi karar verir (Käll, Krogh & Sonnhammer,
*J Mol Biol* 338:1027-1036, 2004; PMID 15111065). Phobius bizim yamadan daha iyi mi?

**Kurulum.** Phobius 1.01 yerel kuruldu (`tools/phobius-1.01/`, 0.6 MB, Perl,
`decodeanhmm` 64-bit derlemesi; akademik kullanım lisansı). Kontrol setleri
**yeniden seçilmedi** — Faz 1'in tam aynı 25 pozitif / 102 negatifi
(`results/tm_gate_calibration.json`). Araç değişirken setin de değişmesi
karşılaştırmayı yorumlanamaz hale getirirdi.

**ÖN KAYITLI ölçüt** (sonuçlara bakılmadan yazıldı): Phobius üretimdeki kuralın
(T3) yerine ancak duyarlılığı ve özgüllüğü **ikisi de düşmeden, en az birinde kesin
artışla** geçerse önerilir. Eşitlik çıkarsa öneri **DEĞİŞTİRME**: aynı sonucu veren
yeni bir dış bağımlılık kurulum, bakım ve atıf yükü getirir, karşılığını vermez.

| kural | duyarlılık | özgüllük |
|---|---|---|
| T1 TMHMM `tm ≤ 1` (kalibrasyon öncesi) | 24/25 %96.0 | 89/102 %87.3 |
| **T3 TMHMM sinyal-duyarlı `≤ 1` (ÜRETİMDEKİ)** | **25/25 %100** | **89/102 %87.3** |
| P1 Phobius `tm ≤ 1` (doğal) | 25/25 %100 | 89/102 %87.3 |
| P2 Phobius `tm ≤ 2` | 25/25 %100 | 88/102 %86.3 |
| P3 Phobius + aynı yama | 25/25 %100 | 89/102 %87.3 |

**Sonuç: DEĞİŞTİRME.** Phobius doğal haliyle (P1) üretimdeki kuralın sonucunu
birebir veriyor — yalnız toplamlar değil, **127 proteinin her birinde kapı kararı
aynı** (T3 ≠ P1 olan protein sayısı: 0). Ön kayıtlı ölçüt "kesin artış" istiyor,
artış yok → Phobius üretime alınmadı. Pipeline değişmedi.

**Ama ölçüm boşa gitmedi: mekanizma bağımsız araçla doğrulandı.** Faz 1'in dayandığı
iddia şuydu — "TMHMM, Gram(+) yüzey proteininin N-ucundaki salgı sinyalini sarmal
sayıyor". Kritik vaka SspA (`A0A0H3N6J4`): TMHMM 2 sarmal diyor, **Phobius aynı
proteinde 1 sarmal + sinyal peptidi(+)** diyor. Yani ortak modelli bağımsız bir
tahminci, bizim elle düştüğümüz sarmalın gerçekten salgı sinyali olduğunu söylüyor.
Aynı desen iki negatif kontrolde de görülüyor (glutamin ABC taşıyıcı,
`A0A0H3MV05`/`A0A0H3N690`: TMHMM 5 → Phobius 3 + sinyal(+)) — ikisi de her iki kuralda
doğru şekilde elenmeye devam ediyor.

**Yan gözlemler (kayda geçiyor, karar değiştirmiyor).**
- TM sayısı birebir uyum yalnız 85/127; Phobius 20 proteinde daha az, 22'sinde daha
  çok sarmal sayıyor. Buna rağmen **kapı kararları birebir aynı** — eşik 1'de olduğu
  için sayımdaki orta-aralık farkları karara yansımıyor. Eşik gevşetilirse (P2)
  bu farklar görünür hale geliyor ve özgüllük 1 puan düşüyor.
- Sinyal peptidi kararında 14/127 fark var, **hepsi aynı yönde**: Phobius SP(+)
  diyor, SignalP-5.0 demiyor (pozitif sette 21/25 vs 17/25). Phobius sinyal
  çağrısında daha liberal. Bu, kapıyı etkilemedi ama SignalP'in ayrı bir kullanımı
  (kesim yeri) varsa akılda tutulmalı.
- Hız: Phobius tam ssuis proteomunu (1922 protein) tek çekirdekte 529 s / 11 MB RSS'te
  bitirdi. Yani hız gerekçesiyle de geçiş gerekmiyor; funnel adımı zaten koşumun
  yalnız %8.5'i (darboğaz MHC paneli, %86.6).

**Phobius kurulu kalıyor** (3.1 MB) — kaldırmak için sebep yok ve bir sonraki TM
sorusunda hazır ikinci görüş. Pipeline onu **çağırmıyor**.
Koşum: `benchmark/compare_phobius_tm.py`, veri: `results/phobius_vs_tmhmm.json`.

**Not — ölçüm hatası ve düzeltmesi.** İlk koşumda `phobius.pl -short` çıktısının SP
kolonu sayı sanılmış (`int()`), oysa `Y`/`0` yazıyor; sonuç olarak **sinyal peptidi(+)
proteinlerin hepsi sessizce düşmüştü** (127 → 94, pozitifler 25 → 4) ve ölçüm yine
tam da sormak istediği sınıfı kaybetmişti. Faz 1'deki `mixture` hatasının aynı
ailesinden: "setin kendi sorusunu içerdiğini doğrula". Parser artık beklenmedik
kolon değerinde sessizce atlamıyor, hata veriyor.
