**İleriye dönük ground-truth: S. suis pilotu**

Kıyas ground-truth'umuzun %59'u Vaxign-ML'in eğitim etiketi
(`results/vaxignml_circularity.json`). Bu dairesellik, karşılaştırmanın genelleme
hakkında söyleyebileceğini sınırlıyor. Bağımsız bir set kurmanın tek yolu, o modelin
eğitim kesiti olan **31 Temmuz 2019'dan sonra** yayımlanmış ve Protegen'e uzak duran
antijenler. Bu pilot, yöntemin işlediğini bir türde gösteriyor ve hiçbir yeni pipeline
koşusu gerektirmedi: sıralar diskteki `outputs/validation_ssuis/run.json` aday
listesinden yeniden hesaplandı.

**Küratörlük (9 antijen, hepsi tam metinden okundu)**

Ana kaynak BacScan (Dong ve ark. 2024, Front Immunol 15:1392456,
doi:10.3389/fimmu.2024.1392456): *S. suis* çekirdek genomunun genom-ölçeğinde T7 fay
sergileme + immünopresipitasyon + NGS taraması. 19 yüksek-immünojenik protein, bunların
8'i BALB/c farede ölümcül i.p. challenge'a karşı %20-100 koruma. Antijenler SC19 locus
etiketleriyle verildiği için eşleme ad üzerinden değil, dizi hizalamasıyla yapıldı;
makalenin kendi verdiği üç çapa (2750=Ide-Ssuis, 5995=PrsA, 6615=SAO) SC19 annotasyonunda
doğrulandı. İkinci kaynak Ide-Ssuis'in domuz denemesi (doi:10.1016/j.jvacx.2024.100590).

Dokuzunun dokuzu BM407 proteomunda ≥%93.6 kimlikte bulundu. Protegen'e göre:

| antijen | hedef aksesyon | Protegen kimlik | sınıf | sıra /488 |
|---|---|---|---|---|
| HIP 5995 (PrsA, peptidilprolil izomeraz) | A0A0H3MYQ8 | %48.7 | homolog | 3 |
| HIP 10530 (S1C serin proteaz) | A0A0H3MYR1 | %36.0 | homolog | 12 |
| Ide-Ssuis (IgM sistein endopeptidaz) | A0A0H3MW10 | %30.8 | homolog | 26 |
| HIP 8455 (FMN bağlayan protein) | A0A0H3MXI4 | %32.8 | homolog | 45 |
| HIP 10025 (Zn/Mn ABC taşıyıcı) | A0A0H3MYH9 | %37.4 | homolog | 110 |
| HIP 2085 (MapZ ailesi) | A0A0H3MTL4 | hit yok | **bağımsız** | 66 |
| HIP 6155 (SAG1386/EF1546 yüzey-ilişkili) | A0A0H3MYM3 | hit yok | **bağımsız** | aday üretmedi |
| HIP 575 (30S ribozomal protein S13) | A0A0H3MSQ3 | hit yok | **bağımsız** | aday üretmedi |
| HIP 3685 (LD-karboksipeptidaz LdcB/DacB) | A0A0H3N4I8 | hit yok | **bağımsız** | aday üretmedi |

**Bulgu 1: dört antijen gerçekten bağımsız.** Prior-art'ı olanlar (PrsA, Ide, 10025)
düşüldükten sonra kesitten sonra ilk kez aşı adayı olarak sınanan ve Protegen'e hiç hit
vermeyen dört antijen kalıyor: 6155, 2085, 575, 3685. Yani yöntem çalışıyor, bir türde
dört bağımsız antijen üretiyor.

**Bulgu 2: bağımsız sette recall 1/4, ve nedeni ölçüldü.** Dört antijenden yalnız 2085
aday üretti (sıra 66/488). Üçünün aday üretmeme nedeni PSORTb ile ölçüldü: 6155
sitoplazmik (7.5), 575 sitoplazmik (9.97), 3685 sitoplazmik (9.2). Aday üreten 2085'in
kararı ise `unknown` (2.5), yani kapı onu yumuşak geçirmiş.

Bu, dün bulduğumuz kısıtın bağımsız bir sette kendiliğinden tekrarlaması: **bağlayıcı
kısıt lokalizasyon tahmini.** Önceki 11 kaybın atfı bizim seçtiğimiz antijenler
üzerindeydi, dolayısıyla seçim yanlılığına açıktı; bu dört antijen tasarımımızı
etkileyemeyecek kadar sonra yayımlandı ve aynı sonucu veriyor.

HIP 6155 ayrıca somut bir lokalizasyon hatası örneği: RefSeq'te "SAG1386/EF1546 ailesi
**yüzey-ilişkili** protein" olarak anote, PSORTb ona 7.5 skorla sitoplazmik diyor.

**Bulgu 3: ters yönde bir uyarı.** Protegen'e homolog olanlar iyi sıralanıyor (3, 12, 26,
45), bağımsız olan tek aday 66'da. Örneklem 9 antijen, çıkarım için az; ama yön,
dairesellik ölçümünde gördüğümüzle aynı (uzaklaştıkça sıra bozuluyor) ve bu kez
VaxForge'un kendi tarafında. Kayıt altına alınır, sonuç olarak sunulmaz.

**Bulgu 4: mevcut ground-truth'ta bir çelişki.** Enolaz (A0A0H3N5B2) bizim ssuis
koruyucu listemizde. Dumesnil ve ark. (doi:10.1186/s12917-019-2196-y, BMC Vet Res
15(1):448, Aralık 2019; **2026-10-12 düzeltmesi:** bu kayıt önceki sürümde yazarı yanlış,
"Bichet ve ark." olarak veriyordu, PubMed künyesiyle doğrulanıp düzeltildi) CD-1
fare modelinde dört farklı adjuvanla SsEno'nun ölümcül serotip 2 challenge'a karşı
**koruma sağlamadığını** gösteriyor (antikor üretilmesine rağmen). Bu, kıyas setinin
kendi belirsizliği olarak kaydedilir; antijen listeden çıkarılmadı, çelişki not edildi.

**Henüz çözülemeyenler.** Üç yayının antijenleri locus etiketi olmadan ada göre
eşlenemiyor ve tam metin erişimi gerekiyor: HK/Gtf-2/PsbP (doi:10.1016/j.micpath.2024.106759
— proteomda 10 histidin kinaz ve 6 glikoziltransferaz var, ad eşlemesi belirsiz), SBP2'
(doi:10.1016/j.rvsc.2021.05.006), V5 karışımının beş virülans faktörü
(doi:10.1093/infdis/jiz352). Ad eşlemesi yapılmadı: "suis" regex'inin 28 *Glaesserella
parasuis* kaydını yakaladığı hata tekrarlanmasın diye eşleme yalnız dizi hizalamasıyla
kabul ediliyor.

Betik: `benchmark/build_prospective_groundtruth.py`, küratörlük kaydı:
`benchmark/prospective/ssuis.json`, veri:
`benchmark/results/prospective_groundtruth_ssuis.json`.

---

## 2026-10-12 — Askıda kalan üç yayın kapatıldı

Önceki sürümün "henüz çözülemeyenler" bölümü üç yayın bırakmıştı. Üçü de kapatıldı: ikisi
çözüldü, biri erişim yokluğundan çözülemedi. Kapanışların hiçbiri ada göre eşlemeye
dayanmıyor; her biri ölçülebilir kanıta bağlandı ve
`benchmark/results/pending_antigen_resolution.json`'a yazıldı (betik:
`benchmark/resolve_pending_antigens.py`).

**Erişim denetimi ilk adım oldu.** Üç PMID'in (38906494 / 34020335 / 31287878) hiçbirinin
PMC karşılığı yok; Europe PMC üçü için de "Subscription required" diyor ve annotations API
yalnız başlık+özeti işliyor. Yani locus etiketleri **birincil** tam metinden okunamıyor.
Çözüm, aynı grupların **açık erişimli ikincil** makalelerinden kimliği doğrulanabilir
özelliklerle çıkarmak oldu.

**V5 karışımı (doi:10.1093/infdis/jiz352) — çözüldü.** Aynı grubun açık erişimli eLife
makalesi (doi:10.7554/eLife.101760, PMC12185102) Yöntemler'de bileşimi harf harf veriyor:
SrtA (aa 50-249), ScpCL (36-535), MRP (283-721), SLY (28-497), SspA (41-852), hepsi
*S. suis* A7 genomundan (NC_017622.1) klonlanmış. Kaynak genom + rezidü aralıkları
verildiği için eşleme ad tahmini değil, uzunluk kontrollü: **SrtA (n=249) ve SLY (n=497)
bildirilen aralığın üst sınırıyla birebir**, MRP'nin annotasyondaki gen adı zaten `mrp`.
A7'de iki S8-ailesi peptidaz var (SspA/ScpCL belirsizliği); *S. pyogenes* C5a peptidazı
ScpA (P15926) ile hizalama ayırdı: SSUA7_RS10845 %38.4 kimlik / 396 aa (hizalama q147-539,
bildirilen 36-535 aralığını kapsıyor) → ScpCL; kalan SSUA7_RS09010 (%25.9, n=1692) → SspA,
**eleme ile**, doğrudan değil. Yan bulgu: ScpCL A7 annotasyonunda `pseudo=true`, yani WP
aksesyonu yok, anahtar locus etiketi olmak zorunda.

**SBP2' (doi:10.1016/j.rvsc.2021.05.006) — çözüldü, dört bağımsız kontrolle.** Aynı grubun
açık erişimli AEM makalesi (doi:10.1128/aem.01267-26, PMC13599487) dört doğrulanabilir
özellik veriyor ve dördü de ölçüldü: (1) srtBCD kümesi ZY05719'da **ada göre değil genomik
komşuluktan** bulundu — üç yan yana class C sortaz (RS09890/09895/09900), küme tam beş
pilin içeriyor; (2) beş pilinden **yalnız** RS11000 YPKN pilin motifini taşıyor; (3) RS11000
*S. pneumoniae* RrgB'nin **N-terminal** bölgesine (rez. 36-243) %31.1 kimlikte hizalanıyor —
yayın %29 bildiriyor; (4) n=286 (~31-33 kDa), yayındaki 33 kDa Western bandıyla uyumlu.
Eşi de doğrulandı: RS10995 RrgB'nin **C-terminal** bölgesine %33.1 (yayın %35) → SBP2''.
Sonuç: SBP2' = ZY05719_RS11000 / WP_012027978.1.

**HK / Gtf-2 / PsbP (doi:10.1016/j.micpath.2024.106759) — çözülemedi, kapatıldı.** Özette
ne aksesyon, ne locus etiketi, ne de kaynak suş var; tam metin abonelik arkasında ve grubun
açık erişimli diğer makaleleri bu üç proteini aksesyonla vermiyor. Belirsizliğin büyüklüğü
hedef proteomda **sayıldı**: histidin kinaz 10 eşleşme, glikoziltransferaz 6 eşleşme → ad
eşlemesi imkânsız. PsbP için ad tek bir proteine uyuyor (A0A0H3MUT8), **ama** kaynak
suş/aksesyon bilinmediği için ölçütün gerektirdiği dizi hizalaması yapılamıyor; ada göre
kabul **edilmedi**. Kurumsal tam metin erişimiyle üçü de bir oturumda çözülür.

### Genişletilmiş setin sonuçları

Küratörlük 9 → **15 antijen**. Beş V5 bileşeni `mixture=true` ile işaretlendi: karışım
olarak uygulandıkları için koruma tek bileşene atfedilemez, bu yüzden prospektif
tek-antijen sayımına **hiçbir koşulda** girmiyorlar (ölçüt listesine bu kural eklendi).
Hedef proteomda 12, prospektif kullanılabilir **4** — yani ana sayı değişmedi, **recall
1/4 duruyor**. Değişen, nedenlerin çözünürlüğü:

| yeni antijen | hedef | Protegen | sınıf | sıra |
|---|---|---|---|---|
| SrtA (V5) | A0A0H3MVC0 | %61.9 | homolog | 30 |
| SLY (V5) | A0A0H3MU02 | %50.9 | homolog | 9 |
| ScpCL (V5) | — (%29.8) | %55.6 | homolog | hedefte yok |
| MRP (V5) | — (hit yok) | %28.3 | bağımsız | hedefte yok |
| SspA (V5) | A0A0H3N6J4 | %27.4 | **bağımsız** | **aday üretmedi** |
| SBP2' | — (hit yok) | hit yok | bağımsız | hedefte yok |

**Bulgu 5: ikinci bir kayıp mekanizması ölçüldü — ve lokalizasyon değil.** SspA hedef
proteomda %100 kimlikte var, PSORTb ona `cell_wall 9.98` diyor, yani Gram(+) izinli
listesinde **geçerli** bir yer. Yine de aday üretmemiş. Nedeni iddia edilmedi, huni o
protein için koşulup trace okundu (`benchmark/measure_funnel_dropout.py` →
`results/funnel_dropout.json`): kapı **TM sayısı**, `2 TM > 1`. Lokalizasyon ve
self-homoloji kapılarını geçiyor. Kontrol proteini (suilizin, koşuda sıra 9) aynı
alt-kümede geçti → alt-küme ölçümünün tam proteomla aynı kararı verdiği sınandı.
Şimdiye kadarki tüm kayıpların bağlayıcı kısıtı lokalizasyondu; bu, ondan farklı
**ikinci** bir mekanizma. Kapsamı ölçüldü, kapı değiştirilmedi — bkz. FILTER_CALIBRATION.md.

**Bulgu 6: hedef suş, prospektif setin tavanını belirliyor.** Altı yeni antijenden üçü
BM407 proteomunda **yok**: MRP ve SBP2' hiç hit vermiyor, ScpCL yalnız %29.8 (aynı türde
"aynı protein" eşiği %90). SBP2' için bunun bağımsız bir açıklaması var: AEM makalesi
srtBCD kümesinin yalnız serotip 14 ve virülan serotip 2 izolatlarında bulunduğunu
söylüyor. MRP'nin suş-değişkenliği zaten literatürde bilinen. Sonuç yöntemseldir: ileriye
dönük set kurarken antijenin yayımlandığı suş ile bizim koştuğumuz suş arasındaki fark,
kullanılabilir antijen sayısının üst sınırını baştan kısıyor — küratörlükten değil, suş
seçiminden gelen bir kayıp.

**Bulgu 7: dairesellik, klasik antijenlerde daha da belirgin.** V5'in Protegen'e
hizalanan üç bileşeni (SrtA %61.9, ScpCL %55.6, SLY %50.9) setin en yüksek örtüşme
değerleri arasında. Yani 2019 öncesinin en çok çalışılmış *S. suis* antijenleri
Vaxign-ML'in eğitim etiketlerine en yakın olanlar — dairesellik ölçümünün (results/
vaxignml_circularity.json) beklediği yön. Kayıt altına alınır, bağımsız sonuç sayılmaz.
