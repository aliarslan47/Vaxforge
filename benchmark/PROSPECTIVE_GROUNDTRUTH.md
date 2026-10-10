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
koruyucu listemizde. Bichet ve ark. (doi:10.1186/s12917-019-2196-y, Aralık 2019) CD-1
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
