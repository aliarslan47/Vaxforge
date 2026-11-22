# VaxForge — DURUM (nerede kaldık)

> Forge ailesi konvansiyonu: bu dosya **ilk okunacak** dosyadır. Kod tarihçesi git'te,
> bilimsel ayrıntı `benchmark/` altında, faz günlüğü `ROADMAP.md`'de.

**Son güncelleme:** 2026-11-22
**Son commit:** `c756f48` — konservasyon + alerjenite atfı + pLDDT + tool_status
**Bu turda (beşinci):** PDF rapor düzeni (her başlık kendi sayfasında, iç içe
girmiş tablo başlıkları, satır ortasından bölünen allel/PMID), IFN-γ'nın
"ölçülmedi"yi "✗" sayması, PDF'te hiç olmayan üç bölüm, raporu koşu
tekrarlanmadan yeniden üretme — aşağıda.
**Faz:** makale yazımı / son rötuş. Pipeline tamam, benchmark tamam, taslak hazır.
**Araç durumu:** 21 araçtan 20'si kurulu. DeepLoc-2.1 ve NetCTL-1.2 **akıştan
tamamen çıkarıldı** (kullanıcı kararı 2026-11-22): ikisi de hiç koşmuyordu, adları
hiçbir katmanda geçmiyor. Lokalizasyon = PSORTb (bakteri) / etiketli heuristik;
antijen işleme = IEDB processing API (tek omurga).

## 2026-11-22 (beşinci tur) — rapor kusurları (Ali'nin PDF incelemesi)

Ali, arayüzden başlattığı koşunun PDF'ini (`run_2026-11-22T20-34-48`) okuyup üç
şey istedi: her başlık kendi sayfasından başlasın, iç içe girmiş tablo başlıkları
düzeltilsin, raporun hataları bulunup kaldırılsın — **koşu tekrarlanmadan**.
Rapor `report_pdf.build()`'den yeniden üretildi; hat hiç çağrılmadı.

**Ölçüm (eski PDF, pypdf ile metin konumları):** tablo başlıklarında **2 gerçek
çakışma** — 'Kaynak antijen' başlığı komşu sütunun başlığının üstüne s1'de
13.7 pt, s2'de 2.4 pt biniyordu. Yeni PDF'te **0**.

Bulunan ve düzeltilen kusurlar:

1. **Her bölüm kendi sayfasında değildi** (9 sayfaya 15 bölüm sıkışmıştı; başlık
   sayfa dibinde kalıp tablosu bir sonraki sayfada açılıyordu). `sect()` artık
   `PageBreak` + **otomatik numara** veriyor; elle yazılan `2b/2c/2d/2e`
   etiketleri kalktı (bir bölüm atlandığında numara atlıyordu). 18 sayfa, 15
   bölüm, numaralar kesintisiz.
2. **Tablo başlıkları çıplak `str` idi.** reportlab çıplak metni hücreye
   sarmaz, sınırdan taşırır → komşu başlığın üstüne biner. Artık BAŞLIK DA
   dahil her hücre `Paragraph`; ayrıca hücre boşluğu 6+6 pt'den 2.2 pt'ye indi.
3. **Sütun genişlikleri içeriğe göre ölçülmüyordu.** `_min_widths`/`_fit_widths`:
   hiçbir sütun, bölünemez en uzun simgesinden dar kalamaz; toplam genişlik
   çerçeveyi tam doldurur; tek sütun çerçevenin %45'inden fazlasını yemez.
4. **Satır ortasından bölünen simgeler.** `wordWrap="CJK"` her karakterde
   bölüyordu: allel `BoLA-DRB3_0120`+`1`, PMID `1116366`+`5`, dizi
   `PAPRVTVYVDPPA`+`YP`. Artık boşluklu hücre KELİME sınırında, boşluksuz
   (peptit/DNA) hücre her yerde bölünüyor + (3) genişliği garanti ediyor.
5. **IFN-γ: ölçüm yokluğu olumsuz bulgu olarak basılıyordu.** `ifngamma.model_host`
   yalnız insan/fare için model döndürür; konak sığır olduğundan metrik HİÇ
   ölçülmedi. Rapor her HTL satırına `✗` yazıyordu ("indüklemiyor" diye okunur)
   ve ★ ölçütü `is True` istediği için **209 HTL epitobunun hiçbiri ★ alamıyordu**
   (tabloda tek yıldız yoktu). Artık üç değerli: `ifn_ok ∈ {True, False, None}`,
   hücre `ölç.`, altında açıklama notu, ★ yalnız ÖLÇÜLÜP negatif çıkınca engellenir.
   Düzeltmeden sonra gösterilen 16 HTL adayının 13'ü ★ alıyor. (Dürüstlük
   kusurunun 10. varyantı.)
6. **PDF'te hiç olmayan üç bölüm** (HTML'de vardı): yapısal değerlendirme
   (ESMFold/pLDDT — bu koşuda pLDDT 47.7 "çok düşük, düzensiz model"; yüzey % ve
   7 konformasyonel yama artık "düşük güven" etiketli ve yama başına pLDDT var),
   in-siliko klonlama, araç durumu (21 araç + "kurulu değil"/"doğrulama düştü"
   ayrımı). Eksik bölüm de bir dürüstlük kusuru: PDF'i okuyan kişi pLDDT'yi hiç
   görmüyordu.
7. **MEV alerjenitesi gerekçesiz "Allergen" yazıyordu.** Artık atıf tablosu:
   üç 6-mer eşleşmesinin üçü de linker/epitop **birleşme yerinde** doğuyor
   (`AAAKSK`/`AAKSKL` ← EAAAK+SKLGAARGY, `AAYTPA` ← AAY+TPAPRVTV), hiçbir epitop
   kendiliğinden alerjen değil, linker değişikliğiyle giderilebilir. Düşük
   karmaşıklık nedeniyle dışlanan 20 epitop da gerekçesiyle tabloda.
8. **Kaynakçada kurulu olmayan araç işaretsizdi** (HTML `[ARAÇ YOK]` yazıyordu,
   PDF yazmıyordu) — artık ikisi aynı.

**Raporu koşu tekrarlanmadan üretme** (bu turun kalıcı kazancı): `peptide_io`
her koşuda `peptides.json` yazıyor (rapor `Peptide` nesnelerinden üretiliyor,
`run.json` yalnız CSV sütunlarını taşıyordu), `scripts/rebuild_report.py` ise
`run.json` + `peptides.json`'dan HTML+PDF'i yeniden üretiyor. Bundan sonra bir
rapor hatası 1 saatlik BoHV-1 koşusu değil, saniyeler sürüyor. (Bu koşunun
`peptides.json`'u, dosya eklenmeden önce üretildiği için `candidates.csv` +
`candidates_full.xlsx`'ten geri kuruldu; doğrulama: eski ve yeni PDF'in veri
simgeleri birebir aynı — tek fark, eskiden ortadan bölünmüş olan simgelerin
artık bütün olması.)

**Bölüm sırası (Ali, aynı gün, ikinci istek):** REFERANSLAR EN SONDA; biyolojik
kapsam/sınırlamalar eşiklerden ÖNCE; IEDB taraması tip-tablolarının HEMEN
ARDINDAN. Yeni sıra: Özet · En iyi adaylar · anchor motifi · tipe göre epitoplar ·
**IEDB** · MEV · immünojenisite · yapı · klonlama · araç durumu · popülasyon ·
**biyolojik kapsam** · **eşikler** · yöntemler · **referanslar**. HTML rapor da
aynı sıraya çekildi (iki çıktı arasındaki fark denetlenmeli — bkz. 11. varyant).
Not: eski kod yorumundaki "Eşikler + Yöntemler DAİMA kaynaklardan SONRA gelir"
kuralı bu istekle değişti; yöntemler tablosundaki [n] artık ileri atıf yapıyor
(makale düzeninde olağan).

Testler: `tests/test_rapor_duzeni.py` (10 test) — sütun genişliği güvenceleri,
IFN-γ üç değerliliği, gerçek PDF üretip her başlığın sayfa başında olduğunun ve
numaraların atlamadığının okunması, `peptides.json` gidiş-dönüşü.

## 2026-11-22 (dördüncü tur) — kalan dört biyolojik ölçüt + bir kayma hatası

Ali'nin emri: "diğerleri biyolojik olarak anlamlı şekilde çözülmeden yeni run
benim emrim dışında yapılmayacak." Dördü de çözüldü, koşu YAPILMADI.

### 1. Konservasyon: %100 nereden geliyordu

**Kök neden ikiydi.** (a) `conservation.run` `agree`/`total` sayaçlarını 1'den
başlatıyordu — aday kendi kendisiyle her pozisyonda aynı. HİÇ ortolog
bulunamadığında profil baştan sona 1.0 ve `conservation_percent` tam 100.0
çıkıyordu: **ölçüm yokluğu "mükemmel korunmuş" diye raporlanıyordu.**
(b) Otomatik suş çekme virüslerde çalışmıyordu: `genome/taxon/.../dataset_report`
(assembly uç noktası) BoHV-1 için **tek** kayıt döndürüyor (`GCF_008777455.1`,
referansın kendisi) — ölçüldü. Yani "3 suş çek" isteği 1 suş getiriyordu.

**Düzeltme:**
- Konservasyon artık ORTOLOGLARIN oranı; aday kendisi paya/paydaya girmiyor.
- `MIN_ORTHOLOGS = 2`: tek ortolog "suşlar arası korunmuşluk" değil, iki dizi
  arası kimliktir → altında kalırsa değer **None** ve yöntem etiketi
  "ölçülemedi — n ortolog bulundu, en az 2 gerekli".
- Aynı genomu iki kez sayma koruması: bir suşun bulunan ortologlarının tamamı
  (>=3) adayla birebir aynıysa o suş kıyaslamadan çıkarılır ve raporda geçer.
- `strain_fetch`: VİRÜS DALI eklendi. `virus/taxon/{taxid}/dataset_report`
  (BoHV-1: 854 kayıt, 71'i tam proteom) + E-utilities `fasta_cds_aa` ile protein
  indirme; parça kayıtlar (protein_count < 50) dışlanır. Datasets'in virüs
  indirme ucu protein FASTA vermiyor (denendi, zip yalnız data_report içeriyor).
- `ingest.accession_from_file`: girdinin kendi accession'ı (BoHV-1: `JX898220`)
  çekimden DIŞLANIR. Suş sayısı 3 → 5.

**Doğrulandı (gerçek NCBI + gerçek diamond, 25 aday protein, 5 izolat:**
BoHV-1.2b/CHN/SZ01/2022, PVNRTGVU/VBT/2023/0761, BC01, BJ202406,
**BHV-1/TR 1305**): konservasyon %83.4-100.0 aralığı, ortalama %97.8, tam %100
çıkan yalnız 3/25. Yani bileşen artık gerçekten ayırt ediyor.

### 2. Konservasyon profili KAYIYORDU (testle bulunan yeni hata)

Hizalayıcı `local` modda; `aln[0]` yalnız HİZALANAN parçayı içeriyor, ama eski
kod kolonları sayarken sayaca 0'dan başlıyordu. Hizalama adayın 1. kalıntısından
başlamadığı her durumda **profil kayıyordu**: `MKVLAAG…` ↔ `AKVLAAG…` örneğinde
gerçek uyuşmazlık (poz. 0) "korunmuş", gerçek eşleşme (son poz.) "korunmamış"
okunuyordu. Profil `epitope_conservation`'ı besliyor, o da adaylık puanına
giriyor — kayma doğrudan epitop sıralamasını bozuyordu.

**Ölçülen yaygınlık:** 74 gerçek hizalamanın **12'si (%16)** adayın 0.
kalıntısından başlamıyor; kayma **23-100 kalıntı** (ortalama 51). O proteinlerde
epitop-penceresi konservasyonu tamamen alakasız pozisyonları okuyordu.
**Düzeltme:** kolon sayımı bırakıldı, Biopython hizalama koordinatları
(`aln.indices`, yedek `aln.aligned`) kullanılıyor (`_aligned_pairs`).

### 3. MEV'in kendi alerjenitesi artık eyleme dönüşür uyarı

Adaylarda `require_non_allergen` SERT filtre, yani her epitop temiz; buna karşın
birleştirilmiş konstrükt FAO/WHO 6-mer ölçütüne göre alerjen çıkıyordu ve rapor
bunu tek kelimeye ("Allergen") indiriyordu. `allergen.matches()` artık TÜM
eşleşmeleri konumuyla veriyor, `allergen._kmer_sources()` hangi alerjene
benzediğini söylüyor (142.320 6-mer → kaynak eşlemesi), `mev.allergen_attribution`
her eşleşmeyi bölgeye atfediyor:

  * `junction` — epitop/linker birleşme yerinde doğan YENİ 6-mer → **linker
    değiştirilerek giderilebilir**, antijenlerin özelliği değil.
  * `linker` — linker dizisinin kendisi (aynı çözüm).
  * `adjuvant` — β-defensin konak kökenli, eşleşme BEKLENİR (not, eleme değil).
  * `epitope` — epitobun içinde → aday filtresi kaçırmış, o epitop çıkarılmalı.

Rapora uyarı kutusu + tablo (6-mer / konum / köken / parçalar / eşleşen alerjen)
ve karar cümlesi geldi. **Test sırasında planlanmamış bir bulgu:** standart
`EAAAK` adjuvan linker'ı ardından E ile başlayan bir epitop gelince `EAAAKE`
alerjen 6-mer'i doğuyor (kaynak Q8J1L4, *Trichophyton schoenleinii*).

### 4. pLDDT güveni alt sonuçlara YAYILIYOR

MEV pLDDT ortalaması 55.6, eşik 70 — altında. Rapor bu sayıyı diğer satırlar
gibi yazıyordu, oysa yapıdan TÜREYEN her şey (RSA/yüzey yüzdesi, konformasyonel
B-hücre yamaları) o modelin güvenine bağlı. `structure.confidence()` +
`structure.region_plddt()` eklendi: rapora uyarı banner'ı, yüzey/konformasyonel
satırlara "(düşük güven)" etiketi, her B yamasına KENDİ ortalama pLDDT'si ve
eşik altındaysa ⚠. Banner ayrıca MEV'in linkerlarla birleşmiş yapay bir kimera
olduğunu, yani düşük pLDDT'nin BEKLENEN olduğunu — kusur olmadığını — yazıyor.

### 5. `tool_status` yük altında yanlış "ARAÇ YOK" demiyor

`verify()` aracı gerçekten koşturduğu için geçici bir kaynak hatasıyla da
düşebiliyordu (ölçülmüştü: 9 paralel netMHCIIpan + pytest altında ToxinPred2).
Artık bir kez RETRY var (`_verify_with_retry`) ve eksikliğin TÜRÜ ayrı:
`missing` (dosya yok, kalıcı) / `verify_failed` (dosya var, doğrulama koşusu
düştü — geçici olabilir) / `probe_error`. Rapor üçünü ayrı yazıyor; "dürüstlük
tek yönlü değil, 'var' demek kadar 'yok' demek de yanlış olabilir."

**Testler:** `tests/test_biyolojik_olcutler.py` (14 test) + tüm paket **75/75**.
**Koşu YAPILMADI.**

## 2026-11-22 (üçüncü tur) — poly-Ala epitopları: iki ölçek hatası + eksik ölçüt

**Şikâyet:** MHC-I top epitop listesinde `AAAAARAA` gibi düşük karmaşıklıklı
peptitler var; "skorlayıcı Ala tekrarlarını ödüllendiriyor olabilir" diye not
düşülmüştü. **Bu hipotez ÖLÇÜLDÜ ve YANLIŞ çıktı** — ama yerine üç gerçek
kusur bulundu.

### Ölçüm (gerçek netMHCpan, sığır BoLA 12 allel, 1500 doğal BoHV-1 9-mer'i)

| peptit | %rank | doğal peptitler içindeki yeri |
|---|---|---|
| `AAAAARAA` | 1.131 | ilk %19.7 |
| `AAAAAAAAA` | 1.240 | ilk %21.3 |
| `GGGGGGGGG` | 47.9 | ilk %93 |
| doğal %1'lik dilim | 0.037 | — |

netMHCpan poly-Ala'yı ödüllendirmiyor; orta sıralarda bırakıyor. Düşük
karmaşıklıklı doğal peptitlerin ortalama rank'i (12.58) yüksek karmaşıklıklı
olanlardan (14.64) anlamlı farklı değil. Yani **skorlayıcı suçsuz**.

### Bulunan üç kusur

1. **`mhc_score = 1 - rank/100` ayırt etme gücünü yok ediyordu.** %rank bir
   yüzdelik dilimdir, biyolojik aralığı 0.01-2'dir; bu doğrusal dönüşüm o
   aralığın tamamını 0.98-1.00'e sıkıştırıyordu. Ölçüldü: eşiği geçen
   peptitlerde skorun **toplam yayılımı 0.0196** — `mhc_score` pratikte sabit.
   `scoring.py` bu bileşeni ağırlıklı ortalamaya kattığı için %rank 0.003'lük
   gerçek bir bağlayıcı ile %rank 1.13'lük sınırda bir peptit adaylık puanında
   ayırt edilemiyordu (0.99997 vs 0.9887).
2. **IEDB omurgası bambaşka bir ölçekte skor veriyordu.** `iedb.predict` rank'ı
   hiç kullanmayıp API'nin kendi `score` alanını döndürüyordu. Yerel araç
   düşüp IEDB devralırsa adaylık puanları sessizce kıyaslanamaz hale geliyordu
   (araç dürüstlüğü notundaki "ÖLÇEK UYUŞMAZLIĞI" varyantının üçüncü örneği).
3. **Hiçbir katman düşük karmaşıklığı ölçüt saymıyordu.** `AAAAARAA` %rank
   1.131 ile zayıf-bağlayıcı eşiğini (2.0) GEÇİYOR, dolayısıyla tabloya,
   ★ işaretine, MEV kasetine ve `top_peptides.fasta`'ya girmeye hak kazanıyordu.

### Düzeltme

- **`vaxforge/mhc_rank.py` (yeni, TEK kaynak):** %rank → skor dönüşümü log10
  ölçekte ve sınıfın KENDİ eşiklerine demirli (`rank <= güçlü/10` → 1.0,
  `rank >= zayıf` → 0.0). Üç omurga da (yerel netMHCpan, IEDB, vekil) bunu
  kullanıyor, `rank_strong` `epitope.run`'dan config'ten geçiriliyor.
  Ölçülen etki: eşik geçen peptitlerde yayılım **0.0196 → 0.9948**;
  `AAAAARAA` 0.9887 → **0.1545**, en iyi doğal bağlayıcı 1.0'da kalıyor.
  Adaylık puanı: poly-Ala 0.41, gerçek bağlayıcılar 0.78-0.79.
- **`vaxforge/complexity.py` (yeni):** dizi karmaşıklığı AYRI bir dışlama
  ölçütü (skor cezası değil). Kural: `maxrun>=4 veya farklı_aa<=3 veya
  entropi<1.5 bit`. **BoHV-1 proteomunda kalibre edildi:** 14/14 homopolimer
  sondayı eliyor, gerçek güçlü bağlayıcıların (%rank<=0.5, n=145) yalnız
  %2.1'ini kaybediyor. Daha agresif kurallar (uniq<=4: %4.8, H<2.0: %5.5)
  reddedildi. Ala'dan ZENGİN olmak eleme sebebi değil: `AGAPAYVAA`
  (%rank 0.008, en güçlü bağlayıcılardan) kuraldan GEÇİYOR — elenen tekrar,
  kompozisyon değil.
  *Biyolojik gerekçe:* homopolimer TCR'a ayırt edici yüzey sunmaz; düşük
  karmaşıklıklı bölgeler konak proteomunda yaygın (tolerans/otoimmünite);
  çoğunlukla düzensiz bölgelerden gelir; peptit sentezinde agregasyon yapar.
- **Dışlama GÖRÜNÜR, sessiz değil** (dürüstlük kuralı): epitop tablosunda
  satır kalıyor ama ★ alamıyor ve altında gerekçe yazıyor; MEV bölümünde
  "düşük karmaşıklık nedeniyle dışlananlar" tablosu; `candidates.csv`'de
  5 yeni sütun; `top_peptides.fasta`'da `>cand` değil `>excluded...
  EXCLUDED=low_complexity|reason=...` olarak dosyanın sonunda.
- `passed` alanına DOKUNULMADI — o bağlanma kararıdır, bu ayrı ölçüt.

**Testler:** `tests/test_complexity.py` (12 test) + tüm paket **59/59 geçiyor**.
**Koşu yapılmadı** (Ali'nin emri: biyolojik sorunlar çözülmeden ve o demeden
yeni koşu yok). Eski `outputs/run_*` klasörleri Ali'nin isteğiyle silindi.

## 2026-11-22 (ikinci tur) — MHC-II sessiz düşüşü, tehlike sütunu polaritesi, plan donması

BoHV-1 raporu incelenirken üç hata bulundu ve düzeltildi:

1. **MHC-II hiç koşmamış, ama rapor "21/23 hazır" diyordu.** `tool_status`
   netMHCIIpan'ı `available` sayıyordu (dosya var) ve araç gerçekten çalışıyordu;
   düşme noktası `netmhc_local.predict`'teki sabit `CHUNK = 20000` idi. Maliyet
   peptit DEĞİL **peptit×allel** ile büyüyor: ölçüldü, bu makinede ~115
   (peptit×allel)/sn. 35.993 15-mer × 11 BoLA alleli = 2 chunk × ~30 dk, yani
   1200 sn timeout'un iki katı → iki chunk da düştü, `except Exception: return ""`
   tek kelime etmeden yuttu, `agg` boş kaldı, `predict` None döndü ve çağıran
   heuristik vekile indi. 407 adayın 262'si (tüm MHC-II) proxy skoruyla
   sıralanmıştı; gerçek araç 1 numaralı adaya %rank 33.6 (bağlanmaz) diyor,
   proxy 1.17 verip listenin başına koymuştu.
   **Fix:** chunk iş miktarına göre bölünüyor (`_UNITS_PER_CHUNK=20.000`;
   11 allelde chunk 1818), düşen chunk `last_diagnostic()`'e yazılıyor, kısmi
   sonuç yöntem etiketinde `· KISMİ n/m peptit` olarak görünüyor.
   **Doğrulandı:** tam havuz 13.6 dk, 35.993/35.993 peptit, 20/20 chunk, hiç
   düşen yok, rank≤10 bağlanan 12.028.
2. **Alerjen/Toksik sütunlarının polaritesi tersti.** Başlık tehlikenin adını
   taşıyor, hücre ise `allergen_ok`/`toxic_ok` basıyordu: "Alerjen" sütununda ✔
   aslında "alerjen DEĞİL" demekti. Üstelik toksisitedeki ✗'lerin çoğu gerçek bir
   iddia değildi — ToxinPred2 yalnız ≥35 aa'da doğrulanmış, epitoplar 8-15 aa, bu
   yüzden eleme bilinçli yapılmıyor; koşuda 407 adayın 348'i eşiğin üstündeydi,
   yani sütun hiçbir şeyi ayırt etmiyordu. **Fix:** başlıklar `Alerjenite` /
   `Toksisite`, hücre durumun kendisini yazıyor (`Temiz` / `ALERJEN` / `TOKSİK` /
   `? alan dışı (skor)`), dört katmanda birden (HTML, PDF, Streamlit, Next.js).
3. **Plan tablosunun 8 adımı da `pending` donuyordu.** `meta["plan"]` koşunun en
   başında bir kez alınıyordu. **Fix:** `pipeline.run` içinde yerel `_ev` her
   olayı Step'e işliyor, plan tazeleniyor, adımın son mesajı `not` sütununa
   geçiyor.

**Bu turda bulunan ama DOKUNULMAYAN açık:** `tool_status` yük altında yanlış
"ARAÇ YOK" diyor. 9 paralel netMHCIIpan + pytest koşarken ToxinPred2 `verify()`
düştü ve araç eksik işaretlendi; tek başına `verify() -> True`. Yani geçici bir
yük hatası rapora "ARAÇ YOK" yazdırabiliyor — ilk düzeltilen hatanın tersi ama
aynı derecede yanıltıcı. Düzeltmesi: bir kez retry + "dosya yok" ile "verify
düştü" ayrımı.

**Yapılmadı:** BoHV-1 gerçek MHC-II skorlarıyla yeniden koşulmadı; README
benchmark sayıları hâlâ eski koşudan.

## Tek bakışta

| Katman | Durum |
|---|---|
| Pipeline (`vaxforge/`, 50 modül) | ✅ çalışıyor, ROADMAP 11/11 faz kapalı |
| Streamlit arayüzü (`app.py`) | ✅ |
| Web yığını (`web/` FastAPI + Next.js 14) | ✅ bu makinede 2026-11-17'de ilk kez koştu (node elle kuruldu) |
| Benchmark (10 bakteri + 9 virüs, 3'lü head-to-head) | ✅ koşuldu, atıf tamam |
| Makale (`benchmark/MANUSCRIPT*.md`) | 🔨 taslak hazır, gönderim öncesi rötuş |
| CI | 🔨 iş akışı yazıldı, HENÜZ KOŞULMADI |

Benchmark başlık metriği (10 bakteri, 48 koruyucu antijen, DEFAULT parametre):
VaxForge 37/48 · Vaxign-ML 47/48 · NERVE 2.0 14/48. MenB exemplar'ında VaxForge 4/4,
NERVE 1/4. Kaçan 11'in tamamı `benchmark/MISS_ATTRIBUTION.md`'de atfedilmiş; ölçekleme
hatası YOK.

## 2026-11-17 — eksik araçlar kapatıldı (11/22 → 21/23), arayüz ilk kez koştu

Önceki turların aksine bu turda **her araç gerçek girdiyle koşturularak** doğrulandı
(`available()` değil: boş olmayan + farklılaşmış çıktı aranır — bkz. 2026-11-13 dersi).

**Arayüz bu makinede hiç çalışmamıştı:** `node`/`npm` kurulu değildi, venv'de
`fastapi`/`uvicorn` yoktu. Node 20.20.2 `~/.local/opt/node`'a (sudo'suz) kuruldu,
backend bağımlılıkları eklendi, `web/guardian/supervise.sh` ile ayağa kalktı
(`:3000` + `:8011`). Next.js derlemesi Google Fonts'u ağdan çekiyor; ilk denemede
düşerse tekrar koşmak yeterli.

**ImmForge `ModuleNotFoundError` — kod değil dal eksikti.** `immforge/integrations/
vaxforge.py` yalnız `phase-4.13b-anchors` dalında; yerel ImmForge shallow + yalnız
`main` klonlanmıştı. `~/repos/ImmForge-phase4` worktree'si açıldı, `IMMFORGE_PATH`
oraya bakıyor (`supervise.sh` artık varsa kendisi ayarlıyor). `calibration.bovine.yaml`
var → sığır human-proxy DEĞİL. Doğrulandı: import + `evaluate_candidate` imzası uyuyor.

**IEDB taraması neden hiç koşmadı:** `organism_taxon` boştu → aday-başına canlı moda
düşüldü → 703 > `LIVE_MAX`=150 → dürüstçe atlandı. Oysa girdi `.gb` dosyasıydı ve
taxon'u `source` feature'ının `db_xref`'inde taşıyor. `ingest.taxon_from_file()`
eklendi, `pipeline.run` kullanıcı vermediyse dosyadan okuyor (`meta.organism_taxon`
kaynağı da yazılır). Doğrulama: `br1.gb` → `NCBITaxon:430066`, `ebola.gb` →
`NCBITaxon:186538`, FASTA → `None`.

BoHV-1 koşusunun eksik kalan sonucu sonradan hesaplandı: taxon 10320 ile **58 deneysel
doğrulanmış lineer epitop** indi, 703 adaydan **7'si eşleşti** (6 validated), recall
4/57 = **%7.0**. Eşleşmeler **US6 = glikoprotein D**'de (BoHV-1'in bilinen koruyucu
antijeni) → gerçek pozitif kontrol. Uyarı: bilinen setin 34/58'i tek referanstan
(`top_ref_share` 0.59), yani payda taraflı.

**VFDB zinciri → UniProt KW-0843 (kullanıcı kararı).** `www.mgc.ac.cn` bu ağdan DNS'te
çözülmüyor; VFDB olmayınca virülans eşleşmesi + ACC modeli + VaxiJen birlikte düşüyordu.
2.992 reviewed bakteriyel virülans proteini (`keyword:KW-0843 AND reviewed AND
taxonomy_id:2`) indirildi, DIAMOND DB kuruldu. `discovery.virulence_source()` tek
doğruluk kaynağı: VFDB varsa o, yoksa UniProt. **Kaynak adı hiçbir yerde sabit
yazılmıyor** — yöntem etiketi, `tool_status` etiketi, `citations` tablosu, `runner.py`
adım grid'i, `app.py` huni tablosu ve `pipeline-diagram.tsx` hepsi koşan kaynaktan
okuyor. Atıf tablosunda yalnız koşan `used=True` (VFDB `False` / UniProt `True`).

Doğrulama (Brucella abortus S19, 1967 protein): **88 DB eşleşmesi**, 103 benzersiz
kimlik değeri, virülans skoru 90 benzersiz değer. Eşleşenler gerçek virülans faktörleri:
BspG, üreaz operonu (ureA-G), Hfq, TcpB, BtaE, VtlR.

**ACC/VaxiJen'de İKİNCİ kusur: negatif set yolu da yanlıştı.** Kod
`tools/db/human_sprot.fasta.gz` arıyordu; diskteki dosya `human_sprot_auto.fasta`
(sıkıştırılmamış, farklı ad). Yani **VFDB indirilmiş olsa bile model asla
eğitilemezdi** — eksiklik tek başına VFDB'ye bağlı değildi. `antigen_acc` artık her iki
seti aday listesinden çözüyor, `.gz` olmayan dosyayı da açıyor, ve hangi pozitif setle
eğitildiğini modelin İÇİNE yazıyor (`model_info()`; rapor etiketi bunu okur).

**ACC modelinin gerçek gücü ölçüldü (ayrı tutulan 300 protein):** dengeli doğruluk
**0.747** — virülans ort. 0.686, insan ort. 0.254, insanda %24 yanlış pozitif.
Çıktı dereceli (%12 doymuş, 150'de 113 benzersiz). Ayrım var ama mütevazı: **ikincil,
zayıf kanıt hattı**. Birincil antijenite IApred kalmalı. (Uyarı: pozitif=bakteri,
negatif=insan olduğu için model "antijenik mi" yerine "bakteriyel mi" ayırıyor
olabilir — bu VaxiJen'in kendi yayınlanmış tasarımının da sınırı.)

**NetCTL-1.2 → IEDB processing API (kullanıcı kararı).** Binary bu çekirdekte satıcının
kendi test dosyasıyla bile segfault veriyor; yerelde çözümü yok, `processing_norm` her
adayda `None` kalıyordu. `iedb_processing.py` yazıldı (`tools_api/processing/`,
`proteasome_score` + `tap_score`, önbellekli, CHUNK=12) + `processing.py` omurga
seçici: NetCTL gerçekten koşuyorsa o, yoksa IEDB. Yöntem etiketi "NetCTL-1.2" DEMEZ.

**ÖLÇEK TUZAĞI — ve kendi iddiamın yanlışlanması.** IEDB'nin `proteasome_score`'u
NetCTL'in `cle`'si değil: 0-1 olasılık değil, log-ölçekli. İlk yorumuma "netctl'in
normalizasyonu her şeyi 1.0'a doyurur" yazdım; **ölçünce yanlış çıktı** (154/157
benzersiz, 0 doymuş). Gerçek etki ölçüldü: 105 cle değerinin **%68'i 1.0'ı aşıyor**
(aralık 0.321–1.722), `netctl.processing_norm` bunları kırpıyor → kesim bileşeni tek
başına 35/105 benzersiz seviyeye düşüyor, sigmoidle 102/105. Nihai skor doymuyor çünkü
TAP terimi kırpılmıyor. Yani kusur "her şey 1.0 olur" değil, **"kesim çözünürlüğünün
üçte ikisi kaybolur"**. Yorum ölçülen sayılarla düzeltildi.

Doğrulama: 3/3 protein döndü, pozisyonlar sınır içinde (61 aa dizide indeks 8–60),
153/157 benzersiz normalize skor.

**Kurulan ve gerçek koşuyla doğrulanan diğer araçlar:**
- **BepiPred-3.0** — `tools/bepipred-3.0` HİÇ yoktu (yalnız 1.0 vardı). `bp3` 0.0.12.7
  (Fold1-5 gömülü) + resmi `bepipred3_CLI.py`. Doğrulama: 145 aa diziye 145 skor,
  145 benzersiz değer, 0.034–0.314, eşik 0.1512'nin iki yanında → ayrım var.
  **Mevcut koşudaki 324 B-hücre adayı (703'ün %46'sı) vekille üretilmişti.**
- **metapredict 3.1.1** — Cython eklentisi derleyici ister; makinede `gcc` DE,
  Python başlıkları DA yoktu. İKİ sistem paketi şart: `build-essential` +
  `python3.12-dev`. Doğrulama: yapısal protein %0 düzensiz (ort 0.112), düşük
  karmaşıklıklı dizi %100 (ort 0.981).
- **IEDB Population Coverage 3.0.2** — doğrulama: 2 bilinen epitop → World %52.43.
  NOT: frekans verisi yalnız insan HLA'sı; sığır/fare konakta dürüstçe boş döner.
- **Protein-Sol** — doğrulama: `population_sol` 0.446 = yayınlanan referans ortalaması
  0.45 → `seq_reference_data.txt` gerçekten okunuyor.
- **S4PRED** — ağırlıklar MD5 `e04ad7d10b61551f7e07a86b65bb88dc` (yayınlanan değerle
  eşleşti). Çalıştırma doğrulaması YAPILMADI.
- **DeepLoc-2.1 — KAPATILDI (kullanıcı kararı).** CLI kuruldu
  (`pip install ./tools/deeploc2_package`), tek eksik ESM1b. Ölçüldü: ESM1b checkpoint
  **7.29 GB** (ESM-2 650M aynı parametre sayısıyla 2.43 GB — ESM1b optimizer state
  taşıyor), `torch.load` dosyayı bütün okuyor, makinede 7.6 GB RAM. Bu makinede
  koşmaz; GPU/servis makinesine ait. Virüs profilinde lokalizasyon zaten yumuşak filtre.

`requirements.txt` kurulum notlarıyla güncellendi (hangi paket, eksikse hangi vekile
düşüyor, elle indirilenlerin adresleri, S4PRED MD5'i, gereken apt paketleri).

**Kalan eksikler (2/23):** DeepLoc-2.1 (karar), NetCTL-1.2 (segfault — yeri IEDB
processing ile dolduruldu, işlev kaybı yok).

**AÇIK KARAR — örtüşen aday kümelenmesi.** 703 adayın **484'ü** aynı proteinden ≥4 aa
kayan-pencere komşusu (birebir duplikasyon 0; iç içe 14; tandem yineleme 561;
homopolimer ≥3 134; düşük karmaşıklık 262). Üst sıradaki rank 4/5/8/9 (UL44) tam bu
düşük-karmaşıklıklı prolin-serin tekrarları. Kodda kümeleme adımı YOK → **aday sayısı
şişkin**. Eklemek makaledeki sayıyı değiştirir, bu yüzden yapılmadı.

**Yeniden koşu gerekiyor:** mevcut `outputs/run_2026-11-16T21-59-16` B-hücre
epitoplarını vekille üretti, IEDB/ImmForge/çözünürlük/düzensizlik/ikincil yapı
adımları boştu. Aynı dosyayı arayüzden tekrar yüklemek yeterli (taxon artık otomatik).

Testler: 40/40 (üç ayrı koşu). Not: `pytest` çıktısını `tail`'e borulamak log dosyasını
boş bırakıyor ve koşu sessizce ölmüş gibi görünüyor — doğrudan dosyaya yazdır.

## 2026-11-12 — kod gözden geçirmesi, acil maddeler kapatıldı

Koşu yapılmadı (kullanıcı talimatı); tüm değişiklikler statik olarak doğrulandı
(`py_compile`, `tomllib`, `yaml.safe_load`). **Henüz commit edilmedi.**

1. **ImmForge yolu sessiz kırıktı** — `immforge_score.py` varsayılan olarak yalnız
   `~/immforge` arıyordu, repo ise `~/repos/ImmForge` altında. `IMMFORGE_PATH`
   verilmeyen her koşuda `calibrated_species()` BOŞ dönüyor, yani dokuz türün
   kalibrasyonu sessizce insan-proxy'sine düşüyordu. Artık aday konum listesi taranıyor.
   Doğrulandı: 9 tür bulunuyor (bovine, fish, goat, horse, human, mouse, pig, rabbit, sheep).
2. **Sessiz hata yutma** — yeni `vaxforge/diagnostics.py`: `setup()` logger kurar,
   `record()` yakalanan istisnayı `meta["diagnostics"]`'e yapısal yazar (aşama, tip,
   mesaj, son traceback karesi) ve loglar. Pipeline'ın altı `except` bloğu bağlandı,
   `meta["diagnostics_summary"]` rapor öncesi yazılıyor, yani run.json'a giriyor.
   Artık "araç yok" ile "araç çöktü" ayırt edilebiliyor.
3. **CI yok** — `.github/workflows/ci.yml`: birim testleri (ağsız, hafif bağımlılık),
   ruff, frontend `tsc`. Smoke testi `pytest.mark.smoke` ile işaretlendi ve CI'dan
   `-m "not smoke"` ile hariç tutuldu (yerel araçlar + ~4 dk).
4. **Lint/paketleme yok** — `pyproject.toml`: proje metadatası, ruff (kademeli kural
   seti), pytest marker'ları. Sürüm sabitleme BİLİNÇLİ boş; uydurma pin yazılmadı,
   gerçek ortamdan üretilmeli: `pip freeze > requirements.lock.txt`.
5. **README gerçeği yansıtmıyordu** — tavuk konağı seçilebilir gösteriliyordu ama
   `hosts.yaml`'da yorum satırında, yani kapalı (gerçek: 4 konak, domuz MHC-II dürüstçe
   devre dışı). Web yığınından hiç söz edilmiyordu. Benchmark sonuçları yoktu. Modül
   tablosunda on üç gerçek adım eksikti (PSORTb, ESMFold/RSA/konformasyonel B,
   glikozilasyon, klonlama, VaxiJen, mimikri, çok-suş konservasyon, konak-DB, IFN-γ,
   popülasyon kapsamı, ImmForge kararı, HLA-hastalık). TR ve EN birlikte düzeltildi.
6. **CORS joker** — `web/backend/main.py` `allow_origins=["*"]` sabitti. Artık varsayılan
   yalnız localhost:3000/3001; LAN için `VAXFORGE_CORS_ORIGINS` ile açıkça verilir.
   Dosya servisi tarafı zaten güvenliydi (beyaz liste + `_safe_run_dir`), dokunulmadı.
7. **Üretilen Word dosyaları izleniyordu** — `benchmark/*.docx` gitignore'a alındı ve
   `git rm --cached` ile izlemeden çıkarıldı (dosyalar diskte duruyor, `md_to_docx.py`
   yeniden üretir).

## 2026-11-12 (2) — bilimsel doğruluk düzeltmeleri

Öncelik "çalışması + biyolojik olarak doğru olması" olarak belirlendi; koşu/test
sonraya bırakıldı. Aşağıdakiler **sonucu değiştirir**, bu yüzden benchmark sayıları
yeniden koşulana kadar ESKİ koşuya aittir (README'ye bu caveat yazıldı).

8. **Gram-pozitif lokalizasyon tutarsızlığı (en önemlisi).** Bakteri izinli
   lokalizasyon listesi tek ve Gram-negatif mantığıyla yazılmıştı. Gram-pozitifte dış
   zar yoktur, sitoplazmik zar en dış sınırdır ve yüzey proteinleri PSORTb'de
   `membrane` alır; liste bu yüzden meşru yüzey antijenlerini kesiyordu (benchmark'ta
   iki kayıp doğrudan buradan). Artık izinli liste **Gram tipine göre** seçiliyor:
   - Gram-negatif: extracellular, outer_membrane, cell_wall, periplasm (değişmedi)
   - Gram-pozitif: extracellular, cell_wall, membrane (YENİ)
   `config_loader.resolve()` artık `gram` alıyor ve `<profil>_gram_<tip>` anahtarını
   tercih ediyor. `pipeline.run` ve Streamlit eşik tablosu gram'ı geçiriyor, yani
   panelde görünen liste koşulanla aynı. Varyant anahtarı tanımlı değilse eski
   davranış aynen korunuyor.
9. **`hard_filter` bayrağı kodda hiç okunmuyordu.** `thresholds.yaml` her araç için
   sert/yumuşak ilan ediyor ve README "her şey config'te" diyordu, ama funnel
   lokalizasyon ile konak-homolojisini koşulsuz sert uyguluyordu; YAML'da
   `hard_filter: false` yazmak hiçbir şeyi değiştirmiyordu. Artık ikisi de bayrağa
   uyuyor. **Varsayılanlar true olduğu için bugünkü davranış birebir aynı.** Yan
   kazanç: bekleyen self-homoloji yumuşak-mod denemesi artık KOD DEĞİŞİKLİĞİ
   GEREKTİRMEZ, tek satır YAML ile açılır.
10. **CAI atıf verilen yöntemle uyuşmuyordu.** Sharp & Li 1987 tek-kodonlu aileleri
    (Met/Trp) CAI'den hariç tutar; kod onları dahil ediyordu. Daha önemlisi
    `cai(back_translate(x))` tanım gereği 1.0'dır (her kalıntıya en sık kodon
    veriliyor), yani raporlanan CAI bir totolojiydi ve "CAI>0.8 iyi ekspresyon"
    notu boştu. Eklenen `expected_cai()` aynı proteinin **doğal E. coli kodon
    dağılımıyla** beklenen CAI'sini deterministik olarak hesaplıyor; 1.0'a mesafe
    optimizasyonun gerçek kazancı. Rapor alanları: `cai`, `cai_baseline`,
    `n_skipped_residues` (ters-translasyonda sessizce atılan tanınmayan kalıntı).
11. **Örtüşen N-glikozilasyon sequon'ları sayılmıyordu.** `N[^P][ST]` deseni
    `finditer` ile örtüşmeyen eşleşme döndürüyordu; "NNSS"te iki sequon varken biri
    bulunuyordu. Bitişik sequon gerçek glikoproteinlerde yaygın ve eksik sayım
    glikan-maskeli B-epitoplarını olduğundan az gösteriyordu. İleri-bakış deseniyle
    düzeltildi; mevcut birim test kriterleri korunuyor.

Denetlenip **doğru bulunanlar** (değişiklik yapılmadı): Tien 2013 maksimum ASA
tablosu, mimikri BLOSUM62 konservatif-ikame eşiği ve matris yokken dürüst atlama,
PSORTb etiket eşlemesi, PSORTb 'Unknown' durumunda ihtiyatlı geçirme, dosya servisi
yol güvenliği.

## 2026-11-12 (3) — temizlik denetimi + kapı doğrulaması

**Temizlik yapılmadı, çünkü silinecek bir şey YOK.** Denetim sonucu:
- `outputs/` boş; sistemde hiçbir `run_*` dizini yok, yani eski koşu çıktısı sıfır.
- `benchmark/results` 1.5 MB ve 72 dosyanın 72'si git'te izleniyor, yani tamamı
  küratörlü makale verisi. Yeniden üretilebilir ham rakip çıktıları (nerve/,
  vaxignml/) daha önce temizlenmiş.
- Disk 953 GB boş, repo 34 MB, `~/repos` toplam 159 MB. Bellek 7.5 GB'ın 6.3'ü boş.
- **Sonuç: darboğaz disk değil.** Gerçek kısıt bellek ve CPU süresi.

**Gram düzeltmesi GERÇEK koşuyla doğrulandı** (ağsız, araçsız, 5 yeni regresyon
testi, 5/5 geçti). `tests/test_modules.py` sonuna eklendi:
1. İzinli lokalizasyon listesi Gram-duyarlı; Gram-pozitifte `membrane` var, dış zar
   ve periplazma yok; Gram-negatifte `membrane` YOK; gram bilinmiyorsa eski liste.
2. Kayıtlı GERÇEK PSORTb çıktısı (`Cytoplasmic Membrane` → `membrane`, TM=1) artık
   Gram-pozitif kapıyı GEÇİYOR, Gram-negatif kapıda eleniyor. Tek engel
   lokalizasyondu: TM kapısı (max 1) bu antijenleri elemiyor.
3. `hard_filter` bayrakları hem config'te sert hem kodda okunuyor.
4. CAI tek-kodonlu aileleri hariç tutuyor, taban çizgisi 0 ile 1 arasında anlamlı.
5. Örtüşen N-sequon'lar sayılıyor, Pro kuralı korunuyor.

Koşucu (pytest kurulamadığı için asgari shim):
`/tmp/claude-1000/-home-ali/805871b6-9343-4b77-a857-fbd3975005e0/scratchpad/run_regression.py`

## 2026-11-12 (4) — iki dillilik kuralı uygulandı

Kural: arayüz/rapor metinleri **DAİMA TR + EN**. Önceki turda eklediğim uyarılar
yalnız Türkçeydi, düzeltildi.

- `i18n.py`: 11 yeni anahtar (hepsi TR+EN) + 4 yeni EN token eşlemesi. Toplam 152.
- `funnel.py`: yumuşak-mod izleri artık `reason` + `reason_en` taşıyor.
- `cloning.py`: `note` yanına `note_en`.
- `report.py`: klonlama bölümüne `cai_optimized` / `cai_baseline` /
  `cai_skipped` satırları + CAI totoloji uyarısı (her zaman) ve atlanan-kalıntı
  uyarısı (varsa). Yeni **tanılama bölümü** `_diagnostics_section`: beklenmeyen
  hata varsa uyarı tablosu, yoksa dürüst temiz notu. Yöntemler bölümüne
  Gram-duyarlı izinli-liste notu. Hepsi iki dilli, TR ve EN ayrı ayrı doğrulandı.
- İki yeni test: her i18n anahtarının TR+EN dolu olması, ve eklenen anahtarların
  gerçekten kullanılması (ölü anahtar yok). Tek dilli anahtar eklenirse test düşer.

Regresyon paketi: **7/7 geçti** (ağsız, araçsız).

## 2026-11-13 — ORTAM KURULDU, iki sessiz araç kusuru bulundu

Lisanslı araçların tamamı kullanıcının harici diskindeki `D:\İmmün Profil  Programları`
klasöründen çıktı (PC geçiş yedeğinde YOKtu). Kurulanlar: NetMHCpan 4.2,
NetMHCIIpan 4.3, SignalP 5.0b, TMHMM 2.0d, BepiPred 1.0c, NetCTL 1.2a, DeepLoc 2.1.
Ayrıca: venv + Biopython/pandas/pytest/streamlit/sklearn, DIAMOND 2.2.8, insan
Swiss-Prot DB (20440 protein), IApred, ToxinPred2 + RF_model, alerjen seti (1026),
IFNepitope2. Sistem paketleri kullanıcı tarafından kuruldu (python3.12-venv,
podman, tcsh, ncbi-blast+, unzip).

**12. DTU sarmalayıcıları yapılandırılmamıştı — ilk koşu GEÇERSİZDİ.**
`netMHCpan`, `netMHCIIpan` ve `netCTL` satıcı YER TUTUCU yoluyla gelir
(`setenv NMHOME /tools/src/...`). Düzeltilmeden araç hiç hesaplama yapmaz.
İlk kuduz koşusu bu yüzden epitopları **heuristik vekille** üretti; dürüstlük
etiketi (`heuristik proxy (netmhcpan yok)`) bunu yakaladı, ama ancak çıktı
incelendiğinde. Şüpheyi doğrulayan ipucu: insan/sığır/fare bağlanan-epitop
sayısı birebir 56'da eşit ve Jaccard örtüşmesi 1.0 (üç ayrı MHC paneli aynı
sonucu veremez; temel koşuda 31/33/26 ve J=0.88). Üç sarmalayıcıdaki yol
düzeltildi (orijinaller `*.orig` olarak saklandı).

**13. `runnable()` yapılandırılmamış aracı 'çalışır' sayıyordu.** Kontrol yalnız
çıktıda "format error" arıyordu; bozuk kurulumun gerçek mesajı ise "no binaries
found ...". Bu yüzden `runnable()=True` dönüyor, `predict()` sessizce vekile
düşüyordu — yani vekil sayılarla tamamlanan bir koşu GERÇEK sanılabiliyordu.
Eklenen `_wrapper_home_ok()` sarmalayıcıdaki NMHOME/CTLHOME yolunu okur ve
dizin yoksa aracı yok sayar; ayrıca "no binaries found" gibi mesajlar da
başarısızlık sayılır. Bu desene sahip olmayan araçlar (TMHMM/SignalP) hakkında
karar vermez.

**14. ImmForge entegrasyonu HİÇ YOK — immünizasyon adımı hiç çalışmamış.**
Yeni tanılama modülü ilk gerçek koşuda bunu yakaladı:
`ModuleNotFoundError: No module named 'immforge.integrations'`
(`immforge_score.py:133 in score_construct`). ImmForge deposunda `integrations`
alt paketi ne git geçmişinde ne PC-geçiş yedeğinde var; `immforge/bridge.py`
duruyor ama `evaluate_candidate` sağlamıyor. Yani pipeline'ın son adımı
(ters-aşı döngü kapanışı) hiçbir zaman karar üretmemiş. `available()` yalnız
DİZİN varlığına baktığı için "araç var" diyip çöküyordu. Artık
`entrypoint_available()` giriş noktasını gerçekten import etmeyi dener ve adım
gerekçeli, dürüst `ARAÇ_YOK` döner. **ImmForge tarafında entegrasyonu yazmak
ayrı bir iş ve karar KULLANICIYA ait — yapılmadı.**

Alel doğrulaması: `hosts.yaml`'daki 93 alelin tamamı kurulu sürümlerde tanınıyor
(kodun kendi `_norm_allele` çevirisi sonrası; 0 eksik). Domuz MHC-II dürüstçe
devre dışı.

Hâlâ eksik: PSORTb (rootless podman subuid aralığı 65536, imaj >1M UID istiyor →
`sudo usermod --add-subuids`), VFDB (mgc.ac.cn bu ortamdan çözülmüyor), IEDB
popülasyon kapsamı, DeepLoc ESM1b ağırlığı (7.3 GB, bakteri dalında gerekmez).

Testler: birim 21/21 yeşil; uçtan uca smoke dahil 21/21 (451 s) — smoke'tan
sonra eklenen 2 test hariç sayım. UYARI: IFNepitope2 modelleri sklearn 0.24 ile
turşulanmış, 1.5.2 ile açılıyor (`InconsistentVersionWarning`) → makale öncesi
`scripts/convert_ifn_models.py` koşulmalı.

Arayüz: `streamlit run app.py` → http://localhost:8501 (127.0.0.1'e bağlı).
Next.js arayüzü için node/npm gerekiyor, kurulu değil.

## Ortam: ESKİ NOT — sudo engelleri (artık aşıldı, PSORTb hariç)

Bu makinede Python ortamı çıplak: `pip` yok, `ensurepip` yok, `venv` oluşturulamıyor,
`tools/` dizini hiç yok. Kurulu tek üçüncü-parti paket pyyaml. Dış araçların hiçbiri
yok (DIAMOND, podman, PSORTb). Ağ ve disk sorunsuz, engel yalnız yetki.

Kullanıcının kendisinin koşması gereken tek komut:
```
sudo apt install -y python3.12-venv podman
```
Sonrasını ben yapabilirim: venv + hafif bağımlılıklar, DIAMOND ikilisi, UniProt
Swiss-Prot konak DB + diamond makedb, PSORTb imajının çekilmesi.

Lisanslı DTU araçları (NetMHCpan, NetMHCIIpan, SignalP-5.0, TMHMM-2.0, BepiPred-1.0)
akademik lisans formu ister; onları YALNIZCA kullanıcı indirebilir. Epitop aşaması
bunlara bağlı, ama **lokalizasyon kapısının doğrulanması için gerekmiyorlar**.

Tam karşılaştırma koşusunun maliyeti (kayıtlı süreler): B. anthracis 14392 s (~4.0 saat),
L. monocytogenes 11914 s (~3.3 saat) — tek bakteri tam proteom için saatler, 10 bakteri
için 35 saatten fazla CPU süresi.

## Sırada bekleyenler — İZİN GEREKTİRİR (koşu/bilimsel sonuç değişikliği)

Bunlar bilinçli olarak YAPILMADI, çünkü ya pipeline koşusu gerektiriyor ya da makaledeki
sayıyı değiştiriyor.

- **Benchmark'ın yeniden koşulması.** Gram düzeltmesi iki Gram-pozitif kaybı
  (CapB, ActA) kurtarma potansiyeli taşıyor; CAI ve sequon düzeltmeleri rapor
  alanlarını değiştiriyor. Yeni sayılar alınana kadar README ve makale ESKİ koşuyu
  gösteriyor (caveat yazıldı).
- **Self-homoloji yumuşak-mod denemesi.** Artık yalnız `thresholds.yaml` →
  `human_homology.hard_filter: false` ile açılıyor, kod değişikliği YOK. Altı kayıp
  bundan geliyor, ama bu makalenin "ihtiyatlı güvenlik kapısı" tezini değiştirir →
  bilimsel karar, senin.
- **Altın-standart regresyon testi.** Demo proteom üzerinde bir referans çıktı dondurup
  karşılaştırmak; pipeline'ı bir kez koşmak gerekiyor.
- **Bölüm kaydı refaktörü.** Her yeni modül şu an beş yere elle işleniyor (rapor,
  runner kataloğu, `results.tsx`, `lib/i18n.ts`, `app.py`). Tek bildirimsel kayıt bunu
  bire indirir; doğrulaması smoke koşusu ister.
- **Rapor/PDF birleştirmesi.** HTML tarafında 18 bölüm fonksiyonu var, PDF tarafında
  tek ~370 satırlık fonksiyon. İki çıktı zamanla ayrışıyor.
- **Sürüm kilidi** (`pip freeze`) ve **CI'ın ilk koşusu**.
- **`git fetch --unshallow`** — repo yüzeysel klon, yerelde tek commit görünüyor.

## Bilinmesi gerekenler

- **ImmForge bağımlılığı dışsaldır ve private.** Kalibrasyon YAML'ları bu repoda değil;
  hakem bu katmanı olduğu gibi üretemez. Gönderim öncesi ya kalibrasyon dosyaları
  buraya kopyalanmalı ya da ek malzeme olarak verilmeli.
- Lisanslı DTU araçları (NetMHCpan, NetMHCIIpan, SignalP, TMHMM, BepiPred-1.0) elle
  akademik indirme ister; yoksa pipeline dürüst fallback ile çalışır.
- `benchmark/` içindeki yedi kabuk betiği tek seferlik orkestrasyondur; reprodüksiyon
  yolu `checkpoint.py` + `PROGRESS.jsonl` üzerinden izlenir.

## 2026-11-13 (akşam) — GERÇEK ARAÇ KOŞUSU + 15. bulgu (domuz allel biçimi)

Gerçek araçlarla kuduz koşusu tamamlandı (2438 s, G protein rank #1/5, candidacy
0.7787). Konak ayrışması biyolojik olarak makul: Jaccard human~bovine 0.75,
human~mouse 0.67 (geçersiz vekil koşusunda hepsi 1.0 idi). `diagnostics_summary`
temiz: total 0, unexpected 0 → ImmForge adımı dürüst ARAÇ_YOK olarak kaydedildi,
çökme değil. 13. bulgunun düzeltmesi ilk gerçek koşuda doğrulandı.

### BULGU 15 — netMHCpan allel doğrulaması YANLIŞ KAYNAĞA yapılıyordu

Belirti: domuz 0 bağlanan epitop verdi (kayıtlı temel koşuda 19).

Kök neden iki katmanlı:

1. `_norm_allele(..., "mhc_i")` yalnız `*` karakterini atıyordu. Bu HLA için
   doğru (`HLA-A*01:01` → `HLA-A01:01`, çünkü `*` ardından zaten kolon var) ama
   SLA için yanlış: `SLA-1*0101` → `SLA-10101`, haneler kaynaşıyor. Doğru biçim
   `SLA-1:0101`.
2. `_supported()` doğrulamayı `-listMHC` çıktısına yapıyordu. O çıktı eş-anlamlı
   yazımları da basar (121 SLA satırı) ve `SLA-10101` İÇİNDE GÖRÜNÜR. Ama `-a`
   bayrağının gerçek kapısı `data/allelenames` (76 SLA kaydı) ve orada YOK.
   Sonuç: geçersiz allel süzgeçten geçiyor, araç `cannot be found in allelenames
   list` hatası veriyor, tüm chunk düşüyor, `agg` boş kalıyor.

ÖNEMLİ — dürüstlük katmanı burada DOĞRU çalıştı: `predict()` boş sonuçta `None`
döndü, çağıran vekile indi ve run.json'a `heuristik proxy (netmhcpan yok)`
yazdı; insan/sığır/fare için `GERÇEK (NetMHCpan yerel, …)`. Yani sistem "gerçek
araç koştu" yanılgısı üretmedi. Kayıp yetenek kaybıydı, yanlış iddia değil.

Bu aynı zamanda önceki allel doğrulamamın neden yanlış "93/93 geçerli" dediğini
açıklıyor: ben de `-listMHC` üst kümesine bakmıştım.

Düzeltme (`vaxforge/netmhc_local.py`):
- `_norm_allele`: `*` ardından kolon yoksa `*` → `:` çevrilir, varsa atılır.
- yeni `_tool_home()` + `_enforced_allele_names()`: sarmalayıcıdaki `NMHOME`
  satırından kurulum dizini bulunur, `data/allelenames` (MHC-I) ya da
  `data/allelelist.txt` (MHC-II) okunur. `_supported()` bunu CLI listesinden
  ÖNCE dener; dosya yoksa eski davranışa düşer (regresyon yok).

Doğrulama: 93 allelin tamamı zorlayıcı listeye karşı geçerli (önce MHC-I'de
48/59, 11 domuz alleli reddediliyordu). Domuz uçtan uca tahmin artık gerçek
ayrışmış sonuç veriyor (5 peptitte 2 bağlanan, en iyi rank 0.116).

Testler: 3 yeni regresyon testi (biçim, doğrulama kaynağı sırası, hosts.yaml'ın
tamamının araçça kabulü) → 24 birim testi geçiyor.

Domuz düzeltmesiyle koşu yeniden başlatıldı →
`benchmark/results/rabies_2026-11-13_domuz_fix/`. Önceki iki koşu kanıt olarak
duruyor (`rabies_2026-11-13/` vekil, `rabies_2026-11-13_gercek/` domuz eksik).

## 2026-11-13 (gece) — BULGU 16-19: "araç var sanılıyor" kusur sınıfının taraması

15. bulgudan sonra AYNI SINIFI sistematik taradım: her `available()`
fonksiyonunun GERÇEKTEN ne doğruladığını tek tek kontrol ettim. Dört kusur daha
çıktı. Hepsinin ortak deseni şu: araç kullanılamaz durumda ama `available()`
True diyor → metrik sessizce boş/NaN kalıyor → rapor bunu "olumsuz bulgu" gibi
gösteriyor, bazen de aracı atıfla birlikte "kullanıldı" diye listeliyor.

### BULGU 16 — BepiPred sarmalayıcısı yer tutucu yolunda, üstelik YEDEĞİ de kapatıyor

`tools/bepipred-1.0/bepipred` içinde `setenv BP /usr/cbs/packages/bepipred/1.0/...`
duruyordu (dizin yok). `bepipred.available()` yalnız sarmalayıcı dosyası + tcsh
varlığına baktığı için True dönüyordu. Araç `cannot find BP in ...` diyip
çıkıyor, `predict_residues()` BOŞ LİSTE dönüyor.

En kötü yanı: `available()` True olduğu için epitope.py'nin çalışan
Kolaskar-Tongaonkar+Parker YEDEĞİNE hiç düşmüyordu. Yani bozuk araç, çalışan
yöntemi devre dışı bıraktı. Gerçek koşuda sonuç: 85 adayın 60'ında
`bcell_score: NaN`.

Düzeltme: sarmalayıcıda `BP` ve `TMP` yolları gerçek kuruluma çevrildi (yedek
`bepipred.orig`), `tmp/` dizini açıldı. Yeni `_configured()` `setenv BP` yolunu
doğruluyor. Araç artık 37 aa için 37 gerçek artık skoru üretiyor.

### BULGU 17 — netCTL TMPDIR yer tutucu (`/scratch`)

Düzeltildi (`tools/netCTL-1.2a/tmp`). Ama asıl sorun 18. bulguydu.

### BULGU 18 — netCTL binary'si SEGFAULT veriyor, buna rağmen atıfla raporlanıyordu

CTLHOME doğru, TMPDIR doğru, sarmalayıcı doğru binary'yi çağırıyor — ve binary
`Segmentation fault (core dumped)` veriyor. Satıcının KENDİ `test/test.fsa`
dosyasıyla da çöküyor. Statik bağlı eski bir binary; yapılandırmayla düzelmez.

`netctl.available()` yalnız dosya + tcsh varlığına baktığı için True dönüyordu
ve `METHOD` sabiti **"GERÇEK (NetCTL-1.2 kesim+TAP, Larsen 2005)"** diyordu.
Sonuç: 85 adayın TAMAMINDA `processing_norm: None`, ama rapor Yöntemler/
Referanslar tablosunda NetCTL-1.2'yi Larsen 2005 atfıyla gösteriyordu.

Düzeltme: `available()` artık GERÇEK bir duman testi yapıyor (minik FASTA ile
çağır, dönüş kodunu ve ayrıştırılabilir satır var mı diye bak). Artık dürüstçe
False. `processing_norm` yine None ama bu sefer rapor da aracın YOK olduğunu
söylüyor.

### BULGU 19 — atıf tablosu KOŞULSUZ besleniyordu (kök neden)

16-18 tek tek kusurdu; bunların RAPORA yanlış yansımasının kök nedeni ayrıydı:
`citations.for_report()` 39 atıf satırını araç durumuna bakmadan döndürüyordu.
Yani kurulu olmayan her araç atfıyla birlikte kullanılmış gibi görünüyordu.
Ayrıca `population.py` araç yokken bile `method` alanına "GERÇEK (IEDB
Population Coverage, Bui 2006)" yazıyordu (yanındaki `available: false` notu
kurtarıyordu ama yöntem etiketini tarayan okuyucu yanılıyordu).

Düzeltme — yeni dürüstlük katmanı:
- YENİ `vaxforge/tool_status.py`: koşunun en başında 22 dış aracın gerçek
  `available()` durumunu ölçer; `available()` istisna atarsa YUTMAZ, False +
  istisna tipi olarak kaydeder. `meta["tool_status"]` + özet olarak run.json'a
  yazılır.
- `citations.for_report(tool_status)`: her satır üç değerli `used` alır
  (True hazır / False ARAÇ YOK / None eşlenmemiş saf-Python yöntem).
- `report.py`: YENİ `_tool_status_section()` iki dilli "Araç durumu" bölümü;
  eksik araçları adlandırır ve **eksik değerin olumsuz bulgu olmadığını** açıkça
  yazar. Referans listesinde eksik araç `[ARAÇ YOK]` / `[TOOL ABSENT]` ile
  işaretlenir.
- `population.py`: yeni `METHOD_UNAVAILABLE` sabiti; araç yokken "GERÇEK"
  etiketi basılmaz.
- `i18n.py`: 8 yeni TR+EN anahtarı (toplam 160, hiçbiri tek dilli değil).

Bu ortamda ölçülen durum: 22 araçtan 10'u hazır. Eksikler: PSORTb (podman
subuid), VFDB (www.mgc.ac.cn çözülmüyor), DeepLoc (7.3 GB ESM1b), NetCTL
(binary çöküyor), IEDB Population Coverage, BepiPred-3.0, ImmForge
entegrasyonu, VaxiJen/antijenite ACC, ikincil yapı, çözünürlük, metapredict.

Testler: 5 yeni regresyon testi → 29 birim testi geçiyor.

### Kuduz koşuları — dört sütunlu karşılaştırma

| ölçüt | temel (kayıtlı) | vekil (geçersiz) | gerçek araç | domuz düzeltmeli |
|---|---|---|---|---|
| süre (s) | 4472 | 262 | 2439 | 2386 |
| G rank | #2/5 | #1/5 | #1/5 | #1/5 |
| candidacy | 0.7102 | 0.8244 | 0.7787 | 0.8246 |
| insan | 31 | 56 | 54 | 55 |
| sığır | 33 | 56 | 53 | 53 |
| fare | 26 | 56 | 47 | 46 |
| domuz | 19 | 30 | **0** | **28** |
| J(insan~sığır) | 0.75 | 1.00 | 0.75 | 0.833 |
| J(insan~domuz) | 0.50 | 0.667 | 0.00 | 0.545 |

Vekil koşusunun imzası J=1.00 (üç konak birebir aynı — biyolojik olarak
imkânsız). Gerçek araç koşusu temel koşunun konak ayrışmasını yeniden üretti.
Domuz düzeltmesi sonrası domuz da temele yakın (28 vs 19; J 0.545 vs 0.50).
Dört koşunun hepsi kanıt olarak duruyor, hiçbiri silinmedi.

Son koşu (`rabies_2026-11-13_tam/`) BepiPred + netCTL + araç durumu
düzeltmeleriyle başlatıldı; `bcell_score` NaN sayısının düşmesi bekleniyor
(önceki koşularda 60/85).

### Kurulum tekrarlanabilirliği — `scripts/configure_dtu_wrappers.py`

Yaptığım sarmalayıcı düzeltmeleri `tools/` altında ve `tools/` .gitignore'da
(lisanslı araçlar repoya girmemeli — doğru karar). Yani bu düzeltmeler sürüm
kontrolünde SAKLANMIYOR; paketi yeniden çıkaran biri (ya da yeni makine) aynı
sessiz kusurla baştan karşılaşır.

Bu yüzden düzeltmeyi betikleştirdim: 4 sarmalayıcıdaki 6 `setenv` yolunu
(NMHOME ×2, CTLHOME, TMPDIR, BP, TMP) gerçek kuruluma çevirir. IDEMPOTENT
(doğru olanı tekrar yazmaz), ilk değişimde `<ad>.orig` yedeği bırakır ve
`--check` modunda hiçbir şeyi değiştirmeden sorun varsa çıkış 1 döner (kurulum
doğrulaması / CI için). Kopya üzerinde test edildi: yer tutucuyu yakalıyor,
düzeltiyor, ikinci çağrıda değişiklik yapmıyor.

Her iki README'ye kurulum adımı olarak eklendi (TR + EN).

Ayrıca `requirements.lock.txt` üretildi (166 paket, `pip freeze`).

## 2026-11-13 (gece, devam) — BULGU 20-22: "hazır" görünen araçları da sınadım

19. bulgudan sonra tersini kontrol ettim: `available()` **True** diyen araçlar
GERÇEKTEN çalışıyor mu? Her birine gerçek girdi verdim. Üç tanesi çalışmıyordu.
Bu en sinsi kategori, çünkü araç durumu tablosu onları "hazır" gösteriyordu.

### BULGU 20 — TMHMM hiç çalışmıyordu (shebang), buna rağmen "GERÇEK" raporlanıyordu

`tools/tmhmm-2.0d/bin/tmhmm` shebang'i `#!/usr/local/bin/perl`; bu makinede perl
`/usr/bin/perl`. Script `cannot execute: required file not found` ile hiç
başlamıyor. `available()` yalnız dosya varlığı + `which("perl")` kontrol ettiği
için True dönüyordu — ve `funnel_methods.transmembran` **"GERÇEK (TMHMM-2.0)"**
yazıyordu. Transmembran helis sayısı hiç ölçülmemişti.

Düzeltme: `tmhmm` ve `tmhmmformat.pl` shebang'leri gerçek perl yoluna çevrildi
(yedek `.orig`). Doğrulama: çözünür protein → PredHel=0, tasarlanmış 3-helisli
TM proteini → PredHel=3. Modül artık `{'cozunur': 0, 'tm': 3}` döndürüyor.

### BULGU 21 — IApred eksik bağımlılıkla çöküyordu, antijenite HER PROTEİNDE 0.0 oluyordu

`IApred.py` import zincirinde `imblearn` ve `seaborn` var; ikisi de kurulu
değildi → `ModuleNotFoundError` → CSV üretilmiyor → `iapred.predict()` boş
sözlük. Ama `available()` (script + models/ dizini) True diyordu ve huni
**"antijenite: GERÇEK (IApred)"** yazıyordu. Gerçek koşu kaydında adayların
antijenite notu `IApred (ham=0.0, )` — yani ayrıştırmanın varsayılan 0.0'ı.
Antijenite bir huni kriteri; hiç ölçülmemişti.

Düzeltme: `imbalanced-learn` + `seaborn` kuruldu. Doğrulama: iki test
proteini ayrışmış skor veriyor (-0.10 "Moderate" / 0.25 "Moderate").

### BULGU 22 — ToxinPred2 üç ayrı nedenle çalışmıyordu; düzeltmek bir KALİBRASYON sorununu açtı

Üç bağımsız kusur:
1. `RF_model` sklearn **0.22** ile kaydedilmiş → `incompatible dtype` (1.3+
   ağaç düğümlerine `missing_go_to_left` ekledi). Ek olarak 1.4+ `monotonic_cst`
   ve `n_features_in_` bekliyor; eski pickle'da yok.
2. `envfile` GELİŞTİRİCİNİN makine yollarını içeriyor
   (`/home/neelam/...`) → Model 2 (Hybrid) BLAST + MERCI adımlarında düşüyor.
3. `toxinpred2.py` pandas 1.x API'si kullanıyor: `DataFrame.append` (2.0'da
   kaldırıldı) ve `df.sum(axis=1)` (artık metin kolonunu atlamıyor).

Üçü de düzeltildi. **Dönüşüm sadakatini ölçtüm**: aracın kendi referans
setinden 150 toksik + 150 toksik-olmayan tam dizi → AUC 1.0000, duyarlılık
1.000, seçicilik 1.000. (Bu EĞİTİM verisi olduğu için mükemmel ayrım
beklenir; genellemeyi değil, modelin DOĞRU yüklendiğini kanıtlar — bozuk bir
dönüşüm bunu vermezdi.)

**Ve sonra asıl sorun çıktı.** Araç çalışınca huniye ne yaptığına baktım:

| ölçüt | tam dizi (>=35 aa) | 9-mer parçalar |
|---|---|---|
| AUC | 1.00 | 0.69 |
| eşik 0.6'da zararsız peptitin toksik denme oranı | — | 66/100 |

ToxinPred2'nin ALTI referans setinin tamamında en kısa dizi **35 aa**
(medyan 98-353). VaxForge ise huniye **8-15 aa** MHC epitopları besliyor — yani
aracın hiç doğrulanmadığı uzunluk alanı. Toksisite `hard_filter: true` ve eşik
0.6 olduğu için araç çalışır hâle gelince **güvenli epitopların üçte ikisi
sessizce elenecekti**. Araç bozuk olduğu için bu kalibrasyon sorunu bugüne
kadar görünmüyordu.

Eşik taraması da kurtarmıyor: 0.7 → seçicilik .77 / duyarlılık .49; 0.9 →
seçicilik .97 / duyarlılık .13. Araç bu uzunlukta ayırt edemiyor.

Düzeltme (eşik seçmek DEĞİL — alan sınırı):
- `toxinpred.MIN_VALIDATED_LENGTH = 35` + `in_validated_domain()`, gerekçe ve
  ölçüm kodda belgeli.
- `survival.py`: alan DIŞI peptitte skor hesaplanır ve raporlanır ama SERT
  FİLTRE uygulanmaz; yöntem etiketi açıkça uyarıyor, özete
  `toksisite_alan_disi` sayacı eklendi, not iki dilli.
- `i18n.py`: 2 yeni TR+EN anahtarı.

Eşiği değiştirmedim ve filtreyi kapatmadım — bu bilimsel bir karar ve senin.
Benim yaptığım, aracı doğrulanmadığı alanda sert filtre olarak kullanmamak.

### Kurulum betikleri (hepsi idempotent, `--check` modlu)

- `scripts/configure_dtu_wrappers.py` — 4 DTU sarmalayıcısında 6 `setenv` yolu.
- `scripts/setup_toxinpred.py` — envfile + 2 pandas yaması + model dönüşümü.
- `scripts/convert_toxinpred_model.py` — RF_model sklearn köprüsü (nodes dtype,
  values normalize, eksik öznitelik tamamlama). Yükleme BAŞARISI yetmez diye
  gerçek `predict_proba` denemesi de yapıyor.
- `scripts/convert_ifn_models.py` — zaten vardı (IFNepitope2).

Bunlar gerekli çünkü `tools/` .gitignore'da: düzeltmeler sürüm kontrolünde
saklanmıyor, paketi yeniden çıkaran aynı sessiz kusurla karşılaşır.

Testler: 4 yeni regresyon testi → **32 birim testi geçiyor**.

### BULGU 23 — ToxinPred2 TEK peptit girdisinde çöküyordu

`prediction()` içinde `np.loadtxt(...)` tek satırlık dosyada 1B dizi döndürüyor;
eski sklearn bunu tolere ediyordu, 1.x `Expected 2D array` ile düşüyor. Sonuç:
girdi 1 peptit → 0 sonuç; 2+ peptit → çalışıyor. Pipeline normalde çok peptit
beslediği için pratikte nadir, ama tek-aday bir koşuda toksisite sessizce
ölçülmemiş kalırdı. `np.atleast_2d` ile düzeltildi (davranış korunur), yama
`scripts/setup_toxinpred.py`'ye eklendi. 1/2/3 peptit girdisiyle doğrulandı.

### KÖK ÇÖZÜM — araç durumunda DERİN kontrol (`verify()`)

20-22. bulgular ortak bir dersi gösterdi: **`available()` yalan söyleyebilir.**
Üçü de yalnız dosya varlığına bakıyordu ve araç hiç çalışmazken True diyordu.
19. bulguda eklediğim araç durumu katmanı bu yalanı aynen devralıyordu.

Bu yüzden modüllere opsiyonel `verify()` ekledim: aracı GERÇEK girdiyle
çalıştırıp sonuç üretip üretmediğine bakar. `tool_status.collect()` `verify()`
varsa onu kullanır ve kontrolü `check: "verified"` olarak işaretler; yoksa
`available()` ile yetinir ve `check: "shallow"` der. Rapor bu ayrımı iki dilli
yazıyor: "'Dosya var' garanti değildir: bu projede TMHMM, IApred ve ToxinPred2
'hazır' görünüp hiç çalışmıyordu."

`verify()` eklenenler: `tmhmm_local` (0.0 s), `iapred` (1.9 s), `toxinpred`
(1.9 s, iki peptitle — tek-peptit kusuru yüzünden). Koşu başına toplam ek
maliyet ~4.4 s. Kalan 19 araç hâlâ sığ kontrolde; rapor bunu saklamıyor.

IFNepitope2 modelleri de dönüştürüldü (`scripts/convert_ifn_models.py`):
`python3 -W error::UserWarning` altında artık uyarı YOK, yani eski-pickle
sessiz-yanlış-sonuç riski kapandı. Ayrışmış çıktı veriyor (0.52 indükleyici /
0.44 değil).

Testler: 3 yeni regresyon testi → **35 birim testi geçiyor** (14.5 s; derin
kontroller gerçek araç çalıştırdığı için süre arttı).

### Koşu durumu

`rabies_2026-11-13_tam` koşusu TÜM düzeltmelerle yeniden başlatıldı. Bir önceki
denemeyi (01:10'da başlayan) 36. dakikada durdurdum: TMHMM/IApred/ToxinPred2
düzeltmelerinden önce başlamıştı, yani zaten aşılmıştı. Yarım çıktı dizini
silindi; **önceki üç tam koşu korunuyor** (temel, vekil, gerçek araç, domuz
düzeltmeli).

ÖNEMLİ: bu koşu, huninin İLK KEZ gerçek TMHMM + gerçek IApred ile çalıştığı
koşu. Transmembran ve antijenite daha önce hiç ölçülmemişti, dolayısıyla huni
çıktısının önceki koşulardan FARKLI olması beklenir — bu bir regresyon değil,
ilk kez doğru ölçüm.

### Kurulumda eksik kalan üç Python bağımlılığı (bulgu 21 ile aynı sınıf)

`requirements.txt` bunları listelemiyordu, dolayısıyla taze venv'de yoktular:

- `imbalanced-learn` + `seaborn` → IApred'in import zinciri (bulgu 21).
- `reportlab` → `report_pdf`. Eksikken pipeline dürüstçe `pdf_error.txt`
  yazıyordu ("PDF üretilemedi: No module named 'reportlab'"), yani bu kusur
  sessiz DEĞİLDİ — iyi tasarım. Ama `benchmark/run_*.py` çıktı kopyalama adımı
  `pdf_error.txt`'yi kopyalamıyor, bu yüzden korunan benchmark dizinlerinde PDF
  yokluğunun nedeni görünmüyor. Kuruldu; artık PDF üretilebiliyor.

Hepsi `requirements.lock.txt`'e girdi (169 paket). Lock dosyası `pip freeze`
çıktısı olduğu için `ruff` gibi geliştirme araçlarını da içeriyor.

### Ölü i18n anahtarı (ruff F601 ile bulundu)

`i18n.py` içinde `"rep_methods"` İKİ KEZ tanımlıyordu (satır 84 ve 264). Python
sözlük literalinde SONRAKİ kazanır, yani ilk tanım ("Kullanılan yöntemler ve
araçlar" / "Methods and tools used") hiçbir zaman kullanılmıyordu — sessizce
gölgelenmiş bir metin. Ölü olanı kaldırdım; görünen başlık DEĞİŞMEDİ
("Kullanılan yöntemler" / "Methods used"). Anahtar sayısı 163.

Ayrıca ruff'ın işaretlediği 4 kullanılmayan import temizlendi. 55 vaxforge
modülünün tamamı import ediliyor. Yeni yazdığım dosyalar (tool_status.py, üç
kurulum betiği, diagnostics.py) ruff'tan temiz geçiyor. Kalan 53 uyarı
önceden var olan stil meseleleri (noktalı virgül, import sırası, zip strict).

### Benchmark kayıtları artık araç durumunu saklıyor

`benchmark/run_rabies.py` ve `run_bacterium.py` özetlerine `funnel_methods`,
`tool_status`, `tool_status_summary` ve `diagnostics_summary` eklendi. Gerekçe:
10-bakteri benchmark'ının (2026-10) özetleri bu alanları İÇERMİYOR, bu yüzden
TMHMM/IApred/ToxinPred2'nin o koşularda çalışıp çalışmadığı dosyalardan
ANLAŞILAMIYOR. Her iki README'ye bu ikinci uyarı eklendi (TR + EN): VaxForge
rakamı benchmark yeniden koşulana kadar GEÇİCİ sayılmalı, rakipler kendi
konteynerlerinde koştuğu için etkilenmez.

## 2026-11-13 (gece, devam 2) — PSORTb AÇILDI, Gram düzeltmesi GERÇEK araçla doğrulandı

### PSORTb engeli sudo olmadan aşıldı

Engel şuydu: `podman pull` çıkarma sırasında `potentially insufficient UIDs or
GIDs available in user namespace (requested 508:1000025)` veriyordu; `/etc/subuid`
aralığı `ali:100000:65536`, yani istenen 1000025 aralık DIŞINDA. Çözüm olarak
`sudo usermod --add-subuids` bekliyordum.

Gerek kalmadı: `podman --storage-opt ignore_chown_errors=true pull` ile image
çekildi. Bu, overlay sürücüsünün belgelenmiş bir seçeneği — aralık dışı UID'ler
için chown hatalarını yok sayar. Çekildikten sonra **sade `podman run` çalışıyor**,
yani projenin `psortb.py` kodu değiştirilmeden işliyor.

Doğrulama: `PSORTb version 3.0` · modül `predict()` per-dizi lokalizasyon +
skor döndürüyor · tek dizi ~2.8 s.

**Sudo gerektiren açık madde KAPANDI.** `usermod` artık gerekmiyor.

### Gram-pozitif düzeltmesi artık KAYITLA DEĞİL, GERÇEK ARAÇLA doğrulandı

Oturumun başında yaptığım Gram düzeltmesini yalnızca kayıtlı çıktı üzerinden
test edebilmiştim. Şimdi iki kaçan antijeni UniProt'tan çekip gerçek PSORTb'den
geçirdim:

| antijen | aksesyon | PSORTb (gerçek) | skor |
|---|---|---|---|
| ActA (*L. monocytogenes*) | P33379 | CytoplasmicMembrane | 10.0 |
| CapB (*B. anthracis*) | P19580 | CytoplasmicMembrane | 10.0 |

`MISS_ATTRIBUTION.md`'nin kaydıyla BİREBİR aynı — bağımsız doğrulama. İkisi de
`membrane` etiketine eşleniyor ve yeni Gram-pozitif izin listesinde bu etiket
var, yani artık geçiyorlar. Gram-negatifte geçmiyorlar (doğru: orada
sitoplazmik membran dış membranla perdeli). Gram bilinmiyorsa temkinli
Gram-negatif listesi uygulanıyor.

Yeni test gerçek araç çıktısını sabit olarak gömüyor, dolayısıyla podman
olmayan ortamda da koşuyor.

### PSORTb'ye de derin kontrol eklendi

`psortb.available()` yalnız `podman image exists` kontrol ediyordu — aynı sığ
sınıf. `verify()` eklendi: gerçek konteyner çağrısı yapar (~3 s). Araç durumu
artık **11/22 hazır, 4'ü derin doğrulanmış** (tmhmm, iapred, toxinpred, psortb),
toplama süresi 8.8 s.

### DeepLoc bu makinede İMKANSIZ (disk değil, RAM)

ESM1b ağırlığı 7.29 GB (HEAD ile ölçtüm), sistemde toplam 7.5 GB RAM var.
`torch.load` arşivin tamamını çözdüğü için tepe bellek dosya boyutu
civarındadır → pipeline'ın geri kalanına yer kalmaz, büyük olasılıkla OOM.
Disk sorun değil (940 GB boş), ağ da değil. Yani DeepLoc'u indirmek çözüm
olmaz; GPU/daha fazla bellekli makine gerekir. Bakteri profilinde zaten birincil
araç PSORTb (artık çalışıyor); DeepLoc virüs/parazit kolunda yedek.

### VFDB — somut alternatif var ama YÖNTEM değişikliği (senin kararın)

`www.mgc.ac.cn` bu ortamda hiç çözülmüyor (github.com ve rest.uniprot.org
çözülüyor). GitHub'da VFDB aynaları var ama hepsi doğrulanmamış üçüncü-taraf
kopya — makale için köken sorunu olur, bu yüzden kullanmadım.

Erişilebilir ve atıf verilebilir alternatif: UniProt `KW-0843` (Virulence)
anahtar kelimesi — **2.992 doğrulanmış bakteriyel protein** (tümü 5.075).
VFDB setB'den küçük ama Swiss-Prot küratörlü. Bunu kurmak Yöntemler'de
"VFDB" → "UniProt virülans anahtar kelimesi" demek olur, yani manuscript'i
etkiler. Karar senin; şu an pipeline paket-içi mini referansa düşüyor ve bunu
dürüstçe "0 VFDB kanıtlı" diye raporluyor.

Testler: 2 yeni regresyon testi → **37 birim testi geçiyor**.

## 2026-11-13 (gece, devam 3) — TAM KOŞU BİTTİ + iki DÜZELTME kendi iddiamda

### Önceki raporumda İKİ fazla-okuma vardı; düzeltiyorum

**(a) "60/85 adayda bcell_score NaN" BİR KUSUR DEĞİLDİ.** Tip dağılımını
kontrol ettim: 25 B + 30 MHC-I + 30 MHC-II = 85. B-hücre skoru yalnız B
epitoplarına verilir; 60 NaN tam olarak MHC epitoplarıdır. 25 B epitobunun
TAMAMINDA gerçek skor var. Tasarım gereği, benim yanlış okumam.

**(b) "Bozuk BepiPred çalışan yedeği DEVRE DIŞI bıraktı" DOĞRU DEĞİL.**
`epitope.py` `use_bp = len(bp) == len(s)` ile karar veriyor; BepiPred boş liste
dönünce `use_bp` False olur ve klasik Kolaskar-Tongaonkar+Parker yedeğine
DÜŞER, üstelik yöntemi dürüstçe etiketler. Yani yedek doğru çalıştı.

Gerçek etki, düzeltmenin DAHA İYİ yöntemi geri getirmesiydi (yalan düzeltmek
değil). Ölçülebilir:

| koşu | B-hücre skoru (en yüksek 4) | kullanılan yöntem |
|---|---|---|
| gercek (BepiPred bozuk) | 0.589 / 0.589 / 0.576 / 0.571 | klasik Kolaskar+Parker |
| tam (BepiPred çalışır)  | 0.916 / 0.914 / 0.913 / 0.909 | GERÇEK BepiPred-1.0 |

16. bulgunun geçerli kısmı şudur: `bepipred.available()` sarmalayıcı yapılandırılmamışken
True dönüyordu, dolayısıyla ARAÇ DURUMU tablosu yalan söylüyordu. `_configured()`
düzeltmesi bu yüzden doğru. Ama aday-düzeyi yöntem etiketi dürüsttü.

### BULGU 24 — B epitobunun yöntem etiketi ÇIKTIYA HİÇ GİRMİYORDU

Yukarıdaki (b)'yi doğrularken buldum: `report.py` aday satırını kurarken
`"method": method_label(p.methods.get("mhc_score"), lang)` yazıyordu — YALNIZ
MHC yöntemi. B epitoplarında `mhc_score` yöntemi olmadığı için `method` alanı
**NaN** kalıyordu.

Yani `run.json` ve `candidates.csv` okuyan biri, B-hücre skorunun gerçek
BepiPred'den mi klasik yedekten mi geldiğini GÖREMİYORDU. Dürüst etiket bellekte
üretiliyor ama dışa aktarılmıyordu. BepiPred'in sessizce yedeğe düşmesini tam
olarak bu boşluk gizledi. Doğrulandı: üç koşunun 25 B adayının hepsinde
`method: nan`; MHC adaylarında etiket dolu.

Düzeltme: `p.methods.get("mhc_score") or p.methods.get("bcell_score")`. Test
eklendi.

### TAM KOŞU — beş sütunlu karşılaştırma

| ölçüt | temel | vekil | gerçek araç | domuz fix | TAM |
|---|---|---|---|---|---|
| süre (s) | 4472 | 262 | 2439 | 2386 | 2493 |
| G rank | #2/5 | #1/5 | #1/5 | #1/5 | #1/5 |
| candidacy | 0.7102 | 0.8244 | 0.7787 | 0.8246 | 0.8246 |
| insan/sığır/fare/domuz | 31/33/26/19 | 56/56/56/30 | 54/53/47/0 | 55/53/46/28 | 55/53/46/28 |
| J(insan~sığır) | 0.75 | 1.00 | 0.75 | 0.833 | 0.833 |
| J(insan~domuz) | 0.50 | 0.667 | 0.00 | 0.545 | 0.545 |
| **MEV antijenitesi** | — | — | **0.0** | **0.0** | **0.664** |
| araç durumu kaydı | yok | yok | yok | yok | **10/22 (derin 3)** |

**En net kanıt MEV antijenitesi satırı.** Önceki koşularda 0.0 çıkıyordu çünkü
IApred eksik bağımlılıkla çöküyordu (bulgu 21). Düzeltmeden sonra 0.664. Yani
MEV konstruktunun antijenitesi bugüne kadar hiç ölçülmemişti ve sıfır olarak
raporlanıyordu.

Toksisite 85 → 85 (0 elendi), beklendiği gibi: peptitler 8-15 aa, yani
ToxinPred2'nin doğrulanmış alanının altında → skor hesaplanır, SERT FİLTRE
uygulanmaz (bulgu 22). Araç artık çalıştığı hâlde güvenli epitopları elemiyor.

Yapısal doğrulama bu koşuda atlandı: ESMAtlas API hatası. Dürüstçe raporlandı,
yerel bir kusur değil.

PSORTb bu koşudan SONRA açıldı, bu yüzden araç durumu 10/22 (şimdi 11/22).
Virüs profili PSORTb kullanmadığı için koşuyu etkilemez.

Testler: 1 yeni regresyon testi → **38 birim testi geçiyor**.
