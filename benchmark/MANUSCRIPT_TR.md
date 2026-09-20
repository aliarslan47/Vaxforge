**VaxForge: konak-özgü epitop çözünürlüğü ile uçtan uca ters aşı bilişimi**

Ali Arslan

Department of Molecular Biology and Genetics, Institute of Natural Sciences, Yıldız Technical University, Istanbul 34220, Turkey

Yazışma: aliarslan47@gmail.com

Taslak: 2026-09-20. Yöntem anlatısı `MANUSCRIPT_METHODS.md`'de; ground-truth antijen kanıtı ve DOI'ler `PROTECTIVE_ANTIGENS.md`'de; aksesyonlar `PATHOGEN_ACCESSIONS.md`'de. Bu belge, İngilizce sürümün (`MANUSCRIPT.md`) Türkçe karşılığıdır.

---

**Özet**

Ters aşı bilişimi araçları bir genomdan aday antijen sıralar; ancak çoğu bu sıralamada durur ve konağı tek bir referans sayar. Bu iki eksik, hem bilinen antijenleri geri getirme başarımını hem de veteriner ve türler-arası hedeflerdeki kullanışlılığı düşürür. Burada VaxForge'u tanıtıyoruz: bir patojeni aday keşfinden yumuşak-filtreli antijen hunisine, ardından beş konak türü için ayrı ayrı sınıf I ve sınıf II epitop tahminine ve karakterize edilmiş çok-epitoplu bir konstrukta taşıyan, in silico immünizasyonla biten bir hat. VaxForge'u iki yerleşik araçla, NERVE 2.0 ve Vaxign-ML ile, varsayılan parametreler ve özdeş girdi altında karşılaştırdık. On beş ruhsatlı ya da preklinik koruyucu antijenin ground-truth alındığı dört bakteriyel proteomda VaxForge 15 antijenin 14'ünü geri getirdi. Bu değer Vaxign-ML'inkiyle (15/15) yarışır ve NERVE 2.0'ın (9/15) belirgin önündedir; üstelik VaxForge koruyucu antijenleri daha üst sıralara koydu (örneğin demir-düzenlenen yüzey belirleyicisi B, 2.889 protein içinde 4. sırada). VaxForge'un kaçırdığı tek antijen olan ActA, bağımsız araçta da en dipte çıktı; yani bir yapaylık değil, gerçekten zor bir hedef. NERVE 2.0 daha yüksek fold-zenginleşmeye ulaştı, ama bunu dört ruhsatlı meningokok antijeninden üçünü her ölçekte eleyerek yaptı; dolayısıyla fold-zenginleşme tek başına yanıltıcıdır. Bir kuduz modeli ve iki virüste VaxForge, türden türe belirgin biçimde değişen konak-özgü glikoprotein epitop repertuarları verdi (ikili Jaccard örtüşmesi birin altında, bir durumda sıfır); bu, karşılaştırılan araçların ulaşamadığı bir boyuttur. Kısacası VaxForge, güçlü bir bağımsız tahmincinin recall'ını yakalarken alanın bugün sunmadığı konak çözünürlüğünü ve uçtan uca bir konstruktu ekler.

---

**Giriş**

Ters aşı bilişimi serogrup B *Neisseria meningitidis* ile doğdu. Pizza ve arkadaşları (2000), aday antijenleri kültür bileşenlerinden değil doğrudan tüm-genom dizisinden çıkararak bu yöntemi kurdu ve ruhsatlı dört-bileşenli meningokok aşısının antijenlerini, yani rekombinant faktör H bağlayıcı protein, Neisseria adhezini A ve Neisseria Heparin Bağlayıcı Antijen ile PorA taşıyan dış-membran vezikülünü ortaya çıkardı (Giuliani vd., 2006). Aradan geçen yirmi yılda yöntem olgun yazılımlara kavuştu. NERVE 2.0 güncel ve yerel kurulabilen bir ters aşı bilişimi ortamı sunar (Conte vd., 2024); Vaxign-ML ise bakteriyel koruyucu antijenleri denetimli makine öğrenmesiyle tahmin eder (Ong vd., 2020).

Bu araçların pratikte önemli üç ortak eksiği var. Birincisi, sert bir virülans ya da lokalizasyon filtresi aday listesinin görünürdeki zenginliğini artırırken gerçek antijenleri sessizce eleyebilir. Örneğin Neisseria Heparin Bağlayıcı Antijen küratörlü virülans veritabanlarında yoktur, ama ruhsatlı bir aşı bileşenidir; virülans-veritabanı isabeti şart koşan her hat onu kaybeder. İkincisi, epitop çıktısı çoğunlukla tek bir insan referansına göre hesaplanır. Oysa sığır, domuz ya da farenin asıl konak olduğu pek çok veteriner ve zoonotik hedef, dürüst bir sunum modeli olmadan kalır. Üçüncüsü, araçların çoğu sıralı bir protein listesi verip durur; bir antijeni immünojene dönüştüren adımlar, yani epitop seçimi, konstrukt kurulumu ve ilk immünolojik kontrol, kullanıcıya bırakılır.

VaxForge bu üç eksiği tek bir deterministik hatta kapatır. Sert virülans filtresini yumuşak ve ağırlıklı bir filtreyle değiştirir; böylece virülans-veritabanı kaydı olmayan gerçek bir yüzey antijeni de geçebilir. Kullanıcının seçtiği konak için sınıf I ve sınıf II bağlanmayı, insan, fare, sığır, domuz ve tavuk allel panellerinden tahmin eder; böylece epitop repertuarı hedeflenen türü yansıtır. Ve sıralamayla yetinmeyip çok-epitoplu bir konstrukta ve in silico immünizasyon çıktısına uzanır. Aşağıda ilk iki katkıyı dört bakteri ve üç virüste NERVE 2.0 ve Vaxign-ML'e karşı ölçüyor, ayrıca çalışma sırasında yakalayıp giderdiğimiz iki ölçek kaynaklı kusurun oluşturduğu hat doğrulamasını sunuyoruz.

---

**Bulgular**

**VaxForge hattı.** VaxForge bir proteomu, kodlayan dizileri ya da ham bir genomu otomatik tanımayla alır; ardından keşif, antijen hunisi, konak-duyarlı epitop tahmini, konstrukt kurulumu ve immünizasyon kontrolü gelir. Keşif adımı proteinleri DIAMOND (Buchfink vd., 2021) ile VFDB'ye (Liu vd., 2022) hizalar. Huni; hücre-altı yerleşimi PSORTb (Yu vd., 2010) ve DeepLoc (Thumuluri vd., 2022), transmembran topolojiyi TMHMM (Krogh vd., 2001), sinyal peptidini SignalP (Almagro Armenteros vd., 2019), antijeniteyi IApred (Miles vd., 2025) ile değerlendirir ve otoimmün güvenlik için zorunlu bir insan-homoloji elemesi uygular; birkaç sert kesim dışında kalanları ağırlıklı bir adaylık puanıyla sıralar.

Huniyi geçen proteinler hattın çekirdeğine, immünolojik katmana girer. NetMHCpan ve NetMHCIIpan (Reynisson vd., 2020) sınıf I ve sınıf II bağlanmayı seçilen konağın allelleri için skorlar; böylece bir T-hücre epitobu ancak hedeflenen türde sunuluyorsa tutulur. T-hücre okumasını üç ek yordayıcı keskinleştirir: sınıf I immünojenite (Calis vd., 2013), proteozomal kesim ile TAP taşınımını birleştiren antijen işleme (Larsen vd., 2005) ve antikor hedefleri için BepiPred-3.0 ile lineer B-hücre epitopları (Clifford vd., 2022). Seçilen epitoplar için VaxForge insan HLA allellerinde popülasyon kapsamını hesaplar; allel frekans verisi bulunmayan konaklarda bu değeri uydurmak yerine "veri yok" olarak bildirir. Tutulan epitoplar, tanımlı linkerlerle (sitotoksik epitoplar arası AAY, yardımcı epitoplar arası GPGPG, B-hücre epitopları için KK) ve bir adjuvanla çok-epitoplu bir konstrukta birleştirilir; konstrukt, antikor ve izotip yanıtını zaman içinde benzeten in silico immünizasyon adımına geçer. Her skorun kaynağı belli olsun diye her yöntem atıfıyla raporlanır.

**VaxForge dört ruhsatlı meningokok antijeninin dördünü de bulur, NERVE 2.0 yalnız birini.** İlk karşılaştırmayı serogrup B *N. meningitidis* (MC58 suşu, proteom UP000000425) üzerinde, dört 4CMenB antijenini ground-truth alarak yaptık. Genişletilmiş 204-protein setinde (dört koruyucu, 200 arka plan, sabit seed) VaxForge dört antijenin dördünü de buldu (recall 4/4; 68 aday; fold 3,12; Fisher p = 0,012); sıralama PorA birinci, faktör H bağlayıcı protein üçüncü, Heparin Bağlayıcı Antijen altıncı, adhezin A sekizinci. Vaxign-ML de dördünü buldu (59 aday); NERVE 2.0 yalnız PorA'yı buldu (recall 1/4). Aynı tablo, kompakt 64-protein setinde ve iki hızlı araç için tam 2.003-protein proteomunda da tekrar etti (Tablo 1; Şekil 1 ve 2). En kritik örnek, virülans-veritabanı negatif olan Heparin Bağlayıcı Antijendir: VaxForge onu yumuşak huniden geçirdi, NERVE 2.0 ise sert seçim adımında eledi. Bağımsız bir makine-öğrenmesi aracının, Vaxign-ML'in, aynı 4/4'e ulaşması, bu davranışın tek araca özgü bir tesadüf değil, yumuşak filtrelemenin bir sonucu olduğunu gösterir.

*Tablo 1. Üç girdi ölçeğinde meningokok karşılaştırması (varsayılan parametreler, özdeş girdi).*

| Set | Araç | Recall | Aday | Fold | Fisher p | AUC |
|---|---|---|---|---|---|---|
| 204 | VaxForge | 4/4 | 68 | 3,12 | 0,012 | 0,99 |
| 204 | Vaxign-ML | 4/4 | 59 | 3,64 | 0,007 | 0,97 |
| 204 | NERVE 2.0 | 1/4 | 3 | 25,0 | 0,058 (ad) | yok |
| Tam 2.003 | Vaxign-ML | 4/4 | 613 | 3,28 | yok | yok |
| Tam 2.003 | NERVE 2.0 | 1/4 | 9 | 62,5 | yok | yok |

**Yüksek fold-zenginleşme yüksek recall demek değildir.** NERVE 2.0 kısa aday listeleriyle yüksek fold verdi (204 setinde 25, tam proteomda 62,5); bu, yaklaşık 16,6'lık yayınlanmış değeriyle uyumlu. Ne var ki aday sayısı çok küçük olduğundan Fisher p'si anlamlı değildi ve araç her ölçekte dört ruhsatlı antijenden üçünü kaçırdı (Şekil 3). Bu üç kayıp, faktör H bağlayıcı protein, adhezin A ve Heparin Bağlayıcı Antijen, NERVE 2.0'ın kendi sürekli skorunun düşüklüğünden gelmedi; tersine bu skor onları yüksek verdi. Kayba yol açan, skordan sonra gelen sert seçim adımıydı (periplazmik yerleşim, 0,5 eşiğinin altında adhezin olasılığı ve tek transmembran heliksli bilinmeyen yerleşim). Demek ki fold-zenginleşme seçiciliği ödüllendirir ve recall düşerken bile yükselebilir; bu yüzden yalnız başına okununca aldatıcıdır. Bir aracın gerçekte ne verdiğini gösteren ölçüt, bilinen antijenlerdeki recall ve bu antijenlerin sıralarıdır.

**Kalıp genom ölçeğinden bağımsız ve bir doğrulama adımı sonuçları güvenilir kıldı.** Tek örnekle sınırlı kalmamak için, artan genom boyutu ve Gram dengesine göre seçtiğimiz üç zoonotik bakteriyi tam proteom ölçeğinde koştuk (Tablo 2). VaxForge *Staphylococcus aureus* için 5/5 (demir-düzenlenen yüzey belirleyicisi B, 2.889 içinde 4. sırada), *Listeria monocytogenes* için 4/5 ve *Salmonella enterica* serovar Typhimurium için 5/5 buldu; toplam 14/15 (Şekil 5). Vaxign-ML 15/15, NERVE 2.0 ise 9/15 bularak meningokoktakiyle aynı türden yüzey antijenlerini yine kaçırdı. Her iki aracın da tam recall'a ulaştığı iki türde VaxForge koruyucu antijenleri Vaxign-ML'den daha üste yerleştirdi. VaxForge'un tek kaybı olan *L. monocytogenes* ActA'sı gerçekten zor bir antijendir; Vaxign-ML listesinde de en dipte çıktı, yani hedefin zorluğunu yansıtır, bir hat kusurunu değil. Bu koşular sırasında, ilk denemelerde VaxForge recall'ını yapay olarak düşüren iki ölçek kaynaklı kusuru saptayıp giderdik (ilk değerler 3/5, 3/5 ve 0/5); lokalizasyon ve MHC adımları tam-proteom ölçeğinde tek ve aşırı büyük bir yığın altında çöküyordu, ikisini de sınırlı ve paralel alt-yığınlara böldük (`MANUSCRIPT_METHODS.md`, bölüm B). Kontrol deneyleri toksisite ve alerjenite filtrelerinin suçsuz olduğunu doğruladı.

*Tablo 2. Tam proteom ölçeğinde üç-bakteri karşılaştırması; her bakteride beş koruyucu antijen.*

| Bakteri (genom) | Araç | Recall | Aday | Fold | AUC |
|---|---|---|---|---|---|
| *S. aureus* (2,8 Mb, 2.889) | VaxForge | 5/5 | 729 | 3,98 | 1,00 |
|  | Vaxign-ML | 5/5 | 838 | 3,46 | 0,94 |
|  | NERVE 2.0 | 3/5 | 34 | 55,8 | yok |
| *L. monocytogenes* (2,9 Mb, 2.844) | VaxForge | 4/5 | 478 | 4,79 | yok |
|  | Vaxign-ML | 5/5 | 917 | 3,11 | 0,95 |
|  | NERVE 2.0 | 3/5 | 73 | 24,3 | yok |
| *S.* Typhimurium (4,9 Mb, 4.533) | VaxForge | 5/5 | 1.368 | 3,32 | 0,99 |
|  | Vaxign-ML | 5/5 | 1.276 | 3,56 | 0,92 |
|  | NERVE 2.0 | 3/5 | 118 | 23,6 | yok |
| **Toplam** | **VaxForge 14/15 · Vaxign-ML 15/15 · NERVE 2.0 9/15** |  |  |  |  |

**Konak seçimi öngörülen epitop repertuarını değiştirir; karşılaştırılan araçların ulaşamadığı boyut.** VaxForge bağlanma tahminlerini beş konak için küratörlü allel panellerinden alır: insan, fare, sığır, domuz ve tavuk (Tablo 3). İnsan panelleri IMGT/HLA'dan, diğerleri IPD-MHC'den derlenmiştir. İnsan, fare ve sığırda hem sınıf I hem sınıf II gerçek NetMHCpan ve NetMHCIIpan modelleriyle koşar; sığırın sınıf II BoLA-DRB3 paneli de gerçektir. Domuzda sınıf I gerçek, sınıf II ise tek doğrulanmış allele dayanır; tavuk ise dizi-tabanlı vekile düşer, çünkü bağlanma yordayıcıları kanatlı allellerini kabul etmez. Paneller tam veritabanları değil temsili alt kümelerdir ve geçersiz allel adları çalışma-zamanında elenir; dolayısıyla bildirilen epitoplar temkinlidir.

Virüs kolunu bu konak-özgü sunumu sınamak için kurduk. Tek yapısal proteini olan glikoprotein G'nin nötralize edici ve koruyucu antikorları uyaran tek protein olduğu (Cox vd., 1977) kuduz virüsünde (RefSeq NC_001542), VaxForge G'yi beş protein içinde ikinci sıraya koydu; daha da önemlisi dört doğal memeli konağın her biri için farklı bağlı-epitop setleri üretti (glikoprotein G: farede 8, insanda 7, sığırda 7, domuzda 5). Glikoprotein G epitoplarının ikili konak örtüşmesi birin belirgin altındaydı (Tablo 4); yani konak seçimi öngörülen repertuarı gözle görülür biçimde değiştirdi. İki ek dört-konaklı virüs bu etkiyi yineledi: kuduzun rabdovirüs akrabası veziküler stomatit virüsü glikoprotein G'yi yine 5 içinde 2. sıraya koydu ve insan ile domuz örtüşmesi 0,33'e kadar indi; ensefalomiyokardit virüsü kapsid VP1'i 11 içinde 3. sıraya koydu ve sığır ile domuz arasında hiç ortak epitop çıkmadı. NERVE 2.0 yalnız bakteriler için çalışır, Vaxign-ML'de konak boyutu yoktur (kuduz glikoprotein G'yi konak ayrımı olmadan 5 içinde 4. sıraya koydu); dolayısıyla bu yetenek üç araç arasında yalnız VaxForge'da vardır. Kuduz tek koruyucu antijen sunduğundan buradaki sıra istatistiği çıkarımsal değil betimseldir (permütasyon p = 0,40, n = 1); popülasyon kapsamı da insan-dışı konaklarda uydurulmadan "veri yok" olarak bildirildi.

*Tablo 3. Epitop tahmininde kullanılan konak allel panelleri (temsili alt kümeler; tam paneller `config/hosts.yaml`).*

| Konak | Sınıf I paneli (n) | Sınıf II paneli (n) | Yordayıcı durumu |
|---|---|---|---|
| İnsan (HLA) | HLA-A ve HLA-B (9) | HLA-DRB1 (6) | NetMHCpan ve NetMHCIIpan, gerçek |
| Fare (H-2) | H-2 K/D (5) | H-2 IA/IE (5) | gerçek |
| Sığır (BoLA) | BoLA sınıf I (8) | BoLA-DRB3 (8) | gerçek |
| Domuz (SLA) | SLA sınıf I (6) | SLA-DRB1 (1) | sınıf I gerçek; sınıf II tek allel |
| Tavuk (BF) | BF2 (1) | BLB2 (1) | yordayıcı kabul etmez, vekil |

*Tablo 4. Üç dört-konaklı virüste konak-özgü epitop çözünürlüğü (VaxForge).*

| Virüs | Ground-truth antijen | Sıra | Konak-çifti epitop örtüşmesi (Jaccard) |
|---|---|---|---|
| Kuduz virüsü | glikoprotein G | 2/5 | insan~sığır 0,75; insan~fare 0,67; insan~domuz 0,50 |
| Veziküler stomatit virüsü | glikoprotein G | 2/5 | insan~sığır 0,78; insan~domuz 0,33 |
| Ensefalomiyokardit virüsü | kapsid VP1 | 3/11 | insan~fare 0,78; sığır~domuz 0,00 |

---

**Tartışma**

VaxForge, bilinen koruyucu antijenleri geri getirmede güçlü bir bağımsız tahminci kadar başarılı oldu ve o tahmincinin sunmadığı iki yeteneği ekledi. Dört bakteride 15 doğrulanmış antijenin 14'ünü buldu; bu, Vaxign-ML için 15/15, NERVE 2.0 için 9/15'e karşılık gelir ve VaxForge antijenleri en az onlar kadar, çoğu yerde daha üst sıraya koydu. NERVE 2.0'ın eksiği belirli bir nedene dayanıyor: sert seçim adımı, kendi sürekli skorunun yüksek verdiği gerçek yüzey antijenlerini eledi; bu da fold-zenginleşmeyi şişirirken recall'ı kesti. Fold-zenginleşme bir hattı olduğundan keskin gösterebilir; bu yüzden bu tür karşılaştırmalarda güvenilecek sayı, bilinen antijenlerdeki recall'dır.

İki özellik VaxForge'u recall'ın ötesinde ayırır. Yumuşak virülans filtresi, virülans-veritabanı kaydı olmayan ruhsatlı bir bileşeni, Neisseria Heparin Bağlayıcı Antijeni, geri kazandı; konak-duyarlı epitop katmanı ise kuduz, veziküler stomatit ve ensefalomiyokardit virüsleri için türe özgü repertuarlar verdi. İkincisi, hedef konağın insan olmadığı geniş veteriner ve zoonotik alan için önemlidir; orada tek bir insan referansı hangi peptitlerin sunulduğunu yanlış gösterir. VaxForge ayrıca sıralamadan sonra konstrukt kurulumuna ve in silico immünizasyon çıktısına geçer; böylece elde bir liste değil, aday bir immünojen kalır.

Birkaç sınırlama bu iddiaları çerçeveler. ActA'yı VaxForge kaçırdı; bağımsız aracın da en dibe koyduğu gerçekten zor bir antijen olsa da, geri kazanımı lineer-epitop bir hattın modellemediği konformasyonel ya da bağlama-özgü sinyalleri gerektirir. Varsayılan toksisite filtresi, toksin antijenlerinden gelen peptitleri eleyebilir; bu, ölçek düzeltmelerinden önce por-açan listeriolizin O'da peptit-düzeyi recall'ını düşürmüştü ve karşılaştırma için ayarlamadığımız temkinli bir seçim olarak kalıyor. AUC değerlerini başlık değil, açık bir küçük-örneklem uyarısıyla ikincil olarak veriyoruz; çünkü patojen başına yalnız dört ila beş pozitifle farklar istatistiksel olarak ayrışmıyor, bu temelde üstün ayrım iddia etmiyoruz. Bakteriyel karşılaştırma bir ölçek asimetrisi de taşır: VaxForge'un tam-proteom meningokok koşuları bugünkü CPU-yalnız donanımda pratik değil, bu yüzden amiral meningokok karşılaştırması araçlar arası adil olan eşlenmiş alt-örneklem setlerini kullandı; tam proteom yalnız hızlı araçlarda koşulabildi. Son olarak, B-hücre epitoplarının kabaca onda dokuzu süreksizdir ve buradaki lineer konstrukt konformasyonel epitopları temsil etmez; bunu çözmek yerine bir kapsam sınırı olarak belirtiyoruz.

İki yön kendiliğinden çıkıyor. Bir GPU arka ucu, tam-proteom VaxForge koşularına ve her patojende tam ölçekli karşılaştırmaya olanak verir. Konak-özgü epitop tahminleri, özellikle viral glikoproteinlerde görülen türler-arası ayrışma, bağlanma ve immünizasyon deneyleriyle sınanabilecek somut hipotezlerdir. Hat ve karşılaştırması, bu iki genişlemenin de yeniden tasarım gerektirmeden mevcut çerçeveye oturacağı biçimde kurulmuştur.

---

**Şekiller**

![Şekil 1](figures/fig1_antigen_recovery.png)

*Şekil 1. VaxForge ve Vaxign-ML dört ruhsatlı meningokok antijeninin dördünü de bulur, NERVE 2.0 yalnız PorA'yı bulur.* Dört 4CMenB antijeni için antijen-araç geri getirme matrisi; dolu hücreler bulunan antijenleri, kırmızı hücreler NERVE 2.0 seçim adımınca elenenleri gösterir. Recall ve aday sayıları varsayılan parametrelerle özdeş girdi üzerinde hesaplandı.

![Şekil 2](figures/fig2_recall_bar.png)

*Şekil 2. Bilinen antijenlerdeki recall, araçları AUC'den daha net ayırır.* Dört meningokok antijeninin araç-başına recall'ı, ikincil AUC ek notuyla; AUC küçük-örneklem uyarısı taşır ve başlık ölçüt olarak kullanılmaz (n = 4 pozitif).

![Şekil 3](figures/fig3_fold_vs_recall.png)

*Şekil 3. Yüksek fold-zenginleşme recall'ı değil seçiciliği izler.* Üç araç için 204-protein setinde (dolu noktalar) ve tam proteomda (içi boş noktalar) fold-zenginleşmeye karşı recall; kesikli çizgi yayınlanmış NERVE 2.0 fold değeri 16,6'yı gösterir. NERVE 2.0 yüksek fold ve düşük recall'da, VaxForge ile Vaxign-ML ise tam recall ve orta fold'da yer alır.

![Şekil 4](figures/fig4_multihost_jaccard.png)

*Şekil 4. Konak seçimi kuduz glikoprotein G epitop repertuarını yeniden biçimlendirir.* İnsan, sığır, fare ve domuz arasında glikoprotein G sınıf I epitoplarının ikili Jaccard örtüşmesi; birin altındaki değerler öngörülen repertuarın seçilen konağa bağlı olduğunu gösterir.

![Şekil 5](figures/fig5_bacteria_panel.png)

*Şekil 5. Recall kalıbı üç genom ölçeğinde de sürer.* Tam proteom ölçeğinde *S. aureus*, *L. monocytogenes* ve *S.* Typhimurium için her aracın beş koruyucu antijendeki recall'ı, AUC ek notuyla; VaxForge toplam 14/15, Vaxign-ML 15/15 ve NERVE 2.0 9/15.

Not: şekil panelleri şu an Türkçe eksen ve başlık metni taşıyor; hedef dergi isterse `make_figures.py`'den İngilizce etiketli sürümler üretilebilir.

---

**Kaynaklar**

Kaynaklar APA 7 biçiminde ve DOI bağlantılıdır; İngilizce sürümle (`MANUSCRIPT.md`) aynıdır.

Almagro Armenteros, J. J., Tsirigos, K. D., Sønderby, C. K., Petersen, T. N., Winther, O., Brunak, S., von Heijne, G., & Nielsen, H. (2019). SignalP 5.0 improves signal peptide predictions using deep neural networks. *Nature Biotechnology, 37*, 420–423. https://doi.org/10.1038/s41587-019-0036-z

Buchfink, B., Reuter, K., & Drost, H. G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods, 18*, 366–368. https://doi.org/10.1038/s41592-021-01101-x

Calis, J. J. A., Maybeno, M., Greenbaum, J. A., Weiskopf, D., De Silva, A. D., Sette, A., Keşmir, C., & Peters, B. (2013). Properties of MHC class I presented peptides that enhance immunogenicity. *PLoS Computational Biology, 9*(10), e1003266. https://doi.org/10.1371/journal.pcbi.1003266

Clifford, J. N., Høie, M. H., Deleuran, S., Peters, B., Nielsen, M., & Marcatili, P. (2022). BepiPred-3.0: Improved B-cell epitope prediction using protein language models. *Protein Science, 31*(12), e4497. https://doi.org/10.1002/pro.4497

Conte, A., Gulmini, N., Costa, F., Cartura, M., Bröhl, F., Patanè, F., & Filippini, F. (2024). NERVE 2.0: boosting the new enhanced reverse vaccinology environment via artificial intelligence and a user-friendly web interface. *BMC Bioinformatics, 25*, 378. https://doi.org/10.1186/s12859-024-06004-0

Cox, J. H., Dietzschold, B., & Schneider, L. G. (1977). Rabies virus glycoprotein. II. Biological and serological characterization. *Infection and Immunity, 16*(3), 754–759. https://doi.org/10.1128/iai.16.3.754-759.1977

Fooks, A. R., Cliquet, F., Finke, S., Freuling, C., Hemachudha, T., Mani, R. S., Müller, T., Nadin-Davis, S., Picard-Meyer, E., Wilde, H., & Banyard, A. C. (2017). Rabies. *Nature Reviews Disease Primers, 3*, 17091. https://doi.org/10.1038/nrdp.2017.91

Giuliani, M. M., Adu-Bobie, J., Comanducci, M., Aricò, B., Savino, S., Santini, L., Brunelli, B., Bambini, S., Biolchi, A., Capecchi, B., Cartocci, E., Ciucchi, L., Di Marcello, F., Ferlicca, F., Galli, B., Luzzi, E., Masignani, V., Serruto, D., Veggi, D., … Pizza, M. (2006). A universal vaccine for serogroup B meningococcus. *Proceedings of the National Academy of Sciences, 103*(29), 10834–10839. https://doi.org/10.1073/pnas.0603940103

Hampson, K., Coudeville, L., Lembo, T., Sambo, M., Kieffer, A., Attlan, M., Barrat, J., Blanton, J. D., Briggs, D. J., Cleaveland, S., Costa, P., Freuling, C. M., Hiby, E., Knopf, L., & Leanes, F. (2015). Estimating the global burden of endemic canine rabies. *PLoS Neglected Tropical Diseases, 9*(4), e0003709. https://doi.org/10.1371/journal.pntd.0003709

Krogh, A., Larsson, B., von Heijne, G., & Sonnhammer, E. L. L. (2001). Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. *Journal of Molecular Biology, 305*(3), 567–580. https://doi.org/10.1006/jmbi.2000.4315

Larsen, M. V., Lundegaard, C., Lund, O., & Nielsen, M. (2005). An integrative approach to CTL epitope prediction: a combined algorithm integrating MHC class I binding, TAP transport efficiency, and proteasomal cleavage predictions. *European Journal of Immunology, 35*(8), 2295–2303. https://doi.org/10.1002/eji.200425811

Liu, B., Zheng, D., Zhou, S., Chen, L., & Yang, J. (2022). VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Research, 50*(D1), D912–D917. https://doi.org/10.1093/nar/gkab1107

Miles, S., Menafra, G., Iriarte, A., & Chabalgoity, J. A. (2025). IApred: A versatile open-source tool for predicting protein antigenicity across diverse pathogens. *ImmunoInformatics*. https://doi.org/10.1016/j.immuno.2025.100014

Ong, E., Wang, H., Wong, M. U., Seetharaman, M., Valdez, N., & He, Y. (2020). Vaxign-ML: supervised machine learning reverse vaccinology model for improved prediction of bacterial protective antigens. *Bioinformatics, 36*(10), 3185–3191. https://doi.org/10.1093/bioinformatics/btaa119

Pizza, M., Scarlato, V., Masignani, V., Giuliani, M. M., Aricò, B., Comanducci, M., Jennings, G. T., Baldi, L., Bartolini, E., Capecchi, B., Galeotti, C. L., Luzzi, E., Manetti, R., Marchetti, E., Mora, M., Nuti, S., Ratti, G., Santini, L., Savino, S., … Rappuoli, R. (2000). Identification of vaccine candidates against serogroup B meningococcus by whole-genome sequencing. *Science, 287*(5459), 1816–1820. https://doi.org/10.1126/science.287.5459.1816

Reynisson, B., Alvarez, B., Paul, S., Peters, B., & Nielsen, M. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research, 48*(W1), W449–W454. https://doi.org/10.1093/nar/gkaa379

Thumuluri, V., Almagro Armenteros, J. J., Johansen, A. R., Nielsen, H., & Winther, O. (2022). DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. *Nucleic Acids Research, 50*(W1), W228–W234. https://doi.org/10.1093/nar/gkac278

Yu, N. Y., Wagner, J. R., Laird, M. R., Melli, G., Rey, S., Lo, R., Dao, P., Sahinalp, S. C., Ester, M., Foster, L. J., & Brinkman, F. S. L. (2010). PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. *Bioinformatics, 26*(13), 1608–1615. https://doi.org/10.1093/bioinformatics/btq249
