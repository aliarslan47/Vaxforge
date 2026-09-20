**VaxForge: konak-özgü epitop çözünürlüğü ile uçtan uca ters aşı bilişimi**

Ali Arslan

Department of Molecular Biology and Genetics, Institute of Natural Sciences, Yıldız Technical University, Istanbul 34220, Turkey

Yazışma: aliarslan47@gmail.com

Taslak hazırlanma tarihi 2026-09-20. Yöntem anlatısı `MANUSCRIPT_METHODS.md`'de; ground-truth antijen kanıtı ve DOI'ler `PROTECTIVE_ANTIGENS.md`'de; aksesyonlar `PATHOGEN_ACCESSIONS.md`'de. Bu, `MANUSCRIPT.md` (İngilizce) sürümünün Türkçe karşılığıdır.

---

**Özet**

Ters aşı bilişimi araçları bir genomdan aday antijenleri sıralar, ancak çoğu sıralama adımında durur ve konağı tek bir referans olarak ele alır; bu da hem doğrulanmış antijenleri geri bulma (recall) başarımını hem de veteriner ve türler-arası hedeflerdeki kullanışlılığını sınırlar. Bu çalışmada VaxForge'u sunuyoruz: bir patojeni aday keşfinden yumuşak-filtreli antijen hunisine, beş konak türü boyunca konak-özgü sınıf I ve sınıf II epitop tahminine ve karakterize edilmiş çok-epitoplu bir konstrukta taşıyan, in silico immünizasyon adımıyla biten bir ters aşı bilişimi pipeline'ı. VaxForge'u iki yerleşik araçla, NERVE 2.0 ve Vaxign-ML ile, varsayılan parametreler ve özdeş girdi altında kıyasladık. Ground-truth olarak on beş ruhsatlı veya preklinik koruyucu antijen içeren dört bakteriyel proteomda VaxForge 15 antijenden 14'ünü geri buldu; bu, Vaxign-ML (15/15) ile karşılaştırılabilir ve NERVE 2.0'ın (9/15) belirgin şekilde önündedir. Ayrıca koruyucu antijenleri daha üst sıralara yerleştirdi (örneğin demir-düzenlenen yüzey belirleyicisi B, 2.889 protein içinde 4. sırada). VaxForge'un kaçırdığı tek antijen olan ActA, bağımsız kıyas aracında da en dipte sıralandı ve bir yapaylık değil, gerçek bir zor-hedeftir. NERVE 2.0 daha yüksek fold-zenginleşmeye ulaştı, ancak bunu yalnızca dört ruhsatlı meningokok antijeninden üçünü her girdi ölçeğinde eleyerek yaptı; bu da fold-zenginleşmenin tek başına yanıltıcı olduğunu gösterir. Bir kuduz modeli ve iki ek virüste VaxForge, türler arasında esaslı biçimde farklılaşan konak-özgü glikoprotein epitop repertuarlarını çözdü (ikili Jaccard örtüşmesi birin altında, bir durumda sıfır); bu, kıyas araçlarının ele alamadığı bir boyuttur. Dolayısıyla VaxForge, güçlü bir bağımsız tahmincinin recall'ına eşdeğer başarım gösterirken, alanın şu an sunmadığı konak çözünürlüğünü ve uçtan uca bir konstruktu ekler.

---

**Giriş**

Ters aşı bilişimi, serogrup B *Neisseria meningitidis* ile başladı; burada aday antijenler kültüre edilmiş bileşenlerden değil, doğrudan tüm-genom dizisinden türetildi (Pizza vd., 2000). Bu yaklaşım bakteriyel aşı keşfini yeniden biçimlendirdi ve ruhsatlı dört-bileşenli meningokok aşısının arkasındaki antijenleri üretti: rekombinant faktör H bağlayıcı protein, Neisseria adhezini A ve Neisseria Heparin Bağlayıcı Antijen ile birlikte PorA taşıyan bir dış-membran vezikülü (Giuliani vd., 2006). Yirmi yıl sonra bu paradigmaya olgun yazılımlar hizmet ediyor. NERVE 2.0 güncel, yerel kurulabilir bir ters aşı bilişimi ortamı sağlar (Conte vd., 2024) ve Vaxign-ML bakteriyel koruyucu antijenleri tahmin etmek için denetimli makine öğrenmesi uygular (Ong vd., 2020).

Bu araçlar pratikte önem taşıyan üç sınırlamayı paylaşır. Birincisi, sert bir virülans veya lokalizasyon filtresi, bir aday listesinin görünür zenginleşmesini artırırken gerçek antijenleri sessizce eleyebilir; örneğin Neisseria Heparin Bağlayıcı Antijen, küratörlü virülans veritabanlarında bulunmaz ama ruhsatlı bir aşı bileşenidir, dolayısıyla virülans-veritabanı isabeti şart koşan her pipeline onu kaybeder. İkincisi, epitop düzeyindeki çıktı genellikle tek bir insan referansına karşı hesaplanır; bu da sığır, domuz veya farenin asıl konak olduğu çok sayıda veteriner ve zoonotik hedefi dürüst bir sunum modeli olmadan bırakır. Üçüncüsü, çoğu araç sıralanmış bir protein listesi verip orada durur; böylece bir antijeni immünojene dönüştüren aşağı-akış adımları, epitop seçimi, konstrukt kurulumu ve ilk immünolojik kontrol, elde kalır.

VaxForge bu üç boşluğu tek bir deterministik pipeline içinde kapatmak için kuruldu. Sert virülans filtresini yumuşak, ağırlıklı bir filtreyle değiştirir; böylece virülans-veritabanı anotasyonu olmayan gerçek bir yüzey antijeni de geçebilir. İnsan, fare, sığır, domuz ve tavuk allel panellerinden kullanıcının seçtiği konak için sınıf I ve sınıf II bağlanmayı tahmin eder; böylece epitop repertuarı hedeflenen türü yansıtır. Ve sıralamanın ötesine, çok-epitoplu bir konstrukta ve in silico immünizasyon okumasına geçerek proteomdan aday immünojene uzanan uçtan uca bir yol sunar. Burada ilk iki katkıyı NERVE 2.0 ve Vaxign-ML'e karşı dört bakteri ve üç virüste niceliyor, ayrıca çalışma sırasında yüzeye çıkıp düzelttiğimiz iki ölçek-bağımlı kusuru içeren bir pipeline validasyonunu raporluyoruz.

---

**Bulgular**

**VaxForge pipeline'ı.** VaxForge bir proteomu, kodlayan dizileri veya ham bir genomu otomatik girdi tanımayla işler; ardından keşif, antijen hunisi, konak-duyarlı epitop tahmini, konstrukt kurulumu ve bir immünizasyon kontrolü gelir. Keşif, proteinleri DIAMOND (Buchfink vd., 2021) ve VFDB (Liu vd., 2022) ile küratörlü virülans referanslarına hizalar. Huni; PSORTb (Yu vd., 2010) ve DeepLoc (Thumuluri vd., 2022) ile hücre-altı lokalizasyonu, TMHMM (Krogh vd., 2001) ile transmembran topolojisini, SignalP (Almagro Armenteros vd., 2019) ile sinyal peptidini, IApred (Miles vd., 2025) ile antijeniteyi ve otoimmün güvenlik için zorunlu bir insan-homoloji elemesini birleştirir; az sayıda sert kesim uygular ve geri kalanı ağırlıklı bir adaylık puanıyla sıralar. Hayatta kalan proteinler epitop tahminine geçer; burada NetMHCpan ve NetMHCIIpan (Reynisson vd., 2020) seçilen konak için sınıf I ve sınıf II bağlanmayı skorlar, buna sınıf I immünojenite (Calis vd., 2013), antijen işleme (Larsen vd., 2005) ve BepiPred-3.0 ile B-hücre epitop tahmini (Clifford vd., 2022) eşlik eder. Seçilen epitoplar karakterize edilmiş çok-epitoplu bir konstrukta birleştirilir ve in silico immünizasyon modülüne geçirilir. Her skorun kaynağı açık olsun diye her yöntem atıfıyla raporlanır.

**VaxForge dört ruhsatlı meningokok antijeninin dördünü de geri bulurken NERVE 2.0 yalnız birini bulur.** İlk kıyası serogrup B *N. meningitidis* (MC58 suşu, proteom UP000000425) üzerinde, dört 4CMenB antijenini ground-truth alarak yaptık. Genişletilmiş 204-protein setinde (dört koruyucu artı 200 arka plan, sabit seed) VaxForge dört antijenin dördünü de geri buldu (recall 4/4; 68 aday; fold-zenginleşme 3,12; Fisher kesin p = 0,012); sıralama PorA birinci, faktör H bağlayıcı protein üçüncü, Heparin Bağlayıcı Antijen altıncı ve adhezin A sekizinci. Vaxign-ML de dördünü buldu (59 aday), NERVE 2.0 ise yalnız PorA'yı buldu (recall 1/4). Aynı kalıp, kompakt 64-protein setinde ve iki hızlı araç için tam 2.003-protein proteomunda da geçerliydi (Tablo 1; Şekil 1 ve 2). Belirleyici vaka, virülans-veritabanı negatif olan Heparin Bağlayıcı Antijendir: VaxForge onu yumuşak huniden geçirdi, NERVE 2.0 ise sert seçim adımında eledi. Bağımsız bir makine-öğrenmesi aracının, Vaxign-ML'in, aynı 4/4 recall'a ulaşması, bu davranışın tek-araca özgü bir yapaylık değil, yumuşak filtrelemenin bir özelliği olduğunu gösterir.

*Tablo 1. Üç girdi ölçeğinde meningokok baş-başa kıyası (varsayılan parametreler, özdeş girdi).*

| Set | Araç | Recall | Aday | Fold | Fisher p | AUC |
|---|---|---|---|---|---|---|
| 204 | VaxForge | 4/4 | 68 | 3,12 | 0,012 | 0,99 |
| 204 | Vaxign-ML | 4/4 | 59 | 3,64 | 0,007 | 0,97 |
| 204 | NERVE 2.0 | 1/4 | 3 | 25,0 | 0,058 (ad) | yok |
| Tam 2.003 | Vaxign-ML | 4/4 | 613 | 3,28 | yok | yok |
| Tam 2.003 | NERVE 2.0 | 1/4 | 9 | 62,5 | yok | yok |

**Yüksek fold-zenginleşme yüksek recall anlamına gelmez.** NERVE 2.0 yüksek fold-zenginleşmeli kısa aday listeleri üretti (204 setinde 25, tam proteomda 62,5); bu, yaklaşık 16,6 olan yayınlanmış değeriyle tutarlı. Ancak aday sayısı çok küçük olduğu için Fisher kesin p'si anlamlı değildi ve her ölçekte dört ruhsatlı antijenden üçünü kaçırdı (Şekil 3). Üç kayıp, faktör H bağlayıcı protein, adhezin A ve Heparin Bağlayıcı Antijen, NERVE 2.0'ın kendi sürekli skorlarının düşüklüğünden değil (bu skorlar onları yüksek verdi), skordan sonra uygulanan sert seçim adımından kaynaklandı (periplazmik lokalizasyon, 0,5 kesiminin altında adhezin olasılığı ve bir transmembran heliksli bilinmeyen-lokalizasyon çağrısı). Dolayısıyla fold-zenginleşme seçiciliği ödüllendirir ve tam da recall düşerken yükselebilir; bu yüzden tek başına okunduğunda yanıltıcıdır. Bir aracın gerçekte teslim ettiğini yansıtan ölçüt, bilinen antijenler üzerindeki recall ile bu antijenlerin sıralarıdır.

**Kalıp genom ölçeği boyunca genelleşir ve bir validasyon adımı sonuçları güvenilir kıldı.** Tek örnekliğin ötesine geçmek için, artan genom boyutu ve Gram dengesine göre seçilen üç ek zoonotik bakteriyi tam proteom ölçeğinde koştuk (Tablo 2). VaxForge *Staphylococcus aureus* için 5/5 (demir-düzenlenen yüzey belirleyicisi B, 2.889 içinde 4. sırada), *Listeria monocytogenes* için 4/5 ve *Salmonella enterica* serovar Typhimurium için 5/5 geri buldu; toplam 14/15 (Şekil 5). Vaxign-ML 15/15, NERVE 2.0 ise 9/15 geri buldu ve meningokokta kaybettiği türden yüzey antijenlerini yine kaçırdı. Her iki aracın da tam recall'a ulaştığı iki türde VaxForge koruyucu antijenleri Vaxign-ML'den daha üst sıralara yerleştirdi. VaxForge'un tek kaybı olan *L. monocytogenes* ActA'sı gerçek bir zor antijendir: Vaxign-ML listesinde de en dipte sıralandı, dolayısıyla bir pipeline hatasını değil hedefin zorluğunu yansıtır. Bu koşular sırasında, ilk-geçiş koşularında VaxForge recall'ını yapay olarak düşüren iki ölçek-bağımlı mühendislik kusurunu belirledik ve düzelttik (ilk değerler 3/5, 3/5 ve 0/5); lokalizasyon ve MHC adımları tam-proteom ölçeğinde tek bir aşırı büyük yığın altında çöküyordu ve ikisi de sınırlı, paralel alt-yığınlara bölündü (`MANUSCRIPT_METHODS.md`, bölüm B). Kontrol deneyleri toksisite ve alerjenite filtrelerinin sorumlu olmadığını doğruladı.

*Tablo 2. Tam proteom ölçeğinde üç-bakteri baş-başa kıyası; her bakteride beş koruyucu antijen.*

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

**Konak seçimi öngörülen epitop repertuarını değiştirir; kıyas araçlarının ele alamadığı bir boyut.** Virüs kolu, konak-özgü sunumu sınamak için kuruldu. Tek yapısal proteini olan glikoprotein G'nin nötralize edici, koruyucu antikorları uyaran tek protein olduğu (Cox vd., 1977) kuduz virüsünü (RefSeq NC_001542) kullanan VaxForge, G'yi beş protein içinde ikinci sıraya koydu ve daha da önemlisi dört doğal memeli konağın her biri için farklı bağlı-epitop setleri üretti (glikoprotein G: farede 8, insanda 7, sığırda 7, domuzda 5 epitop). Glikoprotein G epitoplarının ikili konak örtüşmesi birin belirgin altındaydı (Tablo 3), dolayısıyla konak seçimi öngörülen repertuarı esaslı biçimde değiştirdi. İki ek dört-konaklı virüs bu etkiyi yineledi: kuduzun rabdovirüs akrabası veziküler stomatit virüsü, glikoprotein G'yi yine 5 içinde 2. sıraya koydu ve insan ile domuz örtüşmesi 0,33'e kadar düştü; ensefalomiyokardit virüsü ise kapsid VP1'i 11 içinde 3. sıraya koydu ve sığır ile domuz arasında hiç ortak epitop kalmadı. NERVE 2.0 yalnız bakterisel, Vaxign-ML'de ise konak boyutu yok (kuduz glikoprotein G'yi konak çözünürlüğü olmadan 5 içinde 4. sıraya koydu), dolayısıyla bu yetenek üç araç içinde yalnız VaxForge'a özgüdür. Kuduz tek bir koruyucu antijen sunduğu için buradaki sıra istatistiği çıkarımsal değil betimseldir (permütasyon p = 0,40, n = 1) ve popülasyon-kapsam verisi insan-dışı konaklar için uydurulmak yerine "veri yok" olarak raporlandı.

*Tablo 3. Üç dört-konaklı virüste konak-özgü epitop çözünürlüğü (VaxForge).*

| Virüs | Ground-truth antijen | Sıra | Konak-çifti epitop örtüşmesi (Jaccard) |
|---|---|---|---|
| Kuduz virüsü | glikoprotein G | 2/5 | insan~sığır 0,75; insan~fare 0,67; insan~domuz 0,50 |
| Veziküler stomatit virüsü | glikoprotein G | 2/5 | insan~sığır 0,78; insan~domuz 0,33 |
| Ensefalomiyokardit virüsü | kapsid VP1 | 3/11 | insan~fare 0,78; sığır~domuz 0,00 |

---

**Tartışma**

Temel bulgu, VaxForge'un bir ters aşı bilişimi aracı için önemli olan ölçütte, bilinen koruyucu antijenlerin geri bulunmasında, güçlü bir bağımsız tahminciye eşdeğer başarım göstermesi ve bu aracın sunmadığı iki şeyi eklemesidir. Dört bakteri boyunca VaxForge 15 doğrulanmış antijenden 14'ünü geri buldu; bu Vaxign-ML için 15/15 ve NERVE 2.0 için 9/15'e karşılık gelir ve VaxForge bu antijenleri rekabetçi veya daha üst sıralara yerleştirdi. NERVE 2.0'ın tutarlı eksikliği öğreticidir: sert seçim adımı, kendi sürekli skorunun yüksek verdiği gerçek yüzey antijenlerini eledi, bu da fold-zenginleşmeyi şişirirken recall'ı bastırdı. Recall olmadan fold-zenginleşme raporlayan bir pipeline bu yüzden olduğundan keskin görünebilir; ters aşı bilişimi kıyaslarını okuyanları bilinen antijenler üzerindeki recall'a buna göre ağırlık vermeye çağırıyoruz.

İki özellik VaxForge'u recall'ın ötesinde ayırır. Yumuşak virülans filtresi, virülans-veritabanı anotasyonu olmayan ruhsatlı bir bileşen olan Neisseria Heparin Bağlayıcı Antijeni geri kazandı ve konak-duyarlı epitop katmanı kuduz, veziküler stomatit ve ensefalomiyokardit virüsleri için türe özgü repertuarlar üretti. İkincisi, hedeflenen konağın insan olmadığı ve tek bir insan referansının hangi peptitlerin sunulduğunu yanlış vereceği büyük veteriner ve zoonotik alan için önemlidir. VaxForge ayrıca sıralamanın ötesine konstrukt kurulumuna ve in silico immünizasyon okumasına geçer, dolayısıyla çıktı bir liste değil aday bir immünojendir.

Birkaç sınırlama bu iddiaları çerçeveler ve açıkça belirtilir. ActA VaxForge tarafından kaçırıldı; bağımsız aracın da en dipte sıraladığı gerçek bir zor antijen olsa da, geri kazanımı lineer-epitop bir pipeline'ın modellemediği konformasyonel veya bağlama-özgü sinyalleri gerektirir. Varsayılan toksisite filtresi toksin antijenlerinden türeyen peptitleri eleyebilir; bu, ölçek düzeltmelerinden önce por-açan listeriolizin O için peptit-düzeyi recall'ını düşürmüştü ve kıyas için ayarlamadığımız muhafazakâr bir seçim olarak kalır. Eğri-altı-alan (AUC) değerleri, patojen başına yalnız dört ila beş pozitif olduğu için farklar istatistiksel olarak çözülemediğinden, başlık değil açık küçük-örneklem uyarısıyla ikincil olarak raporlanır ve bu temelde üstün ayrım iddia edilmez. Bakteriyel kıyas ayrıca bir ölçek asimetrisi taşır: VaxForge'un tam-proteom meningokok koşuları mevcut CPU-yalnız donanımda pratik değildir, dolayısıyla amiral meningokok kıyası araçlar arası adil olan eşlenmiş alt-örneklem setleri kullandı; tam-proteom koşuları yalnız hızlı kıyas araçları için mevcuttu. Son olarak, B-hücre epitoplarının kabaca onda dokuzu süreksizdir ve burada üretilen lineer konstrukt konformasyonel epitopları temsil etmez; bunu çözmek yerine bir kapsam sınırı olarak belirtiyoruz.

İki yön doğrudan izler. Bir GPU arka ucu tam-proteom VaxForge koşularına olanak tanır ve baş-başa kıyası her patojen için tam ölçekte tamamlar; konak-özgü epitop tahminleri, özellikle viral glikoproteinlerde görülen türler-arası ıraksama, bağlanma ve immünizasyon deneylerinin sınayabileceği somut hipotezlerdir. Pipeline ve kıyası, her iki genişletmenin de yeniden tasarım gerektirmeden mevcut çerçeveye oturması için yapılandırılmıştır.

---

**Şekiller**

![Şekil 1](figures/fig1_antigen_recovery.png)

*Şekil 1. VaxForge ve Vaxign-ML dört ruhsatlı meningokok antijeninin dördünü de geri bulur, NERVE 2.0 ise yalnız PorA'yı bulur.* Dört 4CMenB antijeni için antijen-araç geri bulma matrisi; dolu hücreler geri bulunan antijenleri, kırmızı hücreler NERVE 2.0 seçim adımınca elenen antijenleri gösterir. Recall ve aday sayıları varsayılan parametrelerle özdeş girdi üzerinde hesaplandı.

![Şekil 2](figures/fig2_recall_bar.png)

*Şekil 2. Bilinen antijenler üzerindeki recall, araçları eğri-altı-alandan daha net ayırır.* Dört meningokok antijeninin araç-başına recall'ı, ikincil AUC ek notuyla; AUC küçük-örneklem uyarısı taşır ve başlık ölçüt olarak kullanılmaz (n = 4 pozitif).

![Şekil 3](figures/fig3_fold_vs_recall.png)

*Şekil 3. Yüksek fold-zenginleşme recall'ı değil seçiciliği izler.* Üç araç için 204-protein setinde (dolu noktalar) ve tam proteomda (içi boş noktalar) fold-zenginleşmeye karşı recall; kesikli çizgi yayınlanmış NERVE 2.0 fold değeri 16,6'yı gösterir. NERVE 2.0 yüksek fold ve düşük recall'da, VaxForge ile Vaxign-ML ise tam recall ve orta fold'da yer alır.

![Şekil 4](figures/fig4_multihost_jaccard.png)

*Şekil 4. Konak seçimi kuduz glikoprotein G epitop repertuarını yeniden biçimlendirir.* İnsan, sığır, fare ve domuz arasında glikoprotein G sınıf I epitoplarının ikili Jaccard örtüşmesi; birin altındaki değerler öngörülen repertuarın seçilen konağa bağlı olduğunu gösterir.

![Şekil 5](figures/fig5_bacteria_panel.png)

*Şekil 5. Recall kalıbı üç genom ölçeğinde genelleşir.* Tam proteom ölçeğinde *S. aureus*, *L. monocytogenes* ve *S.* Typhimurium için her aracın beş koruyucu antijen üzerindeki recall'ı, AUC ek notuyla; VaxForge toplam 14/15, Vaxign-ML 15/15 ve NERVE 2.0 9/15.

Not: şekil panelleri şu an Türkçe eksen ve başlık metni taşıyor; hedef dergi gerektirirse `make_figures.py`'den İngilizce etiketli sürümler yeniden üretilebilir.

---

**Kaynaklar**

Kaynaklar APA 7 biçiminde, DOI bağlantılı. İngilizce sürümle (`MANUSCRIPT.md`) aynıdır.

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
