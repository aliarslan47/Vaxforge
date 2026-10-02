**VaxForge: konak-özgü epitop çözünürlüğü ile uçtan uca ters aşı bilişimi**

Ali Arslan

Department of Molecular Biology and Genetics, Institute of Natural Sciences, Yıldız Technical University, Istanbul 34220, Turkey

Yazışma: aliarslan47@gmail.com

Taslak: 2026-09-20. Yöntem anlatısı `MANUSCRIPT_METHODS.md`'de; ground-truth antijen kanıtı ve DOI'ler `PROTECTIVE_ANTIGENS.md`'de; aksesyonlar `PATHOGEN_ACCESSIONS.md`'de. Bu belge, İngilizce sürümün (`MANUSCRIPT.md`) Türkçe karşılığıdır.

---

**Özet**

Ters aşı bilişimi araçları bir genomdan aday antijen sıralar; ancak çoğu bu sıralamada durur ve konağı tek bir referans sayar. Bu iki eksik, hem bilinen antijenleri geri getirme başarımını hem de veteriner ve türler-arası hedeflerdeki kullanışlılığı düşürür. Burada VaxForge'u tanıtıyoruz: bir patojeni aday keşfinden yumuşak-filtreli antijen hunisine, ardından dört konak türü için ayrı ayrı sınıf I ve sınıf II epitop tahminine ve karakterize edilmiş çok-epitoplu bir konstrukta taşıyan, in silico immünizasyonla biten bir hat. VaxForge'u iki yerleşik araçla, NERVE 2.0 ve Vaxign-ML ile, varsayılan parametreler ve özdeş girdi altında karşılaştırdık. Meningokok örneğinde VaxForge dört ruhsatlı 4CMenB antijeninin dördünü de geri getirerek Vaxign-ML'i yakaladı; NERVE 2.0 yalnız birini getirdi. Tam proteom ölçeğinde koşulan dokuz zoonotik bakteride daha, genom boyutları 2,0-5,5 Mb arası ve kırk dört koruyucu antijenle, VaxForge 33'ünü geri getirdi; yalnız-sekans sınıflandırıcı Vaxign-ML'in (43/44) altında ama NERVE 2.0'ın (13/44) belirgin önünde, üstelik koruyucu antijenleri rekabetçi ya da daha üst sıralara koyarak (örneğin demir-düzenlenen yüzey belirleyicisi B, 2.889 protein içinde 4. sırada). VaxForge, koruyucu antijenleri salgılanan ya da yüzeyde sunulan türlerde Vaxign-ML'i yakaladı ve yalnızca antijenleri hücre yüzeyinde "moonlighting" yapan sitoplazmik proteinler olan iki türde geride kaldı; bunları lokalizasyon filtresi tasarım gereği eler, yani bir yapaylık değil yorumlanabilir bir ödünleşmedir. NERVE 2.0 daha yüksek fold-zenginleşmeye ulaştı ama bunu yüzey antijenlerini sert seçim basamağıyla eleyerek yaptı ve dokuz türün dördünde hiçbirini döndürmedi; dolayısıyla fold-zenginleşme tek başına yanıltıcıdır. Yedi aileye yayılan on virüste VaxForge, koruyucu ya da patentli antijeni sekizinde adayların üst çeyreğinde geri getirdi ve her durumda türden türe değişen konak-özgü epitop repertuarları verdi (ikili Jaccard örtüşmesi birin altında, bazı konak çiftlerinde sıfır); bu, karşılaştırılan araçların ulaşamadığı bir boyuttur. Antijen sıralamasında Vaxign-ML ile baş-başa sonuç dengeliydi (dörder galibiyet, iki beraberlik; NERVE 2.0 virüsleri işleyemiyor). Kısacası VaxForge, güçlü bir bağımsız tahmincinin recall'ını yakalarken alanın bugün sunmadığı konak çözünürlüğünü ve uçtan uca bir konstruktu ekler.

---

**Giriş**

Ters aşı bilişimi serogrup B *Neisseria meningitidis* ile doğdu. Pizza ve arkadaşları (2000), aday antijenleri kültür bileşenlerinden değil doğrudan tüm-genom dizisinden çıkararak bu yöntemi kurdu ve ruhsatlı dört-bileşenli meningokok aşısının antijenlerini, yani rekombinant faktör H bağlayıcı protein, Neisseria adhezini A ve Neisseria Heparin Bağlayıcı Antijen ile PorA taşıyan dış-membran vezikülünü ortaya çıkardı (Giuliani vd., 2006). Aradan geçen yirmi yılda yöntem olgun yazılımlara kavuştu. NERVE 2.0 güncel ve yerel kurulabilen bir ters aşı bilişimi ortamı sunar (Conte vd., 2024); Vaxign-ML ise bakteriyel koruyucu antijenleri denetimli makine öğrenmesiyle tahmin eder (Ong vd., 2020).

Bu araçların pratikte önemli üç ortak eksiği var. Birincisi, sert bir virülans ya da lokalizasyon filtresi aday listesinin görünürdeki zenginliğini artırırken gerçek antijenleri sessizce eleyebilir. Örneğin Neisseria Heparin Bağlayıcı Antijen küratörlü virülans veritabanlarında yoktur, ama ruhsatlı bir aşı bileşenidir; virülans-veritabanı isabeti şart koşan her hat onu kaybeder. İkincisi, epitop çıktısı çoğunlukla tek bir insan referansına göre hesaplanır. Oysa sığır, domuz ya da farenin asıl konak olduğu pek çok veteriner ve zoonotik hedef, dürüst bir sunum modeli olmadan kalır. Üçüncüsü, araçların çoğu sıralı bir protein listesi verip durur; bir antijeni immünojene dönüştüren adımlar, yani epitop seçimi, konstrukt kurulumu ve ilk immünolojik kontrol, kullanıcıya bırakılır.

VaxForge bu üç eksiği tek bir deterministik hatta kapatır. Sert virülans filtresini yumuşak ve ağırlıklı bir filtreyle değiştirir; böylece virülans-veritabanı kaydı olmayan gerçek bir yüzey antijeni de geçebilir. Kullanıcının seçtiği konak için sınıf I ve sınıf II bağlanmayı, insan, fare, sığır ve domuz allel panellerinden tahmin eder; böylece epitop repertuarı hedeflenen türü yansıtır. Ve sıralamayla yetinmeyip çok-epitoplu bir konstrukta ve in silico immünizasyon çıktısına uzanır. Aşağıda ilk iki katkıyı on bakteri ve on virüste NERVE 2.0 ve Vaxign-ML'e karşı ölçüyor, ayrıca çalışma sırasında yakalayıp giderdiğimiz üç ölçek kaynaklı kusurun ve bir virüse-özgü antijen-filtresi kusurunun oluşturduğu hat doğrulamasını sunuyoruz.

---

**Materyal ve Metot**

**Aday keşfi ve antijen hunisi.** VaxForge deterministiktir ve her kesim değerini tek bir yapılandırma dosyasından okur; hiçbir eşik kodda gömülü değildir (`config/thresholds.yaml`). Girdi türü (proteom, kodlayan diziler ya da genom) otomatik tanındı. Proteinler, küratörlü virülans faktörü veritabanı VFDB'ye (Liu vd., 2022) DIAMOND blastp (Buchfink vd., 2021) ile hizalandı; eşikler e-value en çok 1 × 10⁻⁵, dizi kimliği en az %30 ve kapsam en az %50 olarak alındı. Bu kanıt adaylık puanına katkı verdi, kapı olarak kullanılmadı. Hücre-altı yerleşim PSORTb 3.0 (Yu vd., 2010), yedek olarak DeepLoc 2.0 (Thumuluri vd., 2022) ile tahmin edildi; hücre-dışı, dış-membran, hücre-duvarı ve periplazmik proteinler tutuldu, sitoplazmik olanlar elendi (sert filtre). Transmembran topoloji TMHMM-2.0 (Krogh vd., 2001) ile tahmin edildi; birden çok helikse sahip proteinler elendi (sert filtre). İnsan homolojisi, insan proteomuna karşı DIAMOND blastp ile değerlendirildi; %35'in üzerinde kimlik gösteren proteinler otoimmün güvenlik için çıkarıldı (sert filtre). Sinyal peptidi (SignalP-5.0; Almagro Armenteros vd., 2019, olasılık ≥ 0,5), antijenite (IApred; Miles vd., 2025), suşlar arası korunmuşluk (çok suş verildiğinde ≥ %80) ve adezin olasılığı (≥ 0,51) puana katkı veren yumuşak terimlerdi. Geçen proteinler ağırlıklı bir adaylık puanıyla sıralandı.

**Konak-özgü epitop tahmini.** Seçilen konak için sınıf I bağlanma, 8-11 kalıntılık peptitler için NetMHCpan-4.1 (Reynisson vd., 2020) ile tahmin edildi; güçlü ve zayıf bağlayıcılar sırasıyla ≤ %0,5 ve ≤ %2 yüzdelik-sırada tanımlandı. Antijen işleme, NetCTL-1.2'den allel-bağımsız proteozomal kesim ve TAP taşınımı olarak eklendi (Larsen vd., 2005); sınıf I immünojenite IEDB modelinden alındı (Calis vd., 2013). Sınıf II bağlanma, 9-kalıntılık çekirdekli 15-mer'ler için NetMHCIIpan (Reynisson vd., 2020) ile ≤ %2 ve ≤ %10 yüzdelik-sırada tahmin edildi; yardımcı epitoplar ayrıca insan ve farede interferon-γ indükleme açısından IFNepitope2 (skor ≥ 0,49; Dhall vd., 2024) ile tarandı. Lineer B-hücre epitopları BepiPred-3.0 (Clifford vd., 2022) ile tahmin edildi.

**Epitop güvenliği ve popülasyon kapsamı.** Tutulan epitopların, bilinen alerjenlere altı-kalıntılık birebir benzerlikle değerlendirilerek (FAO/WHO, 2001) alerjen olmaması ve makine-öğrenmesi skoru 0,6'nın üzerindeki peptitlerin çıkarıldığı ToxinPred2 ile değerlendirilerek (Sharma vd., 2022) toksik olmaması istendi; ikisi de sert filtredir. Epitopun konak proteomuna kendi-mimikrisi, birebir 9-kalıntılık eşleşme olarak tanımlanıp yumuşak biçimde işaretlendi. Seçilen epitop setinin insan HLA allellerindeki popülasyon kapsamı IEDB popülasyon-kapsamı aracıyla (Bui vd., 2006) hesaplandı; hedef en az %90 olarak alındı.

**Çok-epitoplu konstrukt ve in silico immünizasyon.** Tutulan epitoplar tek bir konstruktta birleştirildi: sitotoksik epitoplar arasında AAY ve B-hücre epitopları çevresinde KK linkerleri (Schubert & Kohlbacher, 2016), yardımcı epitoplar arasında GPGPG linkerleri (Livingston vd., 2002) ve N-ucunda β-defensin adjuvanı (Biragyn vd., 2002) kullanıldı. Konstruktun fizikokimyasal özellikleri ProtParam (Gasteiger vd., 2005) ile hesaplandı; konstrukt, antikor ve izotip yanıtını zaman içinde benzeten in silico immünizasyon modülüne verildi.

**Benchmark patojenleri ve ground-truth antijenleri.** Bakteri kolunda serogrup B *N. meningitidis* MC58 (UniProt proteom UP000000425, 2.003 protein; RefSeq NC_003112.2), ground-truth olarak dört 4CMenB antijeniyle, fHbp (Q9JXV4), NadA (Q9JXK7), NHBA (Q7DD37) ve PorA (P0DH58), kullanıldı. Buna tam proteom ölçeğinde dokuz tür daha eklendi. İlk üçü *S. aureus* NCTC 8325 (UP000008816, 2.889 protein), *L. monocytogenes* EGD-e (UP000000817, 2.844) ve *S.* Typhimurium LT2 (UP000001014, 4.533); her biri beşer koruyucu ya da aday antijenle. Kalan altısı, dört memeli konağının da enfekte edebildiği, genom boyutu ve hücre biyolojisini kapsayacak biçimde seçilen zoonotik bakterilerdi: *Coxiella burnetii* RSA 493 (UP000002671, 1.812 protein; zorunlu hücre-içi), *Streptococcus suis* serotip 2 BM407 (UP000009077, 1.922), *Brucella suis* 1330 (UP000007104, 3.249), *Mycobacterium bovis* AF2122/97 (UP000001419, 4.088), *Leptospira interrogans* serovar Lai 56601 (UP000001408, 3.676; bir spiroket) ve *Bacillus anthracis* Ames Ancestor (UP000000594, 5.493 protein; kromozom artı toksin ve kapsül antijenlerini taşıyan pXO1 ve pXO2 plazmitleri), her biri dört ya da beş koruyucu antijenle. Dokuz türün koruyucu antijenleri, kanıtı ve DOI'leriyle `PROTECTIVE_ANTIGENS.md`'de, aksesyonları `PATHOGEN_ACCESSIONS.md`'de listelidir. Virüs kolunda, dört memeli konağın da doğal olarak enfekte ettiği, yedi aileye yayılan on virüs kullanıldı: Rhabdoviridae'den kuduz virüsü (*Lyssavirus rabies*; NC_001542, glikoprotein G, NP_056796) ve veziküler stomatit Indiana virüsü (*Vesiculovirus indiana*; NC_001560, glikoprotein G); Picornaviridae'den ensefalomiyokardit virüsü (*Cardiovirus A*; NC_001479, kapsid VP1); Togaviridae'den Venezuela at ensefaliti virüsü (NC_001449, E2) ve Getah virüsü (NC_006558, E2 taşıyan yapısal poliprotein); Flaviviridae'den Japon ensefaliti virüsü (NC_001437, zarf E); Phenuiviridae'den Rift Vadisi humması virüsü (NC_014397, glikoprotein Gn/Gc); Reoviridae'den memeli ortoreovirüsü (σ1 tutunma proteini); ve Poxviridae'den vaccinia (NC_006998, L1R) ve cowpox (NC_003663, L1R ortoloğu) virüsleri. Her ground-truth antijen koruyucu, ruhsatlı ya da patentli bir aşı hedefi taşır; patentler ve kanıt `PATENTED_VACCINES.md`'de listelidir. Poliprotein ve segmentli virüslerde olgun peptitler anote olduğu yerde çıkarıldı; iki ortopoksvirüste tek bir RefSeq referansı kullanıldı ve cowpox L1R'si vaccinia L1R'sine ortoloji ile belirlendi. Tüm aksesyonlar NCBI ve UniProt'a karşı doğrulandı (`PATHOGEN_ACCESSIONS.md`). Meningokokta, tam-proteom VaxForge koşuları CPU-yalnız donanımda pratik olmadığından proteomun sabit-seed alt örneklemi (64 ve 204 protein; dört koruyucu antijen artı arka plan) kullanıldı ve aynı alt örneklem her araca verildi.

**Karşılaştırma araçları ve rolleri.** VaxForge, bakteriyel antijen tahmininde güncel, yerel kurulabilir ve en yaygın kullanılan iki ters aşı bilişimi aracıyla karşılaştırıldı; böylece karşılaştırma bir kullanıcının gerçekte koşacağı araçları yansıtır (Tablo 1). Vaxign-ML (Ong vd., 2020), bir protein proteomunu girdi alıp her protein için bir korujenite (protegenicity) olasılığı döndüren, adayları sıralayan ama epitop, konak bilgisi ya da konstrukt üretmeyen denetimli bir makine-öğrenmesi sınıflandırıcısıdır. NERVE 2.0 (Conte vd., 2024), bir bakteriyel proteomu ve Gram etiketini alıp bir dizi sert seçim filtresinden sonra aday protein kısa listesi döndüren, ardından seçenekli bir epitop-tahmin modülüyle insan HLA süpertipleri için lineer MHC-I/II bağlayıcılarını ve popülasyon kapsamını raporlayan (Conte vd., 2024) uçtan uca bir bakteriyel ters aşı bilişimi ortamıdır. VaxForge ise bakteri ya da virüs için proteom, kodlayan diziler veya genom alır ve sıralı antijenlerle birlikte konak-özgü epitoplar, çok-epitoplu bir konstrukt ve in silico immünizasyon çıktısı döndürür. İki karşılaştırma aracı da yerelde konteynerlerde kuruldu ve özdeş girdi üzerinde, eşik ayarı yapılmadan, varsayılan parametrelerle koşuldu; ground-truth antijenlerdeki recall ve antijen sıraları her aracın kendi çıktısından okundu.

**Tablo 1.** Karşılaştırılan üç aracın girdileri, çıktıları ve yetenekleri (✓ var, ✗ yok). Karşılaştırma araçlarının yetenekleri kendi kaynak yayınlarında tanımlandığı gibi verilmiştir (NERVE 2.0, Conte vd., 2024; Vaxign-ML, Ong vd., 2020).

| Yetenek | VaxForge | NERVE 2.0 | Vaxign-ML |
|---|---|---|---|
| Girdi | proteom, CDS ya da ham genom (otomatik tanıma) | bakteriyel proteom + Gram tipi | bakteriyel protein proteomu |
| Desteklenen patojen | bakteri ve virüs | yalnız bakteri | bakteri (virüs yayında tanımlı değil) |
| Yaklaşım | kural + ML hibrit huni + immünolojik katman | kural tabanlı modüler ortam + AI modülleri | denetimli ML sınıflandırıcı (XGBoost) |
| Antijen sıralama | ✓ ağırlıklı adaylık skoru | ✓ Select + skor | ✓ korujenite (protegenicity) skoru |
| Otoimmünite (konak-homoloji) elemesi | ✓ insan, zorunlu | ✓ insan + fare (seçenekli) | ✗ |
| Virülans-negatif antijeni koruma | ✓ yumuşak huni | ✗ sert Select eler | ✓ skor tabanlı |
| T-hücre / MHC epitop tahmini | ✓ NetMHCpan/II + işleme + immünojenite | ✓ MHC-I/II, lineer, insan HLA süpertipleri | ✗ |
| B-hücre epitop tahmini | ✓ BepiPred-3.0 | ✗ | ✗ |
| Konak-özgü MHC allel seçimi | ✓ insan/fare/sığır/domuz | kısmi (epitop yalnız insan HLA; fare yalnız güvenlik) | ✗ |
| Popülasyon kapsamı | ✓ insan HLA (diğer konak "veri yok") | ✓ insan HLA | ✗ |
| Çok-epitoplu (MEV) konstrukt | ✓ linker + adjuvan | ✗ | ✗ |
| In silico immünizasyon | ✓ | ✗ | ✗ |
| Ana çıktı | sıralı antijen + konak-özgü epitop + konstrukt + immünizasyon | PVC kısa listesi + T-hücre epitop / popülasyon kapsamı | protein başına korujenite skoru |

**İstatistiksel analiz.** Ground-truth antijenlerin adaylar arasındaki zenginleşmesi Fisher kesin testiyle, antijen sıraları ise dağılımdan bağımsız bir permütasyon testiyle sınandı. Ayrım gücü ROC-AUC ve bootstrap güven aralıklarıyla, AUC farkları DeLong testiyle özetlendi; her patojende yalnız dört ila beş pozitif olduğundan AUC küçük-örneklem uyarısıyla ikincil olarak verildi. Konak epitop kümelerinin örtüşmesi Jaccard indeksiyle ölçüldü. Tüm testler `benchmark/stats.py`'de uygulanmıştır.

**Hattın ölçek-güvenliği doğrulaması.** Tam-proteom bakteri koşuları sırasında, VaxForge recall'ını yapay olarak düşüren üç ölçek kaynaklı kusur saptandı ve giderildi. Lokalizasyon adımı başlangıçta tüm proteinleri sabit bir zaman aşımı altında tek konteyner çağrısında işliyordu; en büyük proteomda bu çağrı boş dönüp yüzey antijenlerini yanlış atayan bir sezgisel yedeğe düşürüyordu. Adım artık en çok 800 proteinlik alt-yığınlarda koşuyor. Sınıf I adımı ise NetMHCpan'i tüm proteinlerin birleştirilmiş peptitleriyle tek çağrıda (tam-proteom ölçeğinde yaklaşık 1,37 milyon benzersiz 9-mer) çalıştırıyordu; bu neredeyse sıfır kapsamaya çöküyordu. Adım artık peptitleri 20.000'lik parçalar hâlinde çekirdeklerde paralel işliyor. Üçüncü kusur yalnız en büyük proteomlarda ortaya çıktı: paralel sınıf I işçileri tek bir scratch dizinini paylaşıyordu, dolayısıyla eşzamanlı NetMHCpan örnekleri birbirinde kilitlenip bir koşuyu süresiz durdurabiliyordu; bu, her işçiye izole bir geçici dizin ve parça-başına sınırlı zaman aşımı verilerek çözüldü. Kontrol deneyleri toksisite ve alerjenite filtrelerinin suçsuz olduğunu doğruladı; düzeltmelerin ardından etkilenen türler yeniden koşuldu. Ayrıntı `MANUSCRIPT_METHODS.md`, bölüm B'dedir.

**Veri ve kod erişilebilirliği.** Benchmark betikleri, sonuçları ve figürleri `benchmark/` altındadır; hat yöntem atıfları `vaxforge/citations.py`'de, tüm eşikler `config/thresholds.yaml`'dedir.

---

**Bulgular**

**VaxForge hattı.** VaxForge bir proteomu, kodlayan dizileri ya da ham bir genomu otomatik tanımayla alır; ardından keşif, antijen hunisi, konak-duyarlı epitop tahmini, konstrukt kurulumu ve immünizasyon kontrolü gelir. Keşif adımı proteinleri DIAMOND (Buchfink vd., 2021) ile VFDB'ye (Liu vd., 2022) hizalar. Huni; hücre-altı yerleşimi PSORTb (Yu vd., 2010) ve DeepLoc (Thumuluri vd., 2022), transmembran topolojiyi TMHMM (Krogh vd., 2001), sinyal peptidini SignalP (Almagro Armenteros vd., 2019), antijeniteyi IApred (Miles vd., 2025) ile değerlendirir ve otoimmün güvenlik için zorunlu bir insan-homoloji elemesi uygular; birkaç sert kesim dışında kalanları ağırlıklı bir adaylık puanıyla sıralar.

Huniyi geçen proteinler hattın çekirdeğine, immünolojik katmana girer. NetMHCpan ve NetMHCIIpan (Reynisson vd., 2020) sınıf I ve sınıf II bağlanmayı seçilen konağın allelleri için skorlar; böylece bir T-hücre epitobu ancak hedeflenen türde sunuluyorsa tutulur. T-hücre okumasını üç ek yordayıcı keskinleştirir: sınıf I immünojenite (Calis vd., 2013), proteozomal kesim ile TAP taşınımını birleştiren antijen işleme (Larsen vd., 2005) ve antikor hedefleri için BepiPred-3.0 ile lineer B-hücre epitopları (Clifford vd., 2022). Seçilen epitoplar için VaxForge insan HLA allellerinde popülasyon kapsamını hesaplar; allel frekans verisi bulunmayan konaklarda bu değeri uydurmak yerine "veri yok" olarak bildirir. Tutulan epitoplar, tanımlı linkerlerle (sitotoksik epitoplar arası AAY, yardımcı epitoplar arası GPGPG, B-hücre epitopları için KK) ve bir adjuvanla çok-epitoplu bir konstrukta birleştirilir; konstrukt, antikor ve izotip yanıtını zaman içinde benzeten in silico immünizasyon adımına geçer. Her skorun kaynağı belli olsun diye her yöntem atıfıyla raporlanır.

**VaxForge dört ruhsatlı meningokok antijeninin dördünü de bulur, NERVE 2.0 yalnız birini.** İlk karşılaştırmayı serogrup B *N. meningitidis* (MC58 suşu, proteom UP000000425) üzerinde, dört 4CMenB antijenini ground-truth alarak yaptık. Genişletilmiş 204-protein setinde (dört koruyucu, 200 arka plan, sabit seed) VaxForge dört antijenin dördünü de buldu (recall 4/4; 68 aday; fold 3,12; Fisher p = 0,012); sıralama PorA birinci, faktör H bağlayıcı protein üçüncü, Heparin Bağlayıcı Antijen altıncı, adhezin A sekizinci. Vaxign-ML de dördünü buldu (59 aday); NERVE 2.0 yalnız PorA'yı buldu (recall 1/4). Aynı tablo, kompakt 64-protein setinde ve iki hızlı araç için tam 2.003-protein proteomunda da tekrar etti (Tablo 2; Şekil 1 ve 2). En kritik örnek, virülans-veritabanı negatif olan Heparin Bağlayıcı Antijendir: VaxForge onu yumuşak huniden geçirdi, NERVE 2.0 ise sert seçim adımında eledi. Bağımsız bir makine-öğrenmesi aracının, Vaxign-ML'in, aynı 4/4'e ulaşması, bu davranışın tek araca özgü bir tesadüf değil, yumuşak filtrelemenin bir sonucu olduğunu gösterir.

**Tablo 2.** Üç girdi ölçeğinde meningokok karşılaştırması (varsayılan parametreler, özdeş girdi).

| Set | Araç | Recall | Aday | Fold | Fisher p | AUC |
|---|---|---|---|---|---|---|
| 204 | VaxForge | 4/4 | 68 | 3,12 | 0,012 | 0,99 |
| 204 | Vaxign-ML | 4/4 | 59 | 3,64 | 0,007 | 0,97 |
| 204 | NERVE 2.0 | 1/4 | 3 | 25,0 | 0,058 (ad) | yok |
| Tam 2.003 | Vaxign-ML | 4/4 | 613 | 3,28 | yok | yok |
| Tam 2.003 | NERVE 2.0 | 1/4 | 9 | 62,5 | yok | yok |

![Şekil 1](figures/fig1_antigen_recovery.png)

**Şekil 1.** VaxForge ve Vaxign-ML dört ruhsatlı meningokok antijeninin dördünü de bulur, NERVE 2.0 yalnız PorA'yı bulur. Her nokta, bir aracın bir antijen için sonucudur; renk aracı belirtir (lejant, sağ üst). Noktanın yanındaki sayı, o antijenin aracın aday listesindeki sırasıdır; yani baştan konumu (1 = ilk sırada; küçük sayı = daha iyi). "not recovered" (geri getirilmedi) sütunundaki çarpı, aracın hiç döndürmediği antijeni işaretler. Sıralar ve recall varsayılan parametrelerle özdeş girdi üzerinde hesaplandı.

![Şekil 2](figures/fig2_recall_bar.png)

**Şekil 2.** Bilinen antijenlerdeki recall, araçları AUC'den daha net ayırır. Dört meningokok antijeninin araç-başına recall'ı, ikincil AUC ek notuyla; AUC küçük-örneklem uyarısı taşır ve başlık ölçüt olarak kullanılmaz (n = 4 pozitif).

**Yüksek fold-zenginleşme yüksek recall demek değildir.** NERVE 2.0 kısa aday listeleriyle yüksek fold verdi (204 setinde 25, tam proteomda 62,5); bu, yaklaşık 16,6'lık yayınlanmış değeriyle uyumlu. Ne var ki aday sayısı çok küçük olduğundan Fisher p'si anlamlı değildi ve araç her ölçekte dört ruhsatlı antijenden üçünü kaçırdı (Şekil 3). Bu üç kayıp, faktör H bağlayıcı protein, adhezin A ve Heparin Bağlayıcı Antijen, NERVE 2.0'ın kendi sürekli skorunun düşüklüğünden gelmedi; tersine bu skor onları yüksek verdi. Kayba yol açan, skordan sonra gelen sert seçim adımıydı (periplazmik yerleşim, 0,5 eşiğinin altında adhezin olasılığı ve tek transmembran heliksli bilinmeyen yerleşim). Demek ki fold-zenginleşme seçiciliği ödüllendirir ve recall düşerken bile yükselebilir; bu yüzden yalnız başına okununca aldatıcıdır. Bir aracın gerçekte ne verdiğini gösteren ölçüt, bilinen antijenlerdeki recall ve bu antijenlerin sıralarıdır.

![Şekil 3](figures/fig3_fold_vs_recall.png)

**Şekil 3.** Yüksek fold-zenginleşme yüksek recall demek değildir. Her bar, bir aracın aday listesinde bilinen koruyucu antijenlerin fold-zenginleşmesidir; yani bilinen antijenlerin, adaylar arasında rastlantıya göre kaç kat daha yoğun olduğudur. Barın üstündeki etiket, o koşunun recall'ıdır: aracın toplam içinden gerçekten geri getirdiği bilinen antijen sayısı (4/4 = dördü de bulundu, 1/4 = yalnız biri). Kesikli çizgi yayınlanmış NERVE 2.0 fold değeri 16,6'yı gösterir. NERVE 2.0 en yüksek fold'a en düşük recall'la ulaşır.

**Kalıp dokuz bakteri ve genom ölçeği boyunca sürüyor; klasik-olmayan antijenlerde dürüst bir sınırla.** Tek örnekle sınırlı kalmamak için, dördü de memeli konağı doğal olarak enfekte eden dokuz zoonotik bakteriyi tam proteom ölçeğinde koştuk; genom boyutu 2,0-5,5 Mb arası değişen ve Gram tipi ile hücre biyolojisi dengeli seçilmiş bu set, zorunlu hücre-içi *Coxiella burnetii*'den spiroket *Leptospira interrogans*'a ve *Bacillus anthracis*'in plazmit-kodlu toksinlerine uzanır (Tablo 3). Bu dokuz türün kırk dört koruyucu antijeni arasından VaxForge 33'ünü geri getirdi (Şekil 4); yalnız-sekans sınıflandırıcı Vaxign-ML'in (43/44) altında, NERVE 2.0'ın (13/44) ise belirgin önünde. VaxForge dokuz türün yedisinde 4/5 ya da 5/5'e ulaştı; bunlar arasında *Staphylococcus aureus* (demir-düzenlenen yüzey belirleyicisi B, 2.889 içinde 4. sırada), *Salmonella enterica* serovar Typhimurium (5/5), *Brucella suis*, *Mycobacterium bovis*, *Leptospira interrogans* ve şarbonun koruyucu antijeninin geri getirildiği *Bacillus anthracis* var. NERVE 2.0 sert seçim basamağıyla yüzey antijenlerini yine eledi ve dokuz türün dördünde koruyucu antijenlerin hiçbirini döndürmedi; üstelik *B. suis* ve *B. anthracis*'te herhangi bir aracın en yüksek fold-zenginleşmesine ulaşmasına rağmen.

VaxForge'un geride kaldığı iki tür, *Coxiella burnetii* (2/4) ve *Streptococcus suis* (1/5), bir kusuru değil, gerçek ve yorumlanabilir bir sınırı açığa çıkarır. Orada kaçırdığı antijenler, gliseraldehit-3-fosfat dehidrogenaz, enolaz, 6-fosfoglukonat dehidrogenaz ile DnaJ ve GroEL şaperonları, klasik olarak sitoplazmik ama hücre yüzeyinde "moonlighting" yapan proteinlerdir. VaxForge bunları lokalizasyon sert filtresinde tasarım gereği eler; sekansı doğrudan skorlayan ve lokalizasyon kapısı olmayan Vaxign-ML ise tutar. Dolayısıyla bu davranış bir kusur değil, bilinçli bir ödünleşmedir: lokalizasyon filtresi tam-proteom ölçeğinde sitoplazmik yanlış-pozitifleri bastırır ama klasik-olmayan yüzey antijenlerini de dışlar; moonlighting bir antijeni hedefleyen kullanıcı bu filtreyi gevşetirdi. Her iki aracın tam recall'a ulaştığı türlerde VaxForge koruyucu antijenleri rekabetçi ya da daha üst sıralara koydu; tek kaybı olan *L. monocytogenes* ActA'sı Vaxign-ML listesinde de en dipte çıktı, yani bir hat kusurunu değil, zor bir hedefi yansıtır. Bu tam-proteom koşuları sırasında, ilk denemelerde VaxForge recall'ını yapay olarak düşüren üç ölçek kaynaklı kusuru saptayıp giderdik; lokalizasyon ve sınıf I adımları tek ve aşırı büyük bir yığın altında çöküyordu ve en büyük ölçekte paralel sınıf I işçileri ortak bir scratch dizininde kilitlenebiliyordu, üçü de sınırlı ve izole alt-yığınlara bölündü (`MANUSCRIPT_METHODS.md`, bölüm B). Kontrol deneyleri toksisite ve alerjenite filtrelerinin suçsuz olduğunu doğruladı.

**Tablo 3.** Tam proteom ölçeğinde dokuz-bakteri karşılaştırması (varsayılan parametreler, özdeş girdi), genom boyutuna göre sıralı. Her tür dört ya da beş koruyucu antijen taşır.

| Bakteri (genom, protein) | Araç | Recall | Aday | Fold | AUC |
|---|---|---|---|---|---|
| *C. burnetii* (2,0 Mb, 1.812) | VaxForge | 2/4 | 524 | 1,7 | yok |
|  | Vaxign-ML | 4/4 | 477 | 3,8 | 0,91 |
|  | NERVE 2.0 | 0/4 | 9 | yok | yok |
| *S. suis* (2,2 Mb, 1.922) | VaxForge | 1/5 | 288 | 1,3 | yok |
|  | Vaxign-ML | 5/5 | 628 | 3,1 | 0,88 |
|  | NERVE 2.0 | 0/5 | 0 | yok | yok |
| *S. aureus* (2,8 Mb, 2.889) | VaxForge | 5/5 | 729 | 3,98 | 1,00 |
|  | Vaxign-ML | 5/5 | 838 | 3,46 | 0,94 |
|  | NERVE 2.0 | 3/5 | 34 | 55,8 | yok |
| *L. monocytogenes* (2,9 Mb, 2.844) | VaxForge | 4/5 | 478 | 4,79 | yok |
|  | Vaxign-ML | 5/5 | 917 | 3,11 | 0,95 |
|  | NERVE 2.0 | 3/5 | 73 | 24,3 | yok |
| *B. suis* (3,3 Mb, 3.249) | VaxForge | 4/5 | 929 | 2,8 | yok |
|  | Vaxign-ML | 5/5 | 839 | 3,9 | 0,88 |
|  | NERVE 2.0 | 2/5 | 9 | 185,4 | yok |
| *M. bovis* (4,4 Mb, 4.088) | VaxForge | 4/5 | 535 | 6,2 | yok |
|  | Vaxign-ML | 4/5 | 695 | 4,7 | 0,95 |
|  | NERVE 2.0 | 0/5 | 12 | yok | yok |
| *L. interrogans* (4,7 Mb, 3.676) | VaxForge | 4/5 | 1.105 | 2,7 | yok |
|  | Vaxign-ML | 5/5 | 1.211 | 3,0 | 0,88 |
|  | NERVE 2.0 | 0/5 | 116 | yok | yok |
| *S.* Typhimurium (4,9 Mb, 4.533) | VaxForge | 5/5 | 1.368 | 3,32 | 0,99 |
|  | Vaxign-ML | 5/5 | 1.276 | 3,56 | 0,92 |
|  | NERVE 2.0 | 3/5 | 118 | 23,6 | yok |
| *B. anthracis* (5,5 Mb, 5.493) | VaxForge | 4/5 | 1.163 | 3,8 | yok |
|  | Vaxign-ML | 5/5 | 1.297 | 4,2 | 0,94 |
|  | NERVE 2.0 | 2/5 | 10 | 274,4 | yok |
| **Toplam** | **VaxForge 33/44 · Vaxign-ML 43/44 · NERVE 2.0 13/44** |  |  |  |  |

![Şekil 4](figures/fig4_bacteria_panel.png)

**Şekil 4.** On bakteride tam proteom ölçeğinde recall kalıbı, genom boyutuna göre sıralı. Her aracın bilinen koruyucu antijenlerdeki recall'ı, türden türe antijen sayısı değiştiği için kesir olarak, her barın üstünde geri-getirilen/toplam sayısıyla gösterilir; meningokok örneği dokuz genelleme türünün yanında yer alır ve yıldız bu çalışmada eklenen altı zoonotik bakteriyi işaretler. VaxForge, koruyucu antijenleri salgılanan ya da yüzeyde sunulan türlerde Vaxign-ML'i yakalar ve yalnızca *C. burnetii* ile *S. suis*'te geride kalır; bu ikisinde kaçırılan antijenler, lokalizasyon filtresinin tasarım gereği elediği sitoplazmik moonlighting proteinlerdir. NERVE 2.0 dört türde antijenlerin hiçbirini geri getirmez.

**Tablo S1.** Ground-truth koruyucu antijenler, türe göre; aksesyon, ruhsat durumu ve koruyucu/aday statüsünün doğrulanmış kanıtı (tam kanıt ve DOI'ler `PROTECTIVE_ANTIGENS.md`'de). Durum: *Ruhsatlı* = onaylı insan aşısı bileşeni; *Vet-ruhsatlı suş* = cinsin onaylı canlı/tam-hücre veteriner aşısı var ama buradaki subünit adaydır; *Klinik/Preklinik aday* = deneysel koruyucu ya da immunodominant antijen.

| Tür | Antijen (gen) | Aksesyon | Durum | Anahtar kanıt (DOI) |
|---|---|---|---|---|
| *N. meningitidis* B | fHbp | Q9JXV4 | Ruhsatlı, 4CMenB (Bexsero) | Giuliani 2006, 10.1073/pnas.0603940103 |
| | NadA | Q9JXK7 | Ruhsatlı, 4CMenB | Giuliani 2006, 10.1073/pnas.0603940103 |
| | NHBA | Q7DD37 | Ruhsatlı, 4CMenB | Giuliani 2006, 10.1073/pnas.0603940103 |
| | PorA | P0DH58 | Ruhsatlı, 4CMenB (OMV) | Pizza 2000, 10.1126/science.287.5459.1816 |
| *S. aureus* | IsdB | Q2FZF0 | Klinik aday (V710) | Kuklin 2006, 10.1128/iai.74.4.2215-2223.2006 |
| | IsdA | Q2FZE9 | Preklinik aday | Kim 2010, 10.1016/j.vaccine.2010.02.097 |
| | ClfA | Q2G015 | Preklinik aday | Josefsson 2001, 10.1086/324430 |
| | Hla (α-hemolizin) | Q2G1X0 | Preklinik aday | Adhikari 2012, 10.1371/journal.pone.0038567 |
| | SpA | P02976 | Preklinik aday | Kim 2010, 10.1084/jem.20092514 |
| *L. monocytogenes* | LLO | P13128 | Preklinik aday | Phelps 2020, 10.1016/j.vaccine.2020.06.049 |
| | p60 (Iap) | P21171 | Preklinik (CD8 koruma) | Harty 1995, 10.4049/jimmunol.154.9.4642 |
| | InlA | P0DJM0 | Virülans faktörü / aday | Gaillard 1991, 10.1016/0092-8674(91)90009-n |
| | InlB | P0DQD2 | Virülans faktörü / aday | Dramsi 1995, 10.1111/j.1365-2958.1995.tb02297.x |
| | ActA | P33379 | Pasif-koruma adayı | Asano 2016, 10.1038/srep39628 |
| *S.* Typhimurium | OmpD | P37592 | Preklinik aday | Gil-Cruz 2009, 10.1073/pnas.0812431106 |
| | OmpC | P0A263 | Preklinik aday | Secundino 2006, 10.1111/j.1365-2567.2005.02263.x |
| | OmpA | P02936 | Preklinik aday | Domínguez-Medina 2020, 10.1038/s41467-020-14655-9 |
| | FliC | P06179 | Preklinik aday | Bumann 2014, 10.3389/fimmu.2014.00381 |
| | SseB | Q7BVH7 | Preklinik aday | Lee 2017, 10.4049/jimmunol.1601357 |
| *C. burnetii* | Com1 | H7C7D7 | Preklinik (adoptif transfer) | Wang 2011, 10.1186/1471-2172-12-52 |
| | Mip | P51752 | İmmunodominant aday | Wang 2011, 10.1186/1471-2172-12-52 |
| | OmpH | Q83DT1 | İmmunodominant aday | Wang 2011, 10.1186/1471-2172-12-52 |
| | GroEL | P19421 | İmmunodominant (zayıf) | Wang 2011, 10.1186/1471-2172-12-52 |
| *S. suis* | Suilysin (SLY) | A0A0H3MU02 | Preklinik aday | Du 2013, 10.1016/j.vaccine.2013.07.004 |
| | Enolaz | A0A0H3N5B2 | Preklinik (tam koruma) | Feng 2009, 10.1086/644602 |
| | 6PGD | A0A0H3N193 | Preklinik aday | Zhang 2023, 10.3390/vetsci10060383 |
| | GAPDH | A0A0H3MSV3 | Moonlighting aday | Feng 2009, 10.1086/644602 |
| | DnaJ | A0A0H3MXM0 | Aday | Zhang 2023, 10.3390/vetsci10060383 |
| *B. suis* | Omp25 | Q45689 | Vet-ruhsatlı suş / aday | Cloeckaert 2002, 10.1016/s0378-1135(02)00211-0 |
| | Omp31 | P0A3U5 | Preklinik (koruyucu) | Cloeckaert 2002, 10.1016/s0378-1135(02)00211-0 |
| | BP26 | P0A3U9 | İmmunodominant / tanısal | Kim 2013, 10.1016/j.jmb.2013.01.015 |
| | SodC | P66827 | Preklinik (koruyucu) | Oñate 2003, 10.1128/IAI.71.9.4857-4861.2003 |
| | L7/L12 | P0A469 | Preklinik (koruyucu) | Luo 2006, 10.1128/IAI.74.5.2734-2741.2006 |
| *M. bovis* | ESAT-6 | P0A565 | Klinik aday (BCG-ilişkili) | Brandt 2000, 10.1128/IAI.68.2.791-795.2000 |
| | Ag85B | P0C2T2 | Klinik aday (H1/H56) | Davila 2010, 10.1155/2010/208371 |
| | CFP-10 | P0A567 | Aday / tanısal | Davila 2010, 10.1155/2010/208371 |
| | MPB64 | P0A5Q5 | İmmunodominant / tanısal | Brandt 2000, 10.1128/IAI.68.2.791-795.2000 |
| | HspX | P0A5B8 | İmmunodominant aday | Davila 2010, 10.1155/2010/208371 |
| *L. interrogans* | LigA | Q8EYU4 | Vet-ruhsatlı suş / aday | Faisal 2007, 10.1016/j.vaccine.2007.10.029 |
| | LigB | Q8EZS3 | Preklinik aday | Faisal 2007, 10.1016/j.vaccine.2007.10.029 |
| | LipL32 | O34094 | Preklinik (koruyucu) | Feng 2009, 10.1590/s0100-879x2009005000013 |
| | LipL41 | Q8F8E1 | Preklinik (koruyucu) | Feng 2009, 10.1590/s0100-879x2009005000013 |
| | OmpL1 | G1UB30 | Preklinik (koruyucu) | Feng 2009, 10.1590/s0100-879x2009005000013 |
| *B. anthracis* | PA (pagA) | P13423 | Ruhsatlı (BioThrax/AVA, AVP) | Rhie 2005, 10.1016/j.femsim.2005.05.009 |
| | LF (lef) | P15917 | Preklinik (koruyucu) | Baillie 2010, 10.1016/j.vaccine.2010.07.075 |
| | EF (cya) | P40136 | Toksin bileşeni / aday | Baillie 2010, 10.1016/j.vaccine.2010.07.075 |
| | CapB (γ-PGA) | P19580 | Preklinik (koruyucu ek) | Schneerson 2003, 10.1073/pnas.1633512100 |
| | EA1 (S-tabaka) | P94217 | İmmunodominant aday | — |

**Konak seçimi öngörülen epitop repertuarını değiştirir; karşılaştırılan araçların ulaşamadığı boyut.** VaxForge bağlanma tahminlerini dört konak için küratörlü allel panellerinden alır: insan, fare, sığır ve domuz (Tablo 4). İnsan panelleri IMGT/HLA'dan, diğerleri IPD-MHC'den derlenmiştir. İnsan, fare ve sığırda hem sınıf I hem sınıf II gerçek NetMHCpan ve NetMHCIIpan modelleriyle koşar; sığırın sınıf II BoLA-DRB3 paneli de gerçektir. Domuzda sınıf I gerçek, sınıf II ise tek doğrulanmış allele dayanır. Paneller tam veritabanları değil temsili alt kümelerdir ve geçersiz allel adları çalışma-zamanında elenir; dolayısıyla bildirilen epitoplar temkinlidir.

Virüs kolunu bu konak-özgü sunumu sınamak için kurduk ve dört memeli konağın da doğal olarak enfekte ettiği, yedi aileye yayılan (Rhabdoviridae, Togaviridae, Flaviviridae, Phenuiviridae, Picornaviridae, Reoviridae, Poxviridae) on virüse genişlettik. Her virüste ground-truth, koruyucu, patentli ya da ruhsatlı aşı hedefini taşıyan proteindi (örneğin kuduz ve veziküler-stomatit glikoprotein G'si, alfavirüs E2'si, flavivirüs zarf E'si, bunyavirüs Gn/Gc'si, ortopoksvirüs L1R'si ve pikornavirüs kapsid VP1'i; patentler ve kanıt `PATENTED_VACCINES.md`'de). Viral zarf ve kapsid antijenleri aşının asıl hedefidir ama membrana gömülüdür; bakteriye göre ayarlanmış bir transmembran-sayısı filtresi onları eler. Bu yüzden söz konusu filtreyi virüs profilinde yumuşattık ve bir yüzey-maruziyet sinyali ekledik (Yöntemler). Bununla VaxForge, on virüsün sekizinde koruyucu antijeni adayların üst çeyreğinde geri getirdi; Getah, Rift Vadisi humması ve Japon ensefalitinde birinci, kuduz, veziküler stomatit ve Venezuela at ensefalitinde ikinci sıraya koydu (Tablo 5; Şekil 6). Vaxign-ML ile baş-başa sonuç karışıktı ve bunu dürüstçe bildiriyoruz: VaxForge dört virüste hedefi daha üste, Vaxign-ML dört virüste daha üste koydu, ikisi iki virüste berabere kaldı; yani hiçbir genel sıralayıcı viral antijenlerde baskın değil. VaxForge'un tek açık zayıflığı olan reovirüs σ1'in 11 içinde 9. sırası gizlenmedi, açıklandı: σ1 membran çapası olmayan bir fiber proteinidir, bu yüzden yüzey-maruziyet kredisi almaz; Vaxign-ML onu birinci koydu. Bu, viral sinyali membran topolojisine bağlamanın gerçek bir sınırıdır.

VaxForge'u her iki karşılaştırmalı araçtan ayıran özellik ise sıralama değil, konak boyutudur ve bu on virüsün hepsinde geçerliydi. Her virüste ground-truth-antijen epitop seti konak türüne göre farklıydı; ikili konak Jaccard örtüşmesi her durumda birin altında, birkaç virüste bazı konak çiftleri için hiç ortak epitop yok (Tablo 5; Şekil 7). NERVE 2.0 virüslerde koşamaz, Vaxign-ML tek bir konaksız protegenicity skoru döndürür (kuduz glikoprotein G'yi konak ayrımı olmadan 5 içinde 4. sıraya koydu); dolayısıyla bu konak-başına repertuar üç araç arasında yalnız VaxForge'da vardır. Her virüs tek koruyucu antijen sunduğundan sıra istatistiği çıkarımsal değil betimseldir; popülasyon kapsamı da insan-dışı konaklarda uydurulmadan "veri yok" olarak bildirildi.

**Tablo 4.** Epitop tahmininde kullanılan konak allel panelleri (temsili alt kümeler; tam paneller `config/hosts.yaml`).

| Konak | Sınıf I paneli (n) | Sınıf II paneli (n) | Yordayıcı durumu |
|---|---|---|---|
| İnsan (HLA) | HLA-A ve HLA-B (9) | HLA-DRB1 (6) | NetMHCpan ve NetMHCIIpan, gerçek |
| Fare (H-2) | H-2 K/D (5) | H-2 IA/IE (5) | gerçek |
| Sığır (BoLA) | BoLA sınıf I (8) | BoLA-DRB3 (8) | gerçek |
| Domuz (SLA) | SLA sınıf I (6) | SLA-DRB1 (1) | sınıf I gerçek; sınıf II tek allel |

**Tablo 5.** On dört-konaklı virüste konak-özgü epitop çözünürlüğü ve VaxForge–Vaxign-ML baş-başa karşılaştırması. Sıra, ground-truth antijenin her aracın aday listesindeki konumudur; NERVE 2.0 viral proteomları işleyemediği için gösterilmemiştir. Jaccard, ground-truth-antijen epitoplarının konak-çifti örtüşmesinin ortalaması (aralığı); birin altında olması konak seçiminin repertuarı değiştirdiğini gösterir.

| Virüs (tür) | Ground-truth antijen | VaxForge sırası | Vaxign-ML sırası | Konak-çifti Jaccard, ort. (aralık) |
|---|---|---|---|---|
| Kuduz virüsü (*Lyssavirus rabies*) | glikoprotein G | 2/5 | 4/5 | 0,62 (0,50–0,75) |
| Veziküler stomatit Indiana virüsü (*Vesiculovirus indiana*) | glikoprotein G | 2/5 | 3/5 | 0,58 (0,33–0,78) |
| Ensefalomiyokardit virüsü (*Cardiovirus A*) | kapsid VP1 | 3/11 | 1/11 | 0,44 (0,00–0,78) |
| Rift Vadisi humması virüsü (*Phlebovirus riftense*) | glikoprotein Gn/Gc | 1/4 | 1/4 | 0,33 (0,00–1,00) |
| Venezuela at ensefaliti virüsü (*Alphavirus VEE*) | E2 glikoproteini | 2/9 | 5/10 | 0,37 (0,00–1,00) |
| Getah virüsü (*Alphavirus getah*) | yapısal poliprotein (E2) | 1/3 | 1/3 | 0,50 (0,00–1,00) |
| Japon ensefaliti virüsü (*Orthoflavivirus japonicum*) | zarf proteini E | 1/12 | 3/15 | 0,50 (0,00–1,00) |
| Memeli ortoreovirüsü (*Orthoreovirus mammalis*) | σ1 tutunma proteini | 9/11 | 1/11 | 0,33 (0,00–1,00) |
| Vaccinia virüsü (*Orthopoxvirus vaccinia*) | L1R (IMV zar) | 5/193 | 3/223 | 0,50 (0,00–1,00) |
| Cowpox virüsü (*Orthopoxvirus cowpox*) | L1R ortoloğu (CPXV099) | 8/207 | 3/233 | 0,50 (0,00–1,00) |

![Şekil 6](figures/fig6_virus_antigen_rank.png)

**Şekil 6.** On virüste koruyucu ya da patentli antijenin geri getirilmesi, VaxForge'a karşı Vaxign-ML. Her bar çifti bir virüstür; bar uzunluğu, ground-truth antijenin o aracın aday listesindeki sırasının yüzdelik konumudur (kısa = daha iyi), mutlak sıra etiketlidir. Noktalı çizgi üst çeyreği işaretler. NERVE 2.0 viral proteomları işleyemediği için yoktur. Hiçbir araç baskın değil: VaxForge dört virüste, Vaxign-ML dört virüste önde, ikisi iki virüste berabere.

![Şekil 7](figures/fig7_virus_multihost.png)

**Şekil 7.** Konak seçimi her virüste öngörülen epitop repertuarını yeniden biçimlendirir. Her satır bir virüstür; bar, ground-truth-antijen epitoplarının en düşük–en yüksek ikili konak Jaccard örtüşmesini kapsar, nokta ortalamadır. Bir değeri, iki konağın birebir aynı epitop setini paylaştığı anlamına gelir. Her virüs birin altındadır ve birkaçı en az bir konak çifti için sıfıra iner; yani bir konak için öngörülen epitoplar başka bir konağa aktarılmaz. Bu konak boyutu her iki karşılaştırmalı araçta da yoktur.

![Şekil 5](figures/fig5_multihost_jaccard.png)

**Şekil 5.** Konak seçimi kuduz glikoprotein G epitop repertuarını yeniden biçimlendirir (ikili gösterim). Her bar bir konak çiftidir (örneğin insan–domuz) ve değeri, o iki konak için öngörülen glikoprotein G epitoplarının Jaccard örtüşmesidir: 1 iki konağın birebir aynı epitop setini paylaştığı, 0 ise hiç ortak epitop olmadığı anlamına gelir. Her çift birin altındadır; yani bir konak için öngörülen epitoplar başka bir konağa tam aktarılmaz ve konak seçimi öngörülen repertuarı değiştirir.

---

**Tartışma**

Ters aşı bilişimi araçlarının kıyasları çoğunlukla küratörlü pozitif ve negatif setlerde zenginleşme ya da sınıflandırıcı doğruluğu bildirir (Goodswen vd., 2023). Bizim karşılaştırmamız, bir kullanıcı için daha belirleyici bir ölçüt ekler: her aracın bağımsız olarak doğrulanmış kaç koruyucu antijeni gerçekten döndürdüğü. Ruhsatlı 4CMenB antijenlerinde (Pizza vd., 2000; Giuliani vd., 2006) VaxForge ve Vaxign-ML dördünü de, NERVE 2.0 ise yalnız birini geri getirdi; tam proteom ölçeğinde dokuz bakteri daha koşulduğunda VaxForge 33/44'e, Vaxign-ML 43/44'e, NERVE 2.0 13/44'e ulaştı. Bağımsız makine-öğrenmesi sınıflandırıcısının yüksek recall'ı, Ong vd.'nin (2020) Vaxign-ML için bildirdiği başarımla uyumludur; bu, recall tahminlerimizi tek-araca özgü bırakmak yerine güçlendirir. VaxForge iki türde, *C. burnetii* ve *S. suis*, Vaxign-ML'in gerisinde kaldı; bu türlerin koruyucu antijenleri, lokalizasyon filtresinin elediği sitoplazmik moonlighting proteinlerdir. Bu, lokalizasyon-kapılı bir huninin pratik sınırını ve kullanıcının klasik-olmayan bir yüzey antijeni ararken o kapıyı gevşeteceği noktayı işaret eder. NERVE 2.0'ın eksiği sert seçim basamağından kaynaklanır (Conte vd., 2024): kendi skorunun yüksek verdiği antijenler sabit lokalizasyon ve adezin eşikleriyle elendi; bu, fold-zenginleşmeyi yükseltirken recall'ı düşürdü.

En açık örnek, klasik bir virülans faktörü olmayan ve küratörlü virülans veritabanlarında bulunmayan ruhsatlı bir 4CMenB bileşeni olan Neisseria Heparin Bağlayıcı Antijendir (Giuliani vd., 2006). Virülans-veritabanı isabetine kapı koyan her hat onu kaybeder; VaxForge ise bu kanıtı yumuşak ve puana-ağırlıklı biçimde kullandığı için korur. Bu, ters aşı bilişimi metodoloji literatüründe katı filtrelerin duyarlılığı görünürdeki özgüllük uğruna feda ettiği yönündeki yinelenen uyarıyla uyumludur (Goodswen vd., 2023) ve araçlar karşılaştırılırken tek başına fold-zenginleşme değil, bilinen antijenlerdeki recall'ın raporlanmasını gerektirir.

İkinci ayırt edici özellik, konak-özgü epitop tahmini, bakteri-odaklı araçların açık bıraktığı bir boşluğu kapatır. NetMHCpan ve NetMHCIIpan zaten pek çok insan-dışı alleli destekler (Reynisson vd., 2020); ancak ters aşı bilişimi hatları bunları çoğunlukla tek bir insan referansına uygular. On virüsün hepsi incelenen dört memeli konağı enfekte ettiğinden (Cox vd., 1977; Fooks vd., 2017), her biri için dürüst bir konak-başına repertuar hesaplayabildik; ground-truth-antijen epitop setleri her durumda türden türe farklıydı (ortalama ikili Jaccard 0,33-0,62, bazı konak çiftleri hiç ortak epitop paylaşmıyor). Pratik aşı yükünün çoğunu oluşturan veteriner ve zoonotik hedefler için (Fooks vd., 2017; Hampson vd., 2015), yalnız-insan bir epitop katmanı sunumu yanlış verir; dolayısıyla bu yetenek kozmetik değil, somut bir üstünlüktür. Yalnızca antijen sıralamasında ne VaxForge ne Vaxign-ML on virüs genelinde baskındı; bu, viral proteomlara uygulanan iki genel sıralayıcı için dürüst sonuçtur. VaxForge'un eklediği değer tekdüze daha yüksek bir sıra değil, konak boyutu ve sonraki konstrukttur. VaxForge ayrıca sıralamadan sonra konstrukt kurulumuna ve in silico immünizasyon çıktısına geçer; böylece çıktısı bir liste değil, aday bir immünojendir.

Birkaç sınırlama bu iddiaları çerçeveler. *L. monocytogenes* ActA'sı VaxForge tarafından kaçırıldı ve bağımsız araçça da en dibe kondu; yani bir hat kusuru değil, gerçekten zor bir hedeftir. Geri kazanımı, lineer-epitop bir hattın modellemediği konformasyonel ya da bağlama-özgü sinyalleri gerektirir; zaten B-hücre epitoplarının çoğu süreksizdir (Michalik vd., 2016). Varsayılan toksisite filtresi, listeriolizin O gibi toksin antijenlerinden gelen peptitleri eleyebilir; bu, kıyas için ayarlamadığımız temkinli bir seçimdir. AUC'yi yalnız ikincil ölçüt olarak, açık bir küçük-örneklem uyarısıyla veriyoruz; çünkü patojen başına dört ila beş pozitif AUC farklarını ayrıştıramaz ve bu temelde üstün ayrım iddia etmiyoruz. Son olarak, tam-proteom VaxForge koşuları bugünkü CPU-yalnız donanımda pratik olmadığından meningokok karşılaştırması araçlar arası adil eşlenmiş alt-örneklemler kullandı; tam proteomlar yalnız hızlı araçlarda koşuldu.

İki yön öne çıkıyor. Bir GPU arka ucu, tam-proteom VaxForge koşularına ve her patojende tam ölçekli karşılaştırmaya olanak verir; viral glikoprotein epitoplarının türler-arası ayrışması ise bağlanma ve immünizasyon deneyleriyle sınanabilecek somut bir hipotezdir. İki genişleme de yeniden tasarım gerektirmeden mevcut çerçeveye oturur.

---

**Kaynaklar**

Kaynaklar APA 7 biçiminde ve DOI bağlantılıdır; İngilizce sürümle (`MANUSCRIPT.md`) aynıdır.

Almagro Armenteros, J. J., Tsirigos, K. D., Sønderby, C. K., Petersen, T. N., Winther, O., Brunak, S., von Heijne, G., & Nielsen, H. (2019). SignalP 5.0 improves signal peptide predictions using deep neural networks. *Nature Biotechnology, 37*, 420–423. https://doi.org/10.1038/s41587-019-0036-z

Biragyn, A., Ruffini, P. A., Leifer, C. A., Klyushnenkova, E., Shakhov, A., Chertov, O., Coscia, M., Boyd, L. F., Xu, Y., & Kwak, L. W. (2002). Toll-like receptor 4-dependent activation of dendritic cells by β-defensin 2. *Science, 298*(5595), 1025–1029. https://doi.org/10.1126/science.1075565

Buchfink, B., Reuter, K., & Drost, H. G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods, 18*, 366–368. https://doi.org/10.1038/s41592-021-01101-x

Bui, H. H., Sidney, J., Dinh, K., Southwood, S., Newman, M. J., & Sette, A. (2006). Predicting population coverage of T-cell epitope-based diagnostics and vaccines. *BMC Bioinformatics, 7*, 153. https://doi.org/10.1186/1471-2105-7-153

Calis, J. J. A., Maybeno, M., Greenbaum, J. A., Weiskopf, D., De Silva, A. D., Sette, A., Keşmir, C., & Peters, B. (2013). Properties of MHC class I presented peptides that enhance immunogenicity. *PLoS Computational Biology, 9*(10), e1003266. https://doi.org/10.1371/journal.pcbi.1003266

Clifford, J. N., Høie, M. H., Deleuran, S., Peters, B., Nielsen, M., & Marcatili, P. (2022). BepiPred-3.0: Improved B-cell epitope prediction using protein language models. *Protein Science, 31*(12), e4497. https://doi.org/10.1002/pro.4497

Conte, A., Gulmini, N., Costa, F., Cartura, M., Bröhl, F., Patanè, F., & Filippini, F. (2024). NERVE 2.0: boosting the new enhanced reverse vaccinology environment via artificial intelligence and a user-friendly web interface. *BMC Bioinformatics, 25*, 378. https://doi.org/10.1186/s12859-024-06004-0

Cox, J. H., Dietzschold, B., & Schneider, L. G. (1977). Rabies virus glycoprotein. II. Biological and serological characterization. *Infection and Immunity, 16*(3), 754–759. https://doi.org/10.1128/iai.16.3.754-759.1977

Dhall, A., Patyal, S., & Raghava, G. P. S. (2024). A hybrid method for discovering interferon-gamma inducing peptides in human and mouse. *Scientific Reports, 14*(1), 26859. https://doi.org/10.1038/s41598-024-77957-8

FAO/WHO. (2001). *Evaluation of allergenicity of genetically modified foods. Report of a Joint FAO/WHO Expert Consultation.* https://www.fao.org/3/y0820e/y0820e00.htm

Fooks, A. R., Cliquet, F., Finke, S., Freuling, C., Hemachudha, T., Mani, R. S., Müller, T., Nadin-Davis, S., Picard-Meyer, E., Wilde, H., & Banyard, A. C. (2017). Rabies. *Nature Reviews Disease Primers, 3*, 17091. https://doi.org/10.1038/nrdp.2017.91

Gasteiger, E., Hoogland, C., Gattiker, A., Duvaud, S., Wilkins, M. R., Appel, R. D., & Bairoch, A. (2005). Protein identification and analysis tools on the ExPASy server. In J. M. Walker (Ed.), *The Proteomics Protocols Handbook* (pp. 571–607). Humana Press. https://doi.org/10.1385/1-59259-890-0:571

Giuliani, M. M., Adu-Bobie, J., Comanducci, M., Aricò, B., Savino, S., Santini, L., Brunelli, B., Bambini, S., Biolchi, A., Capecchi, B., Cartocci, E., Ciucchi, L., Di Marcello, F., Ferlicca, F., Galli, B., Luzzi, E., Masignani, V., Serruto, D., Veggi, D., … Pizza, M. (2006). A universal vaccine for serogroup B meningococcus. *Proceedings of the National Academy of Sciences, 103*(29), 10834–10839. https://doi.org/10.1073/pnas.0603940103

Goodswen, S. J., Kennedy, P. J., & Ellis, J. T. (2023). A guide to current methodology and usage of reverse vaccinology towards in silico vaccine discovery. *FEMS Microbiology Reviews, 47*(2), fuad004. https://doi.org/10.1093/femsre/fuad004

Hampson, K., Coudeville, L., Lembo, T., Sambo, M., Kieffer, A., Attlan, M., Barrat, J., Blanton, J. D., Briggs, D. J., Cleaveland, S., Costa, P., Freuling, C. M., Hiby, E., Knopf, L., & Leanes, F. (2015). Estimating the global burden of endemic canine rabies. *PLoS Neglected Tropical Diseases, 9*(4), e0003709. https://doi.org/10.1371/journal.pntd.0003709

Krogh, A., Larsson, B., von Heijne, G., & Sonnhammer, E. L. L. (2001). Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. *Journal of Molecular Biology, 305*(3), 567–580. https://doi.org/10.1006/jmbi.2000.4315

Larsen, M. V., Lundegaard, C., Lund, O., & Nielsen, M. (2005). An integrative approach to CTL epitope prediction: a combined algorithm integrating MHC class I binding, TAP transport efficiency, and proteasomal cleavage predictions. *European Journal of Immunology, 35*(8), 2295–2303. https://doi.org/10.1002/eji.200425811

Liu, B., Zheng, D., Zhou, S., Chen, L., & Yang, J. (2022). VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Research, 50*(D1), D912–D917. https://doi.org/10.1093/nar/gkab1107

Livingston, B., Crimi, C., Newman, M., Higashimoto, Y., Appella, E., Sidney, J., & Sette, A. (2002). A rational strategy to design multiepitope immunogens based on multiple Th lymphocyte epitopes. *Journal of Immunology, 168*(11), 5499–5506. https://doi.org/10.4049/jimmunol.168.11.5499

Michalik, M., Djahanshiri, B., Leo, J. C., & Linke, D. (2016). Reverse vaccinology: the pathway from genomes and epitope predictions to tailored recombinant vaccines. In *Vaccine Design* (Methods in Molecular Biology, Vol. 1403, pp. 87–106). Springer. https://doi.org/10.1007/978-1-4939-3387-7_4

Miles, S., Menafra, G., Iriarte, A., & Chabalgoity, J. A. (2025). IApred: A versatile open-source tool for predicting protein antigenicity across diverse pathogens. *ImmunoInformatics*. https://doi.org/10.1016/j.immuno.2025.100014

Ong, E., Wang, H., Wong, M. U., Seetharaman, M., Valdez, N., & He, Y. (2020). Vaxign-ML: supervised machine learning reverse vaccinology model for improved prediction of bacterial protective antigens. *Bioinformatics, 36*(10), 3185–3191. https://doi.org/10.1093/bioinformatics/btaa119

Pizza, M., Scarlato, V., Masignani, V., Giuliani, M. M., Aricò, B., Comanducci, M., Jennings, G. T., Baldi, L., Bartolini, E., Capecchi, B., Galeotti, C. L., Luzzi, E., Manetti, R., Marchetti, E., Mora, M., Nuti, S., Ratti, G., Santini, L., Savino, S., … Rappuoli, R. (2000). Identification of vaccine candidates against serogroup B meningococcus by whole-genome sequencing. *Science, 287*(5459), 1816–1820. https://doi.org/10.1126/science.287.5459.1816

Reynisson, B., Alvarez, B., Paul, S., Peters, B., & Nielsen, M. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research, 48*(W1), W449–W454. https://doi.org/10.1093/nar/gkaa379

Schubert, B., & Kohlbacher, O. (2016). Designing string-of-beads vaccines with optimal spacers. *Genome Medicine, 8*, 9. https://doi.org/10.1186/s13073-016-0263-6

Sharma, N., Naorem, L. D., Jain, S., & Raghava, G. P. S. (2022). ToxinPred2: an improved method for predicting toxicity of proteins. *Briefings in Bioinformatics, 23*(5), bbac174. https://doi.org/10.1093/bib/bbac174

Thumuluri, V., Almagro Armenteros, J. J., Johansen, A. R., Nielsen, H., & Winther, O. (2022). DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. *Nucleic Acids Research, 50*(W1), W228–W234. https://doi.org/10.1093/nar/gkac278

Yu, N. Y., Wagner, J. R., Laird, M. R., Melli, G., Rey, S., Lo, R., Dao, P., Sahinalp, S. C., Ester, M., Foster, L. J., & Brinkman, F. S. L. (2010). PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. *Bioinformatics, 26*(13), 1608–1615. https://doi.org/10.1093/bioinformatics/btq249
