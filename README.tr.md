# Vaxforge

Web arayüzlü, ajan-destekli, in silico **ters aşılama (reverse vaccinology)** pipeline'ı — bir patojen dizisinden sıralanmış çok-epitoplu mRNA aşı yapısına ve tam atıflı bir rapora.

[![type](https://img.shields.io/badge/type-reverse%20vaccinology-0d6b8f)](https://github.com/aliarslan47/Vaxforge)
[![interface](https://img.shields.io/badge/interface-Streamlit%20%2B%20Next.js-2f8f5b)](https://github.com/aliarslan47/Vaxforge)
[![tools](https://img.shields.io/badge/tools-real%20%C2%B7%20cited-c07211)](https://github.com/aliarslan47/Vaxforge)

**Türkçe** · [English](README.md)

## Nedir?

Vaxforge, Forge ailesinin aşı-tasarımı üyesidir — bir patojen dosyasını sıralanmış aşı adaylarına çeviren bir ters aşılama pipeline'ı. Biyolojiyi deterministik bir bilimsel çekirdek yapar; LLM sarmalayıcı yalnızca planlar, yorumlar ve raporlar. İki arayüzle gelir: yerel analiz için tek-dosya Streamlit uygulaması ve çok-koşumlu web dağıtımı için FastAPI + Next.js yığını.

## Ne yapar?

Bir patojen dosyası (FASTA/FASTQ) yükle → sistem girdiyi otomatik algılar → virülans faktörlerini ve aşı hedeflerini tarar → antijen hunisini çalıştırır → seçilebilir konakçılarda B- ve T-hücre epitoplarını tahmin eder → hayatta kalma filtreleri (toksisite, alerjenite) → adayları yapılandırılabilir bir adaylık skoruyla sıralar → çok-epitoplu bir mRNA yapısı kurar → tam araç atıflı, yayın-tarzı bir rapor üretir.

- **Gerçek, yayınlanmış araçlar**; her biri raporda (PDF/HTML/JSON) ve arayüzde atıflı.
- **Konakçı-seçilebilir, çok-organizmalı MHC**: `config/hosts.yaml`'da dört konak etkin (insan, fare, sığır, domuz); tüm aleller yerel NetMHCpan alel listesine karşı valide. Domuz MHC-II bilinçli olarak kapalı, çünkü NetMHCIIpan'da SLA-DR modeli yok; proxy kullanmak sahte epitop üretirdi.
- **Sihirli sayı yok**: eşikler `config/thresholds.yaml`'da (organizma ön-ayarları), arayüzden düzenlenebilir ve her rapora yazılır.
- **Dürüst etiketleme**: gerçek bir araç yoksa açıkça-etiketli bir fallback kullanılır.

## Kurulum

```bash
git clone https://github.com/aliarslan47/Vaxforge.git
cd Vaxforge
python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt
```

Harici araçlar `tools/` altındadır (git-ignore). Ücretsiz araçlar otomatik kurulur (DIAMOND, VFDB, Swiss-Prot, ToxinPred2, IApred); lisanslı DTU HealthTech araçları (NetMHCpan, NetMHCIIpan, SignalP-5.0, TMHMM-2.0, BepiPred-1.0) manuel akademik indirme gerektirir — bkz. `tools/README.md`. Bir araç eksikse pipeline etiketli bir fallback ile yine çalışır.

**Herhangi bir DTU HealthTech paketini çıkardıktan sonra bir kez koşturun:**

```bash
python3 scripts/configure_dtu_wrappers.py          # yolları düzelt
python3 scripts/configure_dtu_wrappers.py --check  # yalnız doğrula (sorun varsa çıkış 1)
python3 scripts/setup_toxinpred.py                 # ToxinPred2: envfile, pandas API, RF modeli
python3 scripts/convert_ifn_models.py              # IFNepitope2 sklearn pickle köprüsü
```

**PSORTb (bakteriyel lokalizasyon) rootless podman konteynerinde koşar.** Pull
`potentially insufficient UIDs or GIDs available in user namespace` hatası
veriyorsa, image'ın istediği UID `/etc/subuid` aralığının dışındadır. `sudo` ya
da `usermod` GEREKMEZ: image'ı bir kez, aralık dışı chown hatalarını yok sayan
overlay sürücü seçeneğiyle çekin; sonrasında sade `podman run` çalışır.

```bash
podman --storage-opt ignore_chown_errors=true pull docker.io/brinkmanlab/psortb_commandline:1.0.2
podman run --rm docker.io/brinkmanlab/psortb_commandline:1.0.2 /usr/local/psortb/bin/psort --version
```

DTU paketleri elle düzenlenmesi gereken satıcı **yer tutucu** yollarıyla gelir
(`setenv NMHOME /tools/src/netMHCpan-4.2`, `setenv BP /usr/cbs/packages/...`).
Düzeltilmezse araç hiç hesap yapmaz; `tools/` git-ignore olduğu için bu
düzeltmeler sürüm kontrolünde saklanmaz — paketi yeniden çıkarmak yer tutucuyu
geri getirir. Pipeline artık yapılandırılmamış sarmalayıcıyı tespit edip aracı
sessizce vekile düşürmek yerine **yok** olarak bildiriyor; ayrıca her koşu
hangi dış araçların gerçekten hazır olduğunu listeleyen bir **Araç durumu**
bölümü kaydediyor.

## Kullanım

Yerel Streamlit uygulaması:

```bash
streamlit run app.py
```

Ya da web yığını, API 8000 portunda, arayüz 3000 portunda:

```bash
cd web/backend && pip install -r requirements.txt && uvicorn main:app --port 8000
cd web/frontend && npm install && npm run dev
```

Bir FASTA/FASTQ dosyası sürükle (ya da örnek seç), patojen profilini ve konakçı(lar)ı seç, eşikleri gözden geçir ve çalıştır. Çıktılar: yayın-tarzı **PDF**, **HTML** panosu, sıralı **CSV**, peptit **FASTA**, **GenBank** mRNA yapısı ve tam-run **JSON**.

## Benchmark

Yayınlanmış iki ters aşılama aracıyla head-to-head; tüm araçlar DEFAULT parametre ve
özdeş girdiyle koştu, eşik tuning yapılmadı. Ground-truth, ruhsatlı ya da literatürle
doğrulanmış koruyucu antijenler. Ayrıntılar, aksesyonlar ve istatistikler `benchmark/`
altında.

On bakteri, tam proteom, 48 koruyucu antijen üzerinden recall:

| Araç | Recall | Oran |
|---|---|---|
| Vaxign-ML | 47/48 | %98 |
| VaxForge | 37/48 (geçici, 2026-10 koşusu) | %77 |
| NERVE 2.0 | 14/48 | %29 |

> **VaxForge satırı eskidir ve HENÜZ DÜZELTİLMEDİ.** Aşağıdaki üç uyarıda sayılan
> düzeltmelerin tamamından ÖNCEye, 2026-10 koşusuna aittir. Elle düzeltilmek yerine
> olduğu gibi bırakıldı: bu sayıyı güncellemenin tek dürüst yolu benchmark'ı yeniden
> koşmaktır, tablodaki hiçbir sayı elle yeniden hesaplanmamıştır. Rakip satırları
> etkilenmez.

Ground-truth'un ruhsatlı dört 4CMenB antijeni olduğu MenB exemplar'ında VaxForge 4/4,
NERVE 2.0 ise 1/4 yakalar.

Tablo 2026-10 koşusunun sonucudur. Aşağıda anlatılan Gram-pozitif lokalizasyon
düzeltmesinden ÖNCEye aittir, dolayısıyla benchmark yeniden koşulduğunda VaxForge
sayısı değişecektir. Rakip sayıları etkilenmez.

**İkinci uyarı, 2026-11'de eklendi.** O koşu, sonradan bulunan bir dizi araç-düzeyi
kusurdan da ÖNCEye aittir: TMHMM, IApred ve ToxinPred2 kendilerini "hazır"
bildirirken hiç çalışmıyordu (bozuk yorumlayıcı yolu, iki eksik Python bağımlılığı,
eski bir model pickle'ı). Bu ortamda sonuç şuydu: transmembran topolojisi ve
antijenite hiç ölçülmedi, toksisite hiçbir şeyi elemedi. 2026-10 koşusunun etkilenip
etkilenmediği **dosyalardan ANLAŞILAMIYOR**, çünkü o dönemin benchmark özetleri hangi
araçların gerçekten koştuğunu kaydetmiyordu. Artık kaydediyor, yani bu soru bundan
sonraki her koşu için yanıtlanabilir. Bu yüzden VaxForge sayısı, benchmark yeniden
koşulana kadar GEÇİCİ sayılmalı; rakipler kendi konteynerlerinde koştuğu için her
hâlükârda etkilenmez.

**Üçüncü uyarı, 2026-11-23'te eklendi.** Aday PUANLAMA yolunun kendisinde de, hepsi
2026-10 koşusundan sonra bulunan bir dizi kusur çıktı: MHC sınıf II tahmini sessizce
yedeğe düşüyordu, bir `%rank` ölçek hatası poli-alanin dizilerinin epitop gibi
puanlanmasına izin veriyordu, konservasyon tek suşa (girdinin kendisine) karşı
hesaplanıyordu — yani hiçbir karşılaştırma yapılmamışken %100 okunabiliyordu — ve
yalnızca bir birim testiyle bulunan sonuncusu: kalıntı-başı konservasyon profili,
yerel hizalama adayın 1. kalıntısından başlamadığı her durumda **KAYIYORDU** (74
gerçek hizalamanın 12'si, 23-100 kalıntı kayma). Bu son kusur `epitope_conservation`'ı,
o da adaylık puanını besliyor; yani protein ortalamaları doğru görünürken sıralamayı
sessizce bozuyordu. 48 antijenlik recall üzerindeki etkinin yönü ve büyüklüğü,
benchmark yeniden koşulana kadar BİLİNMİYOR.

VaxForge'un kaçırdığı on bir antijenin tamamı `benchmark/MISS_ATTRIBUTION.md` içinde
mekanistik olarak atfedilmiştir. Altısı, konakçı homolojisi ve lokalizasyon kapılarının
elediği sitoplazmik moonlighting antijenidir; bu, sitoplazmik yanlış-pozitifi bastıran
bilinçli bir ödünleşmedir. Ölçekleme hatası bulunmamıştır. AUC yalnızca ikincil ve
caveat'lı bir metrik olarak verilir, çünkü pozitif sayısı farkların anlamlı olması için
fazla küçüktür.

### Bilinen kısıt ve o koşudan sonra uygulanan düzeltme

Bakteri izinli lokalizasyon listesi Gram-negatif anatomisiyle yazılmıştı. Gram-pozitif
bakteride dış zar yoktur, dolayısıyla sitoplazmik zar en dış sınırdır ve yüzeye çıkan
proteinler PSORTb'de `membrane` etiketi alır. Tek ortak liste bu yüzden meşru
Gram-pozitif yüzey antijenlerini atıyordu. İzinli liste artık Gram tipine göre
seçiliyor: Gram-pozitif koşu `extracellular`, `cell_wall` ve `membrane` kabul ediyor,
Gram-negatif koşu ise perdelenmiş iç zarı elemeye devam ediyor. Ayrıntı
`config/thresholds.yaml` içindeki notta.

## Modüller

Bilimsel çekirdek gerçek, atıflı araçlarla çalışır. Konstrukt-düzeyi yapı tahmini ESMAtlas API üzerinden etkindir; yalnızca AlphaFold düzeyinde peptit–MHC modellemesi ve moleküler dinamik GPU makinesine ertelenmiştir.

| Adım | Araç / yöntem | Durum |
|---|---|---|
| Girdi otomatik algılama | FASTA/FASTQ · nt/protein · genom/CDS/okuma | yerleşik |
| Keşif (virülans faktörleri) | DIAMOND + VFDB | ✅ gerçek |
| Antijen hunisi — lokalizasyon | PSORTb 3.0 (bakteri) · aksi hâlde Kyte-Doolittle + sinyal heuristiği | ✅ gerçek / etiketli heuristik |
| Antijen hunisi — transmembran | TMHMM-2.0 | ✅ gerçek |
| Antijen hunisi — sinyal peptidi | SignalP-5.0 | ✅ gerçek |
| Antijen hunisi — antijenite | IApred | ✅ gerçek |
| Antijen hunisi — konakçı homolojisi (güvenlik) | DIAMOND vs insan Swiss-Prot | ✅ gerçek |
| Epitoplar — B-hücre | BepiPred-1.0 | ✅ gerçek |
| Epitoplar — MHC-I / MHC-II | NetMHCpan / NetMHCIIpan (yerel + IEDB) | ✅ gerçek |
| Hayatta kalma — toksisite / alerjenite | ToxinPred2 · FAO/WHO 6-mer + UniProt | ✅ gerçek |
| Adaylık skorlaması | ağırlıklı, yapılandırılabilir | ✅ gerçek |
| mRNA yapısı | linker + adjuvan + insan kodon-opt + GC/CAI | ✅ gerçek |
| Antijen hunisi — lokalizasyon (bakteri) | PSORTb 3.0 (podman) | ✅ gerçek |
| Antijen hunisi — antijenite çapraz-kontrol | VaxiJen ACC yöntemi (bağımsız 2. kanıt) | ✅ gerçek |
| Güvenlik — konakçı mimikrisi | tam k-mer + BLOSUM62 konservatif ikame | ✅ gerçek |
| Güvenlik — konakçı DB otomatik sağlama | UniProt Swiss-Prot + DIAMOND makedb | ✅ gerçek |
| Konservasyon — çok-suş | NCBI Datasets otomatik çekme + artık-başına konservasyon | ✅ gerçek |
| Epitoplar — IFN-γ indüksiyonu (HTL) | IFNepitope2 | ✅ gerçek |
| Epitoplar — popülasyon kapsamı | IEDB popülasyon kapsam yöntemi | ✅ gerçek |
| Glikozilasyon — glikan-maskeli epitop | N-X-S/T sequon + O-bölge sezgiseli (yumuşak uyarı) | ✅ gerçek |
| Konstrukt yapısı | ESMFold, ESMAtlas API (≤400 aa) → pLDDT | ✅ gerçek |
| Konstrukt yüzeyi | Shrake-Rupley RSA + konformasyonel B yamaları | ✅ gerçek |
| Üretilebilirlik | E. coli kodon optimizasyonu, CAI, GC%, restriksiyon taraması, pET vektör | ✅ gerçek |
| İmmünizasyon kararı | ImmForge motoru, 9 tür kalibreli | ✅ gerçek |
| HLA–hastalık anotasyonu | küratörlü alel tablosu, satır başına PMID/DOI (yalnız yorumlama) | ✅ gerçek |
| Peptit–MHC yapısı + MD | AlphaFold + docking/MD | ⏸️ ertelendi (GPU) |

Tam atıflar (`vaxforge/citations.py`) ve proje yapısı depoda.

---

Forge ailesi: [RNAForge](https://github.com/aliarslan47/RNAForge) (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bakteri) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virüs/faj) · [MicrobiomeForge](https://github.com/aliarslan47/MicrobiomeForge) (mikrobiyom) · **Vaxforge** (ters aşılama) · [ImmForge](https://github.com/aliarslan47/ImmForge) (bağışıklık simülasyonu) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG üreticisi).
