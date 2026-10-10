# VaxForge — 11 Eksik Kapatma Yol Haritası

Ters-aşı-bilişimi pipeline'ında tespit edilen 11 eksik, önem sırasına göre tek tek
kapatılıyor. Her faz: **pipeline'da çalışır + rapor (HTML/PDF/run.json) + Next.js web +
Streamlit** üçlüsünde görünür, pipeline asla bozulmaz, faz sonunda **commit + smoke testi**.

Plan dosyası: `~/.claude/plans/tender-munching-thacker.md`

## Değişmez kurallar
1. Modül `available()`/`predict()` — araç yoksa dürüst fallback (uydurma yok).
2. Pipeline aşaması `try/except` ile sarılı → hata pipeline'ı düşürmez.
3. Rapor: `report.py` `_xxx_section` + `_html()` montaj.
4. Katalog: `web/backend/runner.py:_PIPELINE_STEPS`.
5. Next.js: `results.tsx` bileşen + `lib/api.ts` tip + `lib/i18n.ts`.
6. Streamlit: `app.py` panel.
7. i18n: `vaxforge/i18n.py` (TR+EN).
8. Faz sonu: smoke geç → commit → bu tabloyu güncelle.

## Durum tablosu

| # | Faz | Eksik | Durum | Commit |
|---|-----|-------|-------|--------|
| 0 | Güvenlik ağı + kayıt (smoke test, ROADMAP) | #10-kısmi | ✅ | f70f4b0 |
| 1 | PSORTb lokalizasyon (bakteri, Gram) | #7 | ✅ | ✓ |
| 2 | DeepLoc2 lokalizasyon (ökaryot/virüs) | #8 | ✅ | ✓ |
| 3 | ESMFold yapı + konformasyonel B + RSA | #1 | ✅ | ✓ |
| 4 | Glikozilasyon + glikan-maskeli epitop | #2 | ✅ | ✓ |
| 5 | Kodon optimizasyonu + in silico klonlama | #3 | ✅ | ✓ |
| 6 | VaxiJen antijenite çapraz-kontrol | #4 | ✅ | ✓ |
| 7 | BLOSUM tabanlı kısmi mimikri | #5 | ✅ | ✓ |
| 8 | Otomatik çok-suş konservasyonu | #6 | ✅ | ✓ |
| 9 | Otoimmünite konak-DB genişletme | #9 | ✅ | ✓ |
| 10 | Türe-özgü ImmForge kalibrasyonu | #11 | ✅ | ✓ |
| 11 | Kapsamlı birim test paketi | #10 | ✅ | ✓ |

Durum: ⬜ bekliyor · 🔨 sürüyor · ✅ bitti

## Kaldığın yer / notlar
- **2026-09-26 Faz 0 ✅**: smoke testi (module-scope fixture, ~4dk) + ROADMAP. Baz commit f70f4b0.
- **2026-09-26 Faz 1 ✅**: PSORTb podman image çekildi (`psortb.available()=True`); funnel
  otomatik GERÇEK'e geçti. Funnel yöntem-özeti artık `meta['funnel_methods']`'e yazılıyor →
  rapor (`_methods_section`) + Next.js (`MethodsUsed` kartı, run+results sayfaları) + Streamlit
  (yöntem tablosu) + i18n TR/EN. Smoke geçti (2 passed, 4:15), tsc temiz.
- **2026-09-26 Faz 2 ✅**: `vaxforge/deeploc.py` (Fast/ESM1b, venv CLI; ESM1b 7.3GB önbellekte).
  funnel: **parazit** → DeepLoc SERT filtre (config allowed), **virüs** → DeepLoc YUMUŞAK
  (virionlar için güvenilmez, sadece etiket). Bakteri PSORTb'de kaldı (DeepLoc ökaryot-yanlı).
  Doğrulandı: parazit koşusu "GERÇEK (DeepLoc-2.0, Fast) [sert:parasite]". Rapor/UI görünürlüğü
  Faz 1'in funnel_methods altyapısından otomatik. runner katalog + EN token (sert/yumuşak) eklendi.
- **2026-09-26 Faz 3 ✅**: yapısal katman — 3 modül: `structure.py` (ESMAtlas fold API,
  ≤400aa, pLDDT 0-100'e normalize), `rsa.py` (Bio.PDB Shrake-Rupley + Tien2013 maxASA →
  RSA, yüzey%), `bcell_conf.py` (RSA+Parker → konformasyonel yüzey B yamaları, DiscoTope-benzeri).
  Pipeline: MEV sonrası "structure" aşaması MEV konstruktu katlar → pLDDT/yüzey%/konf-B.
  Rapor `_structure_section` + Next.js `StructureCard` + Streamlit paneli + runner katalog satırı.
  Ağ yok/>400aa → dürüst ARAÇ_YOK (BepiPred lineer korunur). Chain doğrulandı (pLDDT 42.4).
- **2026-09-27 Faz 4 ✅**: `vaxforge/glyco.py` — N-X-S/T sequon (kesin motif) + O Ser/Thr-yoğun
  heuristik. Pipeline mimikri sonrası 'glycosylation' aşaması B-epitopları işaretler (YUMUŞAK,
  elenmez); meta['glycosylation']. Rapor `_glyco_section` + Next.js GlycoCard + Streamlit paneli
  + runner katalog. tsc temiz. Doğrulandı (N-sequon saptama).
- **2026-09-27 Faz 5 ✅**: `vaxforge/cloning.py` — E. coli K-12 kodon tablosu (Kazusa), en-sık
  kodon ters-translasyon, CAI (Sharp&Li), GC%, Bio.Restriction iç bölge taraması (8 enzim),
  pET-28a(+) His6 vektör/etiket önerisi. Pipeline 'cloning' aşaması (MEV sonrası). Rapor
  `_cloning_section` + Next.js CloningCard + Streamlit + runner katalog. tsc temiz; CAI/GC doğrulandı.
- **2026-09-27 Faz 6 ✅**: `vaxforge/vaxijen.py` — VaxiJen'in yayınlanmış özellik yöntemi
  (antigen_acc: z-ölçek+ACC+lojistik) BAĞIMSIZ 2. antijenite kanıtı; canlı sunucu 403-toleranslı
  (ARAÇ_YOK). funnel her adaya VaxiJen-ACC skoru + IApred ile UYUM sayacı; funnel_methods
  ['antijenite_caprazkontrol']. Rapor/UI methods tablolarına satır eklendi (Faz1 altyapısı) +
  runner katalog. tsc temiz.
- **2026-09-27 Faz 7 ✅**: `mimicry.py` — tam k-mer eşleşme yanında BLOSUM62 **Hamming-1
  konservatif ikame** kısmi mimikri (komşu-üretme → hızlı, brute-force yok). Pipeline mimikri
  bloğu n_exact+n_near sayar; meta['mimicry']. Rapor `_mimicry_section` + Next.js MimicryCard +
  Streamlit + runner katalog. Doğrulandı (I→V BLOSUM 3.0). tsc temiz.
- **2026-09-27 Faz 8 ✅**: `vaxforge/strain_fetch.py` — NCBI Datasets v2'den taxon'a göre N suş
  proteomu (önbellekli, dürüst boş-dönüş). Pipeline: kullanıcı suş vermediyse + taxon varsa +
  ağ varsa otomatik 3 suş çeker → conservation koşar. Rapor `_conservation_section` + Next.js
  ConservationCard + Streamlit + runner katalog. CANLI test geçti (M. genitalium GCF_040556925.1
  219KB indirildi). tsc temiz.
- **2026-09-27 Faz 9 ✅**: `vaxforge/host_db.py` — DB'si olmayan konak için UniProt reviewed
  (Swiss-Prot) proteomu indir + diamond makedb (isim→taxid haritası; Host.taxid önceliklidir).
  funnel _self_homology öncesi eksik DB'leri otomatik sağlar; summary self_homoloji'ye
  'otomatik sağlandı' notu (Faz1 altyapısıyla rapor/UI'de görünür). runner katalog satırı.
  CANLI test: tavşan (9986) → 980 protein Swiss-Prot + diamond DB (539KB). tsc temiz.
- **2026-09-27 Faz 10 ✅**: `immforge_score.py` — kalibrasyon DOSYA-TABANLI (`calibrated_species()`
  = gerçekte var olan calibration.*.yaml; hardcoded liste kaldırıldı). 9 tür kalibre (human/mouse/
  bovine/pig/goat/sheep/rabbit/fish/horse); chicken → human-proxy (DÜRÜST: avian ABC kalibrasyonu
  ayrı literatür-çapalı mini-proje, sayı UYDURULMADI). Per-host calibration_species + provenance;
  rapor immüno tablosuna 'Kalibrasyon' sütunu + Next.js/Streamlit kalibrasyon türü. tsc temiz.
- **2026-09-27 Faz 11 ✅**: `tests/test_modules.py` — 11 birim testi (ağsız, deterministik):
  deeploc eşleme, structure pLDDT-parse+ölçek, rsa matematiği, bcell_conf yamaları, glyco sequon,
  cloning CAI/GC, vaxijen skor, mimicry BLOSUM Hamming-1 (konservatif-only), strain_fetch taxon,
  host_db taxid, immforge dosya-tabanlı kalibrasyon. Ağ/ağır-araç yolları smoke + manuel canlı
  testlerle kapsandı. `pytest tests/` yeşil = kalıcı regresyon güvencesi.

## 🎉 11/11 FAZ TAMAMLANDI (2026-09-27)
Tüm eksikler kapatıldı; her biri pipeline + rapor + Next.js + Streamlit'te görünür + smoke ile
korunuyor. Kurulan araçlar: PSORTb (podman), DeepLoc2, ESMFold (ESMAtlas API), NCBI Datasets
suş-çekme, UniProt konak-DB. Yeni modüller: deeploc, structure, rsa, bcell_conf, glyco, cloning,
vaxijen, strain_fetch, host_db (+ mimicry/immforge_score genişletildi).

## 2026-09-27 — Konak MHC panelleri referans setlere genişletildi + patojenite/HLA-hastalık anotasyonu
- **hosts.yaml genişletildi** (hepsi yerel NetMHCpan `-listMHC/-list`'e karşı %100 valide, atıflı):
  İnsan MHC-I 9→27 (Sidney 2008, DOI 10.1186/1471-2172-9-1), MHC-II 6→15 (Greenbaum 2011,
  DOI 10.1007/s00251-011-0513-0); Fare 5→9 / 5→8; Sığır 8→12 / 8→11 (Nielsen 2017
  DOI 10.1021/acs.jproteome.7b00675 + Fisch 2021 DOI 10.4049/jimmunol.2001409); Domuz-I 6→11.
  **Domuz MHC-II DEVRE DIŞI** (NetMHCIIpan hiç SLA-DR desteklemiyor → eski "II:1" proxy'ydi;
  tavukla tutarlı dürüst devre dışı). `predictors.py`: boş-allel guard (proxy'ye düşüp sahte
  epitop üretmeyi engeller).
- **YENİ patojenite katmanı**: `config/hla_disease.yaml` (7 allel, her satır PubMed PMID+DOI ile
  doğrulanmış: B*57:01 abacavir-advers / B*53:01 sıtma-koruyucu / B*51:01 Behçet / DRB1*15:01 MS /
  DRB1*04:01,*04:05,*01:01 RA) + `vaxforge/hla_disease.py` loader (annotate/summarize; salt-yorumlama,
  skorlamayı etkilemez). `runner.py` /api/config'e allel listeleri + hla_disease özeti eklendi.
- Patojen-tarafı virülans zaten mevcut (discovery.py VFDB) → yeniden inşa yok.
- **UI NOTU**: klinik-not kartı + risk rozetleri kullanıcı isteğiyle KALDIRILDI (bileşen `results.tsx`'te
  duruyor, veri backend'de duruyor; geri açmak tek satır). Panelde sadece I:x·II:y + tooltip kaldı.
- `tests/test_modules.py`: hla_disease birim testi (12/12 yeşil). Plan: ~/.claude/plans/delegated-painting-liskov.md

## 2026-10-10 — Makale yeniden tasarımı (Nature kalıbı) + Gram(+) lokalizasyon düzeltmesi

**Gerekçe (kullanıcı):** benchmark makalesinin yazımı acemi bulundu → hedef-dergi makaleleri
okunup sıfırdan yeniden tasarlanacak. Okunanlar: PLGDL (Zai 2025, *Nat Commun*,
10.1038/s41467-025-67778-2, tam metin), NERVE 2.0 (Conte 2024, 10.1186/s12859-024-06004-0,
tam metin), Vaxign-ML (Ong 2020), Vaxign-DL+ESM (Chen 2024 ön baskı), betimleyicisiz ESM-2
(Podda 2025, 10.1371/journal.pone.0323895). Tüm DOI'ler PubMed'den doğrulandı.

**Yapısal yeniden tasarım** (`MANUSCRIPT.md` + `MANUSCRIPT_TR.md`, docx'ler yeniden üretildi):
- Nature-ailesi sırası: Giriş → **Sonuçlar** → Tartışma → **Yöntemler** (Yöntemler artık sonda).
  Eski sıra Yöntemler'i öne koyuyordu; "acemi" duran en belirgin şey buydu.
- Özet 450+ kelimelik veri dökümünden ~215 kelimeye indirildi, tek tezli.
- Sonuçlar'dan Yöntemler'i tekrarlayan "The VaxForge pipeline" paragrafı yetenek yüzeyine
  (Tablo 1) bağlanan kısa bir mimari paragrafına dönüştü.
- Ek materyal ayrı bölüme çekildi (Tablo S1/S2, Şekil S1/S2); Şekil 5'in Şekil 7'den sonra
  geldiği sıra bozukluğu giderildi (eski fig5 → Ek Şekil 2). Şekil dosya eşlemesi başlıkta.
- **YENİ ana sonuç bölümü:** kaçan-antijen atfı (Tablo 4) artık dipnot değil merkezî argüman.
  Tartışma'da tasarım tercihi olarak savunuluyor (PLGDL'in "deliberate methodological choice"
  retoriği), HSP60 moleküler mimikri atıflarıyla (Matsuura 2009, Lucchese & Flöel 2020) +
  üç dürüstlük kaydı (eşik parametre tercihi / homoloji ≠ patoloji / moonlighting hedefte yükümlülük).
- **Yeni nesil konumlandırma eklendi:** protein-dil-modeli sınıflandırıcıları (Vaxign-DL, PLGDL,
  Podda) farklı soruyu yanıtlıyor → rakip değil, huninin skorlama aşamasının yerine geçecek
  tamamlayıcı. 2026 gönderiminin yalnız 2020 sınıflandırıcısıyla kıyaslanması açığı kapandı.
- NERVE 2.0 ile yakınsama argümanı: kendi kıyasında Vaxign-ML'i fold 1,58 ile zayıf buluyor,
  biz aynı ödünleşmeyi ters uçtan gözlüyoruz; ayrıca NERVE'in ground-truth seti sitozolik/
  bilinmeyen-lokalizasyon antijenlerini dışlamış, bizimki dışlamıyor → recall'ımız ihtiyatlı.
- Limitations açık ve tam (Gram(+) açığı kusur olarak, PSORTb LigA, epitop eşiği, ActA,
  toksisite filtresi, AUC küçük-örneklem, CPU ölçek, deneysel doğrulama yok).

**Okuma sırasında bulunan 3 METİN HATASI düzeltildi** (atıf analizi 2026-10-04 ile çelişiyordu):
1. *C. burnetii*/*S. suis* kayıpları "lokalizasyon filtresi" diye yazılmıştı → gerçekte 6 antijen
   **konak-homoloji** filtresinden düşüyor (Fisher p=3,8e-05).
2. L7/L12 "lokalizasyon eliyor" yazılmıştı → üç sert filtreyi geçip **epitop aşamasında** kayboluyor
   (Brucella aşı-suşu bölümü ve Ek Şekil 1 altyazısı da düzeltildi).
3. ActA "zor hedef" diye geçiştirilmişti → **Gram(+) yapılandırma açığı**.

**Kod düzeltmesi — Gram-duyarlı lokalizasyon:**
- `config/thresholds.yaml`: `localization.allowed_gram_positive = [extracellular, cell_wall, membrane]`.
- `vaxforge/funnel.py`: `profile == "bacteria" and gram == "positive"` → bu listeyi seçer.
- `vaxforge/evaluate.py`: `candidate_rows(..., gram=...)` → rapor da gram-duyarlı
  (`vaxforge/report.py` çağrısı güncellendi). `vaxforge/psortb.py` bayat yorumları düzeltildi.
- Özgüllük güvencesi: `transmembrane.max_helices` gömülü çok-sarmallıları ayıklamaya devam eder;
  `cytoplasm` hiçbir Gram tipinde kabul edilmez.
- Test: `tests/test_modules.py::test_localization_allowed_is_gram_aware`. **15/15 yeşil**
  (13 modül + 2 smoke, uçtan uca pipeline dahil, 241 s).

**KALAN:** Gram(+) dalının yeniden koşulması — yalnız 5 tür (saureus/listeria/ssuis/mbovis/
banthracis); Gram(−) 5 tür etkilenmiyor. Beklenti: CapB + ActA geri gelir → 37/48 → 39/48,
Gram(+) türlerde fold bir miktar düşer. Makale bu nedenle recall'ı açıkça "taban" olarak
çerçeveliyor, yani metin yeniden koşumdan önce de tutarlı.

## 2026-10-10 (II) — Filtre kalibrasyonu: toksisite kapısı yanlış uzunluk sınıfında kullanılıyordu

**Tetikleyici.** Gram(+) düzeltmesinden sonra ActA huniyi geçmeye başladı ama yine
"aday peptit üretmedi" durumunda kaldı; listeria'da LLO da kayboldu. Epitop
aşamasındaki süzmeyi ölçtük.

**Bulgu.** Toksisite kapısı ToxinPred2 Model 2 (hibrit, eşik 0.6) ile çalışıyordu.
Bağımsız kontrol setlerinde (n=1000): insan Swiss-Prot self 9-mer'lerinin **%75.2**'si,
deneysel **doğrulanmış** IEDB 9-mer epitoplarının **%63.1**'i "toksik" işaretlendi.
Influenza M1(58-66) **GILGFVFTL** (ml_score 0.663) dahil. Aynı araç TAM PROTEİNDE
%94.3 duyarlılık / %4.7 yanlış-pozitif veriyor. Sorun araç değil **uzunluk sınıfı**:
ToxinPred2 "toxicity of *proteins*" için eğitilmiştir, AAC özellikleri 9-mer'de
bilgi taşımaz. Epitop aşaması pratikte rastgele bir 4'te-3 kırımıydı; az epitop
üreten proteinler (ESAT-6, L7/L12, ActA, LLO) bu yüzden tümüyle kayboluyordu.

**Düzeltme.**
- `vaxforge/toxinpred3.py` (YENİ): ToxinPred3 sarmalayıcı (Rathore 2024,
  doi:10.1016/j.compbiomed.2024.108926) — peptit düzeyi için eğitilmiş.
  `scripts/convert_toxinpred3_model.py` ile sklearn köprüsü (IFNepitope2 ile aynı
  yamalı unpickler); upstream'in tek-dizi girdide çöken AAC/DPC birleştirmesi
  sarmalayıcıda kukla dizi ile aşıldı.
- `config/thresholds.yaml` → `toxicity.hard_filter: false`; yeni paramlar
  `hard_discard_score` (yüksek-güven eleme) ve `toxicity_penalty` (×0.9).
- `vaxforge/survival.py` yeniden yazıldı: epitopta ToxinPred3 **işaretler**,
  elemez; ToxinPred2'ye DÜŞMEZ. Gerçek araç yoksa kapı ATLANIR — eskiden
  demonstrasyon sezgiselleri (kompozisyon vekilleri) sessizce sert filtre gibi
  davranıp özellikle hidrofobik MHC-I bağlayıcılarını eliyordu.
- `vaxforge/scoring.py`: `toxicity_penalty_factor` yumuşak cezası (mimikriyle aynı desen).
- `vaxforge/mev.py`: konstruktta sert toksisite elemesi ToxinPred2 hibritle KALIYOR
  (orada kalibre). Alerjenite konstruktta artık FAO/WHO **birincil** ölçütüyle
  (80 aa pencerede ≥%35 kimlik, diamond → `allergen.predict_protein`): 6-mer kuralı
  tam proteinde rastgele non-alerjen proteinlerin 8/8'ini işaretliyordu, yeni ölçüt
  0/8 işaretliyor ve gerçek alerjenleri %100 kimlikle yakalıyor.
- `vaxforge/pipeline.py`: sağ kalım özeti tümüyle `meta["survival"]`e yazılıyor.
  (`n_epitope` ÜRETİLEN peptit sayısıdır; eleme yalnız bağlanma eşiğini geçenlere
  uygulanır — ikisi karıştırılınca alerjenite kapısı olduğundan agresif görünüyordu.)

**Kalıcı düzenek.** `benchmark/calibrate_filters.py` + `write_calibration_md.py`:
her kapıyı insan self k-mer'leri, doğrulanmış IEDB epitopları ve uzunluk-eşli
protein setlerinde tarar; `results/filter_calibration.json` + `FILTER_CALIBRATION.md`.
Testler: `test_toxicity_gate_is_soft_and_length_appropriate`,
`test_allergen_criterion_matches_sequence_length`. 15/15 modül testi yeşil.

**Kayıt düzeltmesi.** `MISS_ATTRIBUTION.md`'deki "mühendislik hatası yok" sonucu
geri alındı; 4. kategori (2 kayıp) + ActA + LLO bu kalibrasyon hatasına atfedildi.

**Durdurulan iş.** Gram(+) yeniden koşumu (ssuis, listeria bitmişti; saureus
yarıda) kullanıcı talimatıyla DURDURULDU — bozuk toksisite kapısıyla üretilen
sayıların değeri yoktu. Yeniden koşum kalibre edilmiş kapıyla, tek kod sürümünde
10 türün tamamı için yapılacak.

**Kayıt notu.** Yarım kalan Gram(+) koşumunun yazdığı sonuç dosyaları
(`results/{listeria,ssuis}*`) geri alındı: bozuk toksisite kapısıyla üretildiler ve
kalan 8 tür eski kodla kaldığı için kıyas tablosunu iç tutarsız yapıyorlardı.
Depodaki kıyas verisi şimdilik tümüyle düzeltme-öncesi; yeniden koşum kalibre edilmiş
kapıyla tek kod sürümünde yapılacak. (Ölçülen ham değerler kayıt için: ssuis
aday 288→422, recall 1/5 değişmedi, fold 1.34→0.91; listeria aday 478→538,
recall 4/5→3/5 — ActA huniyi geçti ama toksisite kapısında kaldı, LLO da aynı kapıda.)

**Örnek veri doğrulaması (2026-10-10).** `benchmark/demo_toxicity_gate.py` örnek
proteom koşusunun nihai aday listesi (26 peptit) üzerinde iki kapıyı yan yana
değerlendirir. Eski kapı (ToxinPred2 Model 2 hibrit, t=0.6, sert) bu adayların
16'sını (%61.5) elerdi; yeni kapı (ToxinPred3 ExtraTrees, işaretleme ≥0.6, eleme
≥0.85) hiçbirini işaretlemedi. Aynı koşunun sağ kalım özeti: 51 girdi peptit →
51 peptit (0 alerjen elenen, 0 toksik işaretlenen). Veri:
`benchmark/results/demo_toxicity_gate.json`. Betik `--run-json` ile var olan bir
koşu üzerinden tekrar koşulabilir (pipeline'ı yeniden çalıştırmaz).

## 2026-10-10 (III) — Otoimmünite kapısının granülaritesi

**Bulgu.** Konak-homoloji kapısı (`funnel.py`, `min_identity=35`, sert) protein
BÜTÜNÜNDE kimlik ölçüp tüm antijeni eliyordu. S. suis'in kaçan 4 antijeni
(enolaz, GAPDH, 6-fosfoglukonat dehidrogenaz, DnaJ) tam bu kapıda, epitop
tahmini hiç koşmadan düşüyordu. Ölçüm (`benchmark/measure_homology_granularity.py`):

| Antijen | protein-bütünü kimlik | 9-mer | self-tam eşleşen | oran | güvenli k-mer |
|---|---|---|---|---|---|
| gnd | %52.3 | 467 | 14 | %3.00 | 453 |
| Enolaz | %49.8 | 427 | 8 | %1.87 | 419 |
| GAPDH | %43.9 | 328 | 7 | %2.13 | 321 |
| DnaJ | %37.0 | 370 | 5 | %1.35 | 365 |
| Suilizin (yakalanan) | %0 | 489 | 0 | %0 | 489 |

Yani protein-bütünü kimlik, epitop-düzeyi otoimmünite riskinin vekili değil:
enolazda 8 riskli 9-mer için 419 self-olmayan epitop da atılıyordu. Toksisite
kapısıyla aynı hata ailesi — orada yanlış uzunluk sınıfı, burada yanlış
granülarite.

**Düzeltme (toksisitedekinin aynadaki görüntüsü).**
- `human_homology`: `hard_filter: false`. ≥%35 → işaretlenir + skor ×0.75;
  eleme yalnız `hard_discard_identity=90` üstünde (konak proteininin kendisi /
  kontaminasyon güvenliği). Ceza `funnel.py`'de annotasyona yazılır,
  `pipeline.py` bunu peptitlere taşır, `scoring.py` uygular.
- `molecular_mimicry`: `hard_filter: true`. Konak self ile TAM eşleşen epitop
  elenir; BLOSUM-yakın eşleşme elenmez, skor ×0.9. Kalibrasyon dayanağı:
  n=1000'de insan self k-mer'lerinde %100 yakalama, DENEYSEL DOĞRULANMIŞ IEDB
  epitoplarında %1.3 (9-mer) / %1.7 (15-mer) yanlış-pozitif.
- Rapor katmanları hizalandı (`evaluate.py` iki satır, `runner.py` iki adım,
  `mimicry.py` modül açıklaması).
- Yeni test: `test_autoimmunity_gate_is_at_epitope_granularity` +
  `test_homology_granularity_measurement_recorded`. 17/17 yeşil.

**Beklenen etki ve riski.** Bu 4 antijen artık huniden geçebilir; ama geçmesi
garanti değil — lokalizasyon/TM/antijenite kapıları hâlâ önlerinde ve skorları
×0.75 ceza alıyor. Karşı taraf dürüstçe: tam eşleşme otoimmünitenin tamamı
değil, TCR dejenerasyonu nedeniyle birebir olmayan peptit de çapraz-reaktif
olabilir (anti-enolaz / anti-GAPDH otoantikorları literatürde gerçek). Bu yüzden
yakın eşleşme cezası korunuyor ve protein-düzeyi işaret raporda görünür kalıyor.
