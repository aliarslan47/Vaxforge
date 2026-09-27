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
