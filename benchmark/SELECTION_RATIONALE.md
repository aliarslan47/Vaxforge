# Benchmark patojen & konak seçimi — biyolojik/seçilim gerekçeleri (DOI-atıflı)

Bu doküman, VaxForge benchmark'ında neden **MenB** (bakteri kolu) + **rabies** (virüs kolu) ve
neden **insan / sığır / fare / domuz** konaklarının seçildiğini bilimsel gerekçeye bağlar. Yayın
Methods → "Rationale for benchmark selection" bölümüne birebir girecek. Aksesyonlar:
`PATHOGEN_ACCESSIONS.md`. Sonuçlar: `BENCHMARK_SUMMARY.md`.

**KURAL:** Her iddia gerçek, web-doğrulanmış DOI'ye bağlı — uydurma YOK. Doğrulama tarihi 2026-08-20
(başlık + yazar + yıl + dergi eşleşmesi WebSearch/WebFetch ile teyit edildi). Doğrulanamayan hiçbir
atıf konmadı (bu turda tümü doğrulandı).

---

## 1. Neden MenB (bakteri kolu) — çok-antijen recall zemini

Serogrup B *Neisseria meningitidis* (MenB), **ters aşı bilişiminin doğduğu** patojendir: Pizza ve
ark. (2000) tüm-genom dizilemesinden aday antijenleri hesaplayarak bu paradigmayı kurdu
[Pizza 2000]. Dolayısıyla MenB, bir RV/RV-benzeri aracın "altın standart" test zeminidir —
alan bu patojen üzerinde tanımlandı.

Ruhsatlı **4CMenB (Bexsero)** aşısının koruyucu bileşenleri deneysel olarak doğrulanmıştır:
rekombinant **fHbp, NadA, NHBA** + PorA taşıyan OMV [Giuliani 2006]. Bu, tek bir "doğru cevap"
değil, **çok-antijenli, deneysel-doğrulanmış bir ground-truth seti** verir → bir aracın
**recall@antijen** (bilinen koruyucuların kaçını geri buluyor) ölçütü için ideal. VaxForge tezi
tam da burada sınanır: **NHBA VFDB-negatiftir** ve sert virülans-filtreli araçlarca elenme
riski taşır; yumuşak-filtre onu korumalıdır (bkz `BENCHMARK_SUMMARY.md` §1, kilit deney).

## 2. Neden rabies (virüs kolu) — tek-antijen ama çok-konak sahnesi

Rabies virüsünün (RABV) tek yüzey proteini **glikoprotein G**'dir ve **virüs-nötralize edici
antikorları indükleyen + koruma sağlayan TEK yapısal proteindir** — Cox, Dietzschold & Schneider
(1977) saflaştırılmış G'nin bu özelliği taşıdığını, diğer yapısal proteinlerin taşımadığını
gösterdi [Cox 1977]. Bu, **kesin ve tek** bir ground-truth (G = koruyucu antijen) verir →
virüs kolunda araç, "doğru antijeni en üste koyuyor mu" diye net sınanır (VaxForge G #2/5).

Rabies'in benchmark'a asıl katkısı **çok-konak** boyutudur: RABV neredeyse **tüm memelileri**
enfekte eder ve konak yelpazesi tür-ötesidir [Fooks 2017]. Bu yüzden **insan, sığır, fare, domuz**
dördü de RABV'nin **gerçek doğal memeli konaklarıdır** — her biri için MHC gösterimi biyolojik
olarak dürüsttür (uydurma konak yok). Influenza (sığırda marjinal), FMDV (insanı enfekte etmez)
veya SARS-CoV-2 (sığır/domuz doğal değil) bu dört-konak matrisini dürüstçe dolduramazdı; rabies
doldurur. Böylece VaxForge'un **konak-spesifik MHC epitop repertuarı** yeteneği (NERVE/Vaxign'de
yok) adil bir zeminde gösterilir: G-epitopları host×host Jaccard <1.0 → konak seçimi öngörülen
epitopları esaslı değiştirir (`BENCHMARK_SUMMARY.md` §2).

### Konak seçimi notu (dürüst kapsam)
Dört konak da RABV'nin doğal memeli konağı olduğu için seçildi [Fooks 2017]; sığır (veteriner/
ekonomik yük), domuz (tarımsal + insana yakın), fare (standart rabies challenge/lab modeli) ve
insan (halk sağlığı hedefi) pratikte de anlamlıdır. Bu dört-konak iddiası **tek bir otoriter
derlemeye** (Fooks 2017) dayandırılmıştır; konak-başına ayrı epidemiyolojik DOI eklenmesi yazım
sırasında opsiyoneldir — mevcut hâliyle hiçbir uydurma atıf yoktur.

## 3. Neden bu iki rakip — NERVE 2.0 & Vaxign-ML

Karşılaştırma araçları da atıfla gerekçelidir: **NERVE 2.0** güncel, AI-destekli, lokal
kurulabilir bir RV ortamıdır [Conte 2024]; **Vaxign-ML** bağımsız, denetimli-ML tabanlı bakteriyel
koruyucu-antijen tahmincisidir [Ong 2020]. İkisi de default parametre + özdeş girdiyle koşuldu
(adil koşul); bağımsız-ML olan Vaxign-ML'in de 4/4 recall vermesi, VaxForge'un soft-filter tezinin
tek araca özgü bir yapaylık olmadığını gösterir.

---

## Atıf listesi (Methods-hazır, hepsi DOI-doğrulanmış 2026-08-20)

1. **Pizza M, Scarlato V, Masignani V, et al.** Identification of vaccine candidates against
   serogroup B meningococcus by whole-genome sequencing. *Science*. 2000;287(5459):1816–1820.
   **doi:10.1126/science.287.5459.1816** · PMID 10710308.
2. **Giuliani MM, Adu-Bobie J, Comanducci M, et al.** A universal vaccine for serogroup B
   meningococcus. *PNAS*. 2006;103(29):10834–10839. **doi:10.1073/pnas.0603940103** · PMID 16825336.
3. **Cox JH, Dietzschold B, Schneider LG.** Rabies virus glycoprotein. II. Biological and
   serological characterization. *Infect Immun*. 1977;16(3):754–759.
   **doi:10.1128/iai.16.3.754-759.1977** · PMID 408269.
4. **Fooks AR, Cliquet F, Finke S, et al.** Rabies. *Nat Rev Dis Primers*. 2017;3:17091.
   **doi:10.1038/nrdp.2017.91** · PMID 29188797.
5. **Conte A, Gulmini N, Costa F, et al.** NERVE 2.0: boosting the new enhanced reverse vaccinology
   environment via artificial intelligence and a user-friendly web interface. *BMC Bioinformatics*.
   2024;25:377. **doi:10.1186/s12859-024-06004-0** · PMID 39695945 · PMC11654298.
6. **Ong E, Wang H, Wong MU, Seetharaman M, Valdez N, He Y.** Vaxign-ML: supervised machine learning
   reverse vaccinology model for improved prediction of bacterial protective antigens.
   *Bioinformatics*. 2020;36(10):3185–3191. **doi:10.1093/bioinformatics/btaa119** · PMID 32096826.

> İlişkili kayıtlar: `PATHOGEN_ACCESSIONS.md` (aksesyonlar), `BENCHMARK_SUMMARY.md` (§5 buraya
> çapraz-referans verir). Çekirdek ilke: her şey literatüre atıflı, uydurma yok.

---

## Bakteriyel kol genişletmesi — 3 bakteri seçim gerekçesi (2026-08)

MenB'nin tek-örnekliğinden çıkmak için 3 ek bakteri seçildi. Ölçütler:
1. **4 memeli konağını da (insan/sığır/fare/domuz) gerçekten enfekte etme** — virüs kolundaki çok-konak
   mantığının bakteriyel karşılığı; hepsi zoonotik/geniş-konaklı. (Konak-aralığı DOI'leri yazımda eklenecek.)
2. **Genom-boyu yayılımı** (2.8 → 4.9 Mb) — pipeline'ın farklı proteom ölçeklerinde tutarlılığını sınar
   (ve tam-proteom ölçek-bug'larını açığa çıkardı → BENCHMARK_SUMMARY §metodoloji).
3. **Gram dengesi** (2 Gram+ / 1 Gram−) — PSORTb Gram-özgü lokalizasyonunu iki sınıfta da test eder.
4. **Bilinen koruyucu/aday antijenler** — ölçülebilir recall@antijen için literatürde yerleşik hedefler.

Seçilenler + ground-truth antijen gerekçesi (DOI'ler yazımda doğrulanacak, [[vaxforge-literature-rule]]):
- **S. aureus NCTC 8325** (Gram+, 2.8 Mb): IsdB, ClfA, Hla, IsdA, SpA — klinik/preklinik aşı adayları
  (IsdB = Merck V710; ClfA/Hla/SpA-mutant çok-bileşenli adaylar).
- **L. monocytogenes EGD-e** (Gram+, 2.9 Mb): LLO (immünodominant CD8 hedefi), p60, InlA, InlB, ActA.
- **S. Typhimurium LT2** (Gram−, 4.9 Mb): OmpD, OmpC, OmpA, FliC(flagellin), SseB.
