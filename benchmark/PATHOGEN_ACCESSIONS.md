# Benchmark patojenleri — veri kaynağı & aksesyon numaraları (yayın Methods için)

Tüm aksesyonlar 2026-08-19'da NCBI/UniProt'tan **doğrulandı** (uydurma yok).
Yayında Methods → "Data sources" bölümüne birebir girecek.

## Bakteri kolu — MenB (çok-antijen recall + 3'lü head-to-head)
- **Organizma:** *Neisseria meningitidis* serogroup B, strain ATCC BAA-335 / **MC58**
- **Taxonomy ID:** 122586 (MC58 strain) · 487 (species, IEDB için kullanılan)
- **Proteom kaynağı:** UniProt/Swiss-Prot, proteome **UP000000425** (2003 protein) → `data/validation/menb_mc58_proteome.faa`
- **NCBI RefSeq genom:** **NC_003112.2** (kromozom, "Neisseria meningitidis MC58, complete sequence")
- **NCBI assembly:** **GCF_000008805.1** (ASM880v1)
- **Ground-truth koruyucu antijenler (UniProt, 4CMenB/Bexsero):**
  - **Q9JXV4** — fHbp (Factor H binding protein, NMB1870)
  - **Q9JXK7** — NadA (Neisseria adhesin A, NMB1994)
  - **Q7DD37** — NHBA (Neisserial Heparin Binding Antigen, NMB2132) ← VFDB-negatif kilit deney
  - **P0DH58** — PorA (Major outer membrane protein P.IA, OMV bileşeni)
- **Test seti:** 4 koruyucu + 60 arka plan = 64 protein, SEED=42 → `data/validation/menb_testset.faa`

## Virüs kolu — Rabies (çok-konak MHC gösterimi)
- **Organizma:** *Lyssavirus rabies* (Rabies virus, RABV)
- **Taxonomy ID:** **11292**
- **NCBI RefSeq genom:** **NC_001542.1** (tam genom, 5 protein) → `data/validation/rabies_proteome.faa`
- **Protein aksesyonları (RefSeq):**
  - **NP_056793.1** — nucleoprotein N (GeneID 1489853)
  - **NP_056794.1** — phosphoprotein P/M1 (GeneID 1489854)
  - **NP_056795.1** — matrix M2 protein (GeneID 1489855)
  - **NP_056796.1** — transmembrane **glycoprotein G** (GeneID 1489856) ← **GROUND-TRUTH antijen**
  - **NP_056797.1** — L protein / polymerase (GeneID 1489857)
- **Konaklar (çok-konak MHC):** insan, sığır (bovine), fare (mouse), domuz (pig) — dördü de
  gerçek konak (biyolojik gerekçe + DOI'ler yazım aşamasında eklenecek).

## Bakteriyel kol — 3 bakteri × TAM PROTEOM (2026-08, çok-antijen recall head-to-head)
4 memeli konağı (insan/sığır/fare/domuz) enfekte eden, genom-boyu artan 3 zoonotik bakteri.
Her araç DEFAULT + özdeş çıplak-aksesyon TAM PROTEOM girdi. Config: `benchmark/bacteria/<bug>.json`.

- **Staphylococcus aureus** NCTC 8325 · UniProt proteome **UP000008816** (2889 protein) · taxid 93061 (suş)/1280 (tür)
  · RefSeq **NC_007795.1** · asm **GCF_000013425.1** · Gram+ · `data/validation/saureus_prepared.faa`
  - Koruyucu antijenler (UniProt, GN= ile doğrulandı): **Q2FZF0** IsdB · **Q2G015** ClfA · **Q2G1X0** Hla(α-hemolizin)
    · **Q2FZE9** IsdA · **P02976** SpA
- **Listeria monocytogenes** EGD-e · UniProt **UP000000817** (2844) · taxid 169963/1639 · RefSeq **NC_003210.1**
  · asm **GCF_000196035.1** · Gram+ · `data/validation/listeria_prepared.faa`
  - Antijenler: **P13128** LLO(listeriolizin O) · **P21171** p60(iap) · **P0DJM0** InlA · **P0DQD2** InlB · **P33379** ActA
- **Salmonella enterica** ser. Typhimurium LT2 · UniProt **UP000001014** (4533) · taxid 99287/90371 · RefSeq **NC_003197.2**
  · asm **GCF_000006945.2** · Gram− · `data/validation/salmonella_prepared.faa`
  - Antijenler: **P37592** OmpD · **P0A263** OmpC · **P06179** FliC(flagellin) · **Q7BVH7** SseB · **P02936** OmpA
  - NOT: Vaxign-ML için selenosistein(U) içeren 3 protein sanitize edildi (U→C) → `salmonella_prepared_vaxign.faa`
    (VaxForge kendi içinde sanitize eder; 5 koruyucu antijenin hiçbiri U içermez → sonuç etkilenmez).

## (İleride, opsiyonel referans) SARS-CoV-2
- UniProt/Swiss-Prot, `data/validation/sarscov2_proteome.faa`; Spike = **P0DTC2** (SPIKE_SARS2),
  taxid 2697049. Şu an benchmark'ta kullanılmıyor.

---
NOT: Patojen + konak **seçilim/biyolojik-kullanım gerekçeleri** DOI-doğrulanmış makalelerle
desteklenecek (kullanıcı kuralı 2026-08-19; yazım aşamasında). Bkz [[vaxforge-competitor-benchmark]],
[[vaxforge-literature-rule]], [[convir-citation-policy]].

## Genişletilmiş bakteri paneli — +6 zoonotik bakteri (2026-10-01, 10'a tamamlama)
Mevcut 4'e (MenB, S.aureus, Listeria, Salmonella) ek olarak 4 memeli konağı (insan/sığır/fare/domuz)
enfekte eden 6 zoonotik bakteri. Tümü UniProt referans proteomu + NCBI assembly ile **doğrulandı**
(2026-10-01). Koruyucu antijen aksesyonlarının hepsi ilgili proteomda bulunuyor (prep_proteome 100%).
Config: `benchmark/bacteria/<bug>.json`. **Genom küçükten büyüğe koşu sırası.**

1. **Coxiella burnetii** RSA 493/Nine Mile · UP000002671 (1812) · taxid 227377/777 · GCF_000007765.2 · 2.03 Mb · Gram−
   - H7C7D7 Com1 · P51752 Mip · Q83DT1 OmpH · P19421 GroEL
2. **Streptococcus suis** ser.2 BM407 · UP000009077 (1922) · taxid 568814/1307 · GCF_000026745.1 · 2.17 Mb · Gram+
   - A0A0H3MU02 Suilysin(SLY) · A0A0H3N5B2 Enolaz · A0A0H3MSV3 GAPDH(plr) · A0A0H3N193 6PGD(gnd) · A0A0H3MXM0 DnaJ
   - NOT: referans P1/7 (UP000002055) UniProtKB'de boşaltılmış → aynı serotip-2 klinik izolat BM407.
3. **Brucella suis** biovar 1 · 1330 · UP000007104 (3249) · taxid 204722/29461 · GCF_000223195.1 · 3.32 Mb · Gram−
   - Q45689 Omp25 · P0A3U5 Omp31 · P0A3U9 BP26 · P66827 SodC · P0A469 L7/L12(RplL)
   - Aşı suşları (KOŞULDU; virülan-vs-aşı antijen recall — Tablo S2 + Şekil S1):
     - **B.abortus RB51** (rough, canlı-atenüe) · GCF_000366025.1 · 2971 protein · taxid 235 · 4 antijen (Omp31 YOK, B.abortus delesyonu; Cloeckaert 2002) · recall 3/4
     - **B.abortus S19** (smooth, canlı-atenüe) · GCF_000018725.1 · 2964 protein · taxid 235 · 4 antijen (Omp31 YOK) · recall 3/4
     - **B.suis S2** (oral, Çin, canlı-atenüe) · GCF_000209635.2 · 2979 protein · taxid 29461 · 5 antijen (Omp31 VAR) · recall 4/5
     - Her üçünde de tek kayıp L7/L12 (sitoplazmik → lokalizasyon filtresi, tasarım gereği).
4. **Mycobacterium bovis** AF2122/97 · UP000001419 (4088) · taxid 233413/1765 · GCF_000195835.3 · 4.35 Mb · Gram+(asit-dirençli)
   - P0A565 ESAT-6(EsxA) · P0A567 CFP-10(EsxB) · P0C2T2 Ag85B(FbpB) · P0A5Q5 MPB64 · P0A5B8 HspX(Acr)
5. **Leptospira interrogans** serovar Lai 56601 · UP000001408 (3676) · taxid 189518/173 · GCF_000092565.1 · 4.7 Mb · Gram−(spiroket)
   - Q8EYU4 LigA · Q8EZS3 LigB · O34094 LipL32 · Q8F8E1 LipL41 · G1UB30 OmpL1
6. **Bacillus anthracis** Ames Ancestor · UP000000594 (5493, kromozom+pXO1+pXO2) · taxid 261594/1392 · GCF_000008445.1 · 5.5 Mb · Gram+
   - P13423 PA(pagA,pXO1) · P15917 LF(lef,pXO1) · P40136 EF(cya,pXO1) · P19580 CapB(pXO2) · P94217 EA1 S-tabaka
   - NOT: PA/LF/EF yalnız pXO1'de → plazmit-içeren proteom (UP000000594) seçildi; UP000000427'de toksin yok.

Koşu planı: her bakteri için ÖNCE rakipler (Vaxign-ML + NERVE2, `run_competitors.sh`), SONRA VaxForge
(`run_bacterium.py`). Genom sırası: coxiella → ssuis → brucella → mbovis → leptospira → banthracis.
