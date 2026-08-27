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
