# Patentli / patent-başvurulu peptit-protein aşıları — benchmark patojenleri

Tarama: 2026-09-29 (web/Google Patents/USPTO). Antijen = protein/peptit hedefi; tip = subunit vs
canlı-attenüe/rekombinant-virüs ayrımı. **DÜRÜSTLÜK:** patent no'ları arama sonuçlarından; ground-truth
antijenle birebir eşleşen kritik olanlar yazıma girmeden Google Patents'ten tek tek doğrulanacak
(uydurma yasağı — VaxForge literatür kuralı).

## BAKTERİ

| Patojen | Antijen (protein/peptit) | Patent / başvuru no | Tip | Kaynak |
|---|---|---|---|---|
| N. meningitidis B | fHbp (modifiye) | US 11,708,394 · US 12,312,384 "Modified meningococcal fHbp polypeptides" | Rekombinant protein | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11708394 |
| N. meningitidis B | NHBA+fHbp+NadA (Bexsero/4CMenB, ruhsatlı) | ürün onaylı (DailyMed) | Rekombinant fusion protein | https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=f70cf2fc-6e6d-4a74-9f7a-db8fec072fd7 |
| S. aureus | ClfA + IsdB (± diğer) | CA 2758490 A1 "Compositions for immunising against S. aureus" | Rekombinant protein komb. | https://patents.google.com/patent/CA2758490A1/en |
| S. aureus | Hla + ClfA + SpA + kapsül | US 2022/0395566 A1 · US 2021/0023192 A1 · EP 3150222 B1 | İmmünojenik kompozisyon | https://patents.justia.com/patent/20220395566 · https://patents.google.com/patent/EP3150222B1/en |
| L. monocytogenes | LLO (non-hemolitik trunke) | US 6,855,320 "Fusion of truncated LLO to antigens" | Protein füzyon platformu | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6855320 |
| Salmonella | OmpD, FliC (+ PltB/CdtB) | US 2019/0192646 A1 "Salmonella vaccines" | Subunit | https://patents.justia.com/patent/20190192646 |
| Salmonella | Tip-III sekresyon proteinleri | US 11,464,843 | Koruyucu protein aşısı | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11464843 |

## VİRÜS

| Patojen | Antijen (protein/peptit) | Patent / başvuru no | Tip | Kaynak |
|---|---|---|---|---|
| Rabies | G glikoprotein peptidi | CN 102174110 B · US 5,066,593 · US 4,707,356 (sentetik peptit) | Peptit subunit | https://patents.google.com/patent/CN102174110B/en · https://patents.google.com/patent/US5066593A/en |
| Rabies | G glikoprotein (gen/rekombinant) | EP 0324157 A1 · US 2016/0166711 A1 | Rekombinant protein | https://patents.google.com/patent/EP0324157A1/en |
| VSV | G glikoprotein (mutant) | US 12,264,180 B2 · US 12,091,434 B2 · US 2020/0165302 A1 · WO 2019/057974 A1 | Modifiye glikoprotein | https://patents.google.com/patent/US12264180B2/en |
| EMCV | kapsid (VP1 dahil) | US 5,213,795 A "Encephalomyocarditis virus vaccine" | Aşı (VP1 DNA literatürde) | https://patents.google.com/patent/US5213795A/en |
| VEEV | rekombinant (E2 içeren) | US 6,936,257 / RE41745 · WO 1999/063098 A1 | Rekombinant aşı | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6936257 |
| JEV | E (envelope) proteini | US 5,494,671 (trunke flavivirüs E) · US 5,021,347 (vaccinia-E) | Rekombinant subunit | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/5021347 |
| RVFV | Gn + Gc glikoprotein | US 9,791,445 · WO 2014/117185 A1 · EP 2948552 A1 · CA 2899659 A1 | Subunit | https://patents.google.com/patent/WO2014117185A1 |
| Reovirus | sigma-1 / sigma-C | JP 2024512575 A (sigma-1 platform) · WO 2009/093251 A2 (sigma-C peptit) | Peptit/protein | https://patents.google.com/patent/JP2024512575A/en |
| Getah | E2 subunit'e özel patent BULUNAMADI | EP 0562136 A1 (canlı aşı, subunit değil) | Canlı-attenüe | https://pubchem.ncbi.nlm.nih.gov/patent/EP-0562136-A1 |
| Vaccinia/Cowpox | L1R + A33R (+ B5R/A27L) | NIH E-168-2003 (başvuru 60/451,337, TERK EDİLMİŞ) | Rekombinant subunit | https://www.techtransfer.nih.gov/patent/e-168-2003-0-us-01 |

## Dürüst notlar
- **Getah**: E2 subunit'e özel patent yok; sadece canlı-attenüe (EP 0562136). Ground-truth E2 seçilirse literatürle desteklenecek.
- **Vaccinia** L1R/A33R subunit patenti terk edilmiş (aktif değil).
- VEEV (US6936257), EMCV (US5213795) canlı/rekombinant virüs patenti — saf peptit/protein subunit değil.
