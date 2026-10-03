# Ground-truth koruyucu antijenler — literatür kanıtı & DOI atıfları

Bu doküman, VaxForge benchmark'ında **recall@antijen** ölçütü için "doğru cevap" kabul edilen
her ground-truth antijenin **neden koruyucu/aday antijen** sayıldığını literatüre bağlar.
Yayın Methods → "Definition of ground-truth protective antigens" bölümüne birebir girer.

**KURAL ([[vaxforge-literature-rule]]):** Her atıf gerçek ve doğrulanmış bir DOI'ye bağlıdır —
**uydurma YOK**. Doğrulama tarihi **2026-09-20**; her DOI ya yayıncı sayfasından ya da Crossref
API'sinden (başlık + yazar + yıl + dergi + cilt/sayfa eşleşmesi) teyit edildi. Doğrulanamayan
hiçbir atıf konmadı.

**Dürüstlük notu:** Kanıt gücü antijene göre değişir ve aşağıda AÇIKÇA belirtilir — bazıları
doğrudan aktif-immünizasyon koruma verisi (ör. IsdB, LLO, OmpD), bazıları pasif-immünizasyon
koruması (ör. ActA) ve bazıları invazyon-virülans faktörü + aday-antijen statüsündedir
(ör. InlA/InlB). Fazla-iddia yok. Aksesyonlar: `PATHOGEN_ACCESSIONS.md`. Sonuçlar:
`BENCHMARK_SUMMARY.md`. Patojen/konak seçim gerekçesi: `SELECTION_RATIONALE.md`.

---

## 1. *Neisseria meningitidis* serogrup B (MenB)

| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **fHbp** | Ruhsatlı 4CMenB (Bexsero) rekombinant koruyucu bileşeni; bakterisidal antikor indükler | Giuliani 2006 — **10.1073/pnas.0603940103** |
| **NadA** | Ruhsatlı 4CMenB rekombinant koruyucu bileşeni; adezin, bakterisidal | Giuliani 2006 — **10.1073/pnas.0603940103** |
| **NHBA** | Ruhsatlı 4CMenB rekombinant koruyucu bileşeni (**VFDB-negatif** → kilit deney) | Giuliani 2006 — **10.1073/pnas.0603940103** |
| **PorA** | 4CMenB'deki OMV bileşeninin baskın serosubtip-spesifik bakterisidal hedefi; RV'nin doğuş patojeni | Pizza 2000 — **10.1126/science.287.5459.1816** |

> fHbp/NadA/NHBA üçü de Giuliani ve ark.'nın "universal vaccine" çalışmasında deneysel olarak
> koruyucu (bakterisidal) gösterildi. NHBA'nın VFDB-negatifliği VaxForge soft-filter tezinin
> kilit deneyidir (bkz `BENCHMARK_SUMMARY.md` §1).

## 2. Rabies virüsü (RABV)

| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **Glikoprotein G** | Virüs-nötralize edici antikor indükleyen + koruma sağlayan **TEK** yapısal protein | Cox 1977 — **10.1128/iai.16.3.754-759.1977** |

## 3. *Staphylococcus aureus* (NCTC 8325)

| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **IsdB** | Aktif immünizasyon murin sepsis modelinde sağkalımı artırır; Merck V710 adayı | Kuklin 2006 — **10.1128/iai.74.4.2215-2223.2006** |
| **IsdA** | Anti-IsdA antikorları apse oluşumu + letal challenge'a karşı korur (hem-demir kısıtı) | Kim 2010 — **10.1016/j.vaccine.2010.02.097** |
| **ClfA** | Rekombinant ClfA immünizasyonu deneysel artrit/sepsise karşı korur | Josefsson 2001 — **10.1086/324430** |
| **Hla (α-hemolizin)** | Yapısal-tasarlanmış toksoid aşı bakteriyemi + pnömoniye karşı korur | Adhikari 2012 — **10.1371/journal.pone.0038567** |
| **SpA (protein A)** | Non-toksijenik SpA_KKAA aşısı MRSA challenge'a karşı korur | Kim 2010 — **10.1084/jem.20092514** |

## 4. *Listeria monocytogenes* (EGD-e)

| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **LLO (listeriolizin O)** | LLO alt-birim aşısı *L. monocytogenes*'e karşı koruyucu; immünodominant CD8 hedefi (LLO91–99) | Phelps 2020 — **10.1016/j.vaccine.2020.06.049** |
| **p60 (Iap)** | Sekrete p60'a özgü CD8 T-lenfositleri enfeksiyona karşı korur (p60 217–225) | Harty & Pamer 1995 — **10.4049/jimmunol.154.9.4642** |
| **InlA** | İnternalin; epitel-hücre invazyonunun prototip virülans faktörü + aday antijen | Gaillard 1991 — **10.1016/0092-8674(91)90009-n** |
| **InlB** | İnternalin; hepatosit invazyon virülans faktörü + immünojenik aday antijen | Dramsi 1995 — **10.1111/j.1365-2958.1995.tb02297.x** |
| **ActA** | Aktin-tabanlı hareket faktörü; anti-ActA pasif immünizasyon enfeksiyona karşı korur | Asano 2016 — **10.1038/srep39628** |

> **Dürüst kapsam:** InlA/InlB için kanıt "yerleşik invazyon virülans faktörü + araştırılan aday
> antijen" düzeyindedir (tek-başına ruhsatlı koruyucu aşı değil); ActA koruması **pasif**
> immünizasyon verisidir. Bu, `BENCHMARK_SUMMARY.md` §1b'deki ActA'nın "meşru zor-antijen"
> yorumuyla tutarlıdır.

## 5. *Salmonella enterica* serovar Typhimurium (LT2)

| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **OmpD** | Nontifoidal *Salmonella*'da koruyucu B1b hücre antikor yanıtının kilit hedefi | Gil-Cruz 2009 — **10.1073/pnas.0812431106** |
| **OmpC** | Porin; antikor-aracılı koruyucu bakterisidal bellekten sorumlu ana protein | Secundino 2006 — **10.1111/j.1365-2567.2005.02263.x** |
| **OmpA** | Yüzey porini; koruyucu antikor hedeflemesinde tanımlı (OMP boyutu + LPS O-antijeni) | Domínguez-Medina 2020 — **10.1038/s41467-020-14655-9** |
| **FliC (flagellin)** | Sistemik salmonelloz challenge'ında koruyucu antijen (orta düzey); SseB ile sinerjik | Bumann 2014 — **10.3389/fimmu.2014.00381**; Lee 2017 — **10.4049/jimmunol.1601357** |
| **SseB** | SPI-2 translokon alt-birimi; en umut verici koruyucu antijenlerden; CD4-aracılı koruma | Bumann 2014 — **10.3389/fimmu.2014.00381**; Lee 2017 — **10.4049/jimmunol.1601357** |

---

## Atıf listesi (Methods-hazır, hepsi DOI-doğrulanmış 2026-09-20)

**MenB / Rabies (mevcut — `SELECTION_RATIONALE.md` ile ortak):**
1. Pizza M, et al. Identification of vaccine candidates against serogroup B meningococcus by
   whole-genome sequencing. *Science*. 2000;287(5459):1816–1820. **doi:10.1126/science.287.5459.1816**
2. Giuliani MM, et al. A universal vaccine for serogroup B meningococcus. *PNAS*.
   2006;103(29):10834–10839. **doi:10.1073/pnas.0603940103**
3. Cox JH, Dietzschold B, Schneider LG. Rabies virus glycoprotein. II. *Infect Immun*.
   1977;16(3):754–759. **doi:10.1128/iai.16.3.754-759.1977**

**S. aureus:**
4. Kuklin NA, et al. A novel *Staphylococcus aureus* vaccine: iron surface determinant B induces
   rapid antibody responses in rhesus macaques and specific increased survival in a murine
   *S. aureus* sepsis model. *Infect Immun*. 2006;74(4):2215–2223.
   **doi:10.1128/iai.74.4.2215-2223.2006**
5. Kim HK, DeDent A, Cheng AG, et al. IsdA and IsdB antibodies protect mice against
   *Staphylococcus aureus* abscess formation and lethal challenge. *Vaccine*.
   2010;28(38):6382–6392. **doi:10.1016/j.vaccine.2010.02.097**
6. Josefsson E, Hartford O, O'Brien L, Patti JM, Foster T. Protection against experimental
   *Staphylococcus aureus* arthritis by vaccination with clumping factor A, a novel virulence
   determinant. *J Infect Dis*. 2001;184(12):1572–1580. **doi:10.1086/324430**
7. Adhikari RP, et al. Novel structurally designed vaccine for *S. aureus* α-hemolysin: protection
   against bacteremia and pneumonia. *PLoS One*. 2012;7(6):e38567.
   **doi:10.1371/journal.pone.0038567**
8. Kim HK, Cheng AG, Kim HY, Missiakas DM, Schneewind O. Nontoxigenic protein A vaccine for
   methicillin-resistant *Staphylococcus aureus* infections in mice. *J Exp Med*.
   2010;207(9):1863–1870. **doi:10.1084/jem.20092514**

**L. monocytogenes:**
9. Phelps CC, Vadia S, Boyaka PN, et al. A listeriolysin O subunit vaccine is protective against
   *Listeria monocytogenes*. *Vaccine*. 2020;38(36):5803–5813. **doi:10.1016/j.vaccine.2020.06.049**
10. Harty JT, Pamer EG. CD8 T lymphocytes specific for the secreted p60 antigen protect against
    *Listeria monocytogenes* infection. *J Immunol*. 1995;154(9):4642–4650.
    **doi:10.4049/jimmunol.154.9.4642**
11. Gaillard JL, Berche P, Frehel C, Gouin E, Cossart P. Entry of *L. monocytogenes* into cells is
    mediated by internalin, a repeat protein reminiscent of surface antigens from gram-positive
    cocci. *Cell*. 1991;65(7):1127–1141. **doi:10.1016/0092-8674(91)90009-n**
12. Dramsi S, Biswas I, Maguin E, Braun L, Mastroeni P, Cossart P. Entry of *Listeria
    monocytogenes* into hepatocytes requires expression of InlB, a surface protein of the
    internalin multigene family. *Mol Microbiol*. 1995;16(2):251–261.
    **doi:10.1111/j.1365-2958.1995.tb02297.x**
13. Asano K, et al. Passive immunization with anti-ActA and anti-listeriolysin O antibodies protects
    against *Listeria monocytogenes* infection in mice. *Sci Rep*. 2016;6:39628.
    **doi:10.1038/srep39628**

**S. Typhimurium:**
14. Gil-Cruz C, et al. The porin OmpD from nontyphoidal *Salmonella* is a key target for a
    protective B1b cell antibody response. *PNAS*. 2009;106(24):9803–9808.
    **doi:10.1073/pnas.0812431106**
15. Secundino I, et al. *Salmonella* porins induce a sustained, lifelong specific bactericidal
    antibody memory response. *Immunology*. 2006;117(1):59–70.
    **doi:10.1111/j.1365-2567.2005.02263.x**
16. Domínguez-Medina CC, et al. Outer membrane protein size and LPS O-antigen define protective
    antibody targeting to the *Salmonella* surface. *Nat Commun*. 2020;11:851.
    **doi:10.1038/s41467-020-14655-9**
17. Bumann D. Identification of protective antigens for vaccination against systemic salmonellosis.
    *Front Immunol*. 2014;5:381. **doi:10.3389/fimmu.2014.00381**
18. Lee SJ, Benoun J, Sheridan BS, et al. Dual immunization with SseB/flagellin provides enhanced
    protection against *Salmonella* infection mediated by circulating memory cells. *J Immunol*.
    2017;199(4):1353–1361. **doi:10.4049/jimmunol.1601357**

> İlişkili kayıtlar: `SELECTION_RATIONALE.md`, `BENCHMARK_SUMMARY.md`, `PATHOGEN_ACCESSIONS.md`.
> Çekirdek ilke: her şey literatüre atıflı, uydurma yok ([[vaxforge-literature-rule]]).
</content>

---

# YENİ 6 ZOONOTİK BAKTERİ (2026-10-02; DOI'ler PubMed'den doğrulandı)

**Kaynak:** PubMed (atıflar aşağıda DOI ile). Kanıt gücü antijen-başına açık (koruma-challenge /
aday-antijen / immunodominant). Aksesyonlar `PATHOGEN_ACCESSIONS.md`.

## 6. *Coxiella burnetii* (RSA 493)
| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **Com1** (H7C7D7) | Koruyucu antijen; Com1-yüklü DH adoptif transferi naif farede koruma (Th1/Tc1) | Wang 2011 — **10.1186/1471-2172-12-52** |
| **Mip** (P51752) | Peptidil-prolil izomeraz; immunodominant/aday antijen (Q-ateşi subünit aday paneli) | Wang 2011 — **10.1186/1471-2172-12-52** (bağlam) |
| **OmpH** (Q83DT1) | Dış membran; immunodominant/aday antijen | Wang 2011 — **10.1186/1471-2172-12-52** (bağlam) |
| **GroEL/HtpB** (P19421) | Hsp60; immunodominant ama Wang 2011'de tek-başına koruma göstermedi → zayıf/aday | Wang 2011 — **10.1186/1471-2172-12-52** |
> Dürüstlük: yalnız **Com1** doğrudan koruyucu (adoptif transfer); Mip/OmpH immunodominant aday; GroEL aynı çalışmada koruma veremedi. Q-VAX (tam-hücre) Avustralya'da ruhsatlı; bu proteinler subünit aday.

## 7. *Streptococcus suis* serotip 2 (BM407)
| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **Suilysin (SLY)** (A0A0H3MU02) | Non-hemolitik rSLY(P353L) aday aşı; anti-SLY pasif immünizasyon fareyi ölümden korur | Du 2013 — **10.1016/j.vaccine.2013.07.004** |
| **Enolaz** (A0A0H3N5B2) | Yüzeyde sunulan **koruyucu antijen**; farede tam koruma (moonlighting) | Feng 2009 — **10.1086/644602** |
| **6-fosfoglukonat DH (6PGD)** (A0A0H3N193) | Çok-epitoplu aday aşının (MVHP6) üç immunojen bileşeninden biri | Zhang 2023 — **10.3390/vetsci10060383** |
| **GAPDH** (A0A0H3MSV3) | Yüzey moonlighting adezin; immunojenik/aday antijen | Feng 2009 — **10.1086/644602** (moonlighting bağlam) |
| **DnaJ** (A0A0H3MXM0) | Şaperon; immunojenik/aday antijen | Zhang 2023 — **10.3390/vetsci10060383** (bağlam) |
> SLY/enolaz doğrudan koruyucu; 6PGD/GAPDH/DnaJ aday/moonlighting. Bunlar makalede VaxForge'un lokalizasyon süzgeciyle kaçırdığı sitoplazmik-moonlighting antijenler.

## 8. *Brucella suis* (1330)
| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **Omp25** (Q45689) | Grup-3 OMP; Omp25-mutantları atenüe + Rev.1'e benzer/üstün koruma → aşı adayı | Cloeckaert 2002 — **10.1016/s0378-1135(02)00211-0** |
| **Omp31** (P0A3U5) | Grup-3 OMP; *B. ovis* fare modelinde **koruyucu antijen** (immunodominant) | Cloeckaert 2002 — **10.1016/s0378-1135(02)00211-0** |
| **BP26/Omp28** (P0A3U9) | İmmunodominant; tanı markörü + brusellozda aşılama için kullanılır | Kim 2013 — **10.1016/j.jmb.2013.01.015** |
| **SodC (Cu-Zn SOD)** (P66827) | DNA aşısı BALB/c'de **koruyucu** (RB51 düzeyinde), Th1 | Oñate 2003 — **10.1128/IAI.71.9.4857-4861.2003** |
| **L7/L12 (RplL)** (P0A469) | Ribozomal; divalan L7/L12 DNA aşısı virülan 544 challenge'a karşı **koruyucu** | Luo 2006 — **10.1128/IAI.74.5.2734-2741.2006** |
> RB51/S19/Rev.1 ruhsatlı **veteriner** canlı-atenüe aşılar (tam-hücre); bu subünitler koruyucu-aday.

## 9. *Mycobacterium bovis* (AF2122/97)
| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **ESAT-6 (EsxA)** (P0A565) | Subünit aşı (MPL+DDA) **BCG düzeyinde koruyucu** T-hücre yanıtı | Brandt 2000 — **10.1128/IAI.68.2.791-795.2000** |
| **Ag85B (FbpB)** (P0C2T2) | Ag85B-ESAT-6 füzyon subünit aşısı hayvan modelinde **koruyucu** | Davila 2010 — **10.1155/2010/208371** |
| **CFP-10 (EsxB)** (P0A567) | ESAT-6 ile RD1 ikilisi; koruyucu/immunodominant aday (Ag85B-ESAT6 füzyon bağlamı) | Davila 2010 — **10.1155/2010/208371** |
| **MPB64** (P0A5Q5) | *M. bovis* salgı; immunodominant/tanısal + aday antijen | Brandt 2000 — **10.1128/IAI.68.2.791-795.2000** (bağlam) |
| **HspX/Acr** (P0A5B8) | Latensi α-kristalin; immunodominant/aday antijen | Davila 2010 — **10.1155/2010/208371** (bağlam) |
> BCG (canlı *M. bovis*) ruhsatlı TB aşısı; ESAT-6/Ag85B bunlar klinik subünit adayları (H1/H56).

## 10. *Leptospira interrogans* serovar Lai (56601)
| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **LigA** (Q8EYU4) | LigA DNA aşısı hamster'da virülan challenge'a karşı **anlamlı koruma** | Faisal 2007 — **10.1016/j.vaccine.2007.10.029** |
| **LigB** (Q8EZS3) | Lig ailesi (Ig-benzeri tekrar); koruyucu-aday (LigA ile aynı aile) | Faisal 2007 — **10.1016/j.vaccine.2007.10.029** (aile bağlamı) |
| **LipL32** (O34094) | Rekombinant LipL32-41-OmpL1 **koruyucu antijen** DNA/protein aşısı | Feng 2009 — **10.1590/s0100-879x2009005000013** |
| **LipL41** (Q8F8E1) | Çok-bileşenli LipL32-41-OmpL1 aşısında anti-LipL41 koruyucu yanıt | Feng 2009 — **10.1590/s0100-879x2009005000013** |
| **OmpL1** (G1UB30) | Çok-bileşenli LipL32-41-OmpL1 aşısında anti-OmpL1 koruyucu yanıt | Feng 2009 — **10.1590/s0100-879x2009005000013** |
> Leptospiroz'da tam-hücre bakterin veteriner aşılar ruhsatlı; bu subünitler koruyucu-aday.

## 11. *Bacillus anthracis* (Ames Ancestor)
| Antijen | Rol / kanıt | Atıf (DOI-doğrulanmış) |
|---|---|---|
| **PA (pagA)** (P13423) | **RUHSATLI** ABD (AVA/BioThrax) + UK (AVP) insan antraks aşılarının ana koruyucu bileşeni | Rhie 2005 — **10.1016/j.femsim.2005.05.009**; Baillie 2010 — **10.1016/j.vaccine.2010.07.075** |
| **LF (lef)** (P15917) | LF/LF-D1 tam **koruma** sağlar (B.anthracis spor challenge); AVP koruyuculuğuna katkı | Baillie 2010 — **10.1016/j.vaccine.2010.07.075** |
| **EF (cya)** (P40136) | Ödem toksini bileşeni; toksin-nötralizasyon hedefi/aday (PA+LF ile tripartit) | Baillie 2010 — **10.1016/j.vaccine.2010.07.075** (bağlam) |
| **CapB / γ-PGA kapsül** (P19580) | γ-DPGA-protein konjugatı koruyucu IgG + opsonofagositoz; antraks aşısına ek | Schneerson 2003 — **10.1073/pnas.1633512100** |
| **S-tabaka EA1 (eag)** (P94217) | Yüzey S-tabaka; plazmit-EA1 kalıcı tip-1/2 yanıt + T-epitop; PA-dışı aday antraks antijeni | Zhang 2007 — **10.1016/j.vaccine.2007.11.072** |
> **PA ruhsatlı insan aşısı bileşeni** (BioThrax/AVA, AVP). LF/kapsül koruyucu-ek; EF/EA1 aday.
> **EA1 atıfı eklendi (2026-10-03):** Zhang 2007, plazmit-EA1 kalıcı bağışık yanıt + T-epitop.

## Atıf listesi — yeni 6 bakteri (Methods-hazır, DOI-doğrulanmış 2026-10-02, PubMed)
14. Wang Y, et al. Efficient activation of T cells by HMDCs pulsed with *Coxiella burnetii* Com1. *BMC Immunol*. 2011;12:52. **doi:10.1186/1471-2172-12-52**
15. Du H, et al. The genetically modified suilysin rSLY(P353L) provides a candidate vaccine ... *Streptococcus suis*. *Vaccine*. 2013;31(38):4209-15. **doi:10.1016/j.vaccine.2013.07.004**
16. Feng Y, et al. *Streptococcus suis* enolase functions as a protective antigen displayed on the bacterial cell surface. *J Infect Dis*. 2009;200(10):1583-92. **doi:10.1086/644602**
17. Zhang Y, et al. Development of a universal multi-epitope vaccine candidate against *S. suis* (HP0197/PepO/6PGD). *Vet Sci*. 2023;10(6):383. **doi:10.3390/vetsci10060383**
18. Cloeckaert A, et al. Major outer membrane proteins of *Brucella* spp.: past, present and future. *Vet Microbiol*. 2002;90(1-4):229-47. **doi:10.1016/s0378-1135(02)00211-0**
19. Oñate AA, et al. A DNA vaccine encoding Cu,Zn SOD of *Brucella abortus* induces protective immunity. *Infect Immun*. 2003;71(9):4857-61. **doi:10.1128/IAI.71.9.4857-4861.2003**
20. Luo D, et al. Protective immunity by a divalent DNA vaccine encoding L7/L12 and Omp16 of *Brucella abortus*. *Infect Immun*. 2006;74(5):2734-41. **doi:10.1128/IAI.74.5.2734-2741.2006**
21. Kim D, et al. *Brucella* immunogenic BP26 forms a channel-like structure. *J Mol Biol*. 2013;425(7):1119-26. **doi:10.1016/j.jmb.2013.01.015**
22. Brandt L, et al. ESAT-6 subunit vaccination against *Mycobacterium tuberculosis*. *Infect Immun*. 2000;68(2):791-5. **doi:10.1128/IAI.68.2.791-795.2000**
23. Davila J, et al. Genetic diversity of esxA/esxH/fbpB and implications for Ag85B-ESAT-6 / Ag85B-TB10.4 subunit vaccines. *J Biomed Biotechnol*. 2010;2010:208371. **doi:10.1155/2010/208371**
24. Faisal SM, et al. Evaluation of protective immunity of *Leptospira* LigA DNA vaccine in hamsters. *Vaccine*. 2007;26(2):277-87. **doi:10.1016/j.vaccine.2007.10.029**
25. Feng CY, et al. Immune strategies using LipL32 and recombinant LipL32-41-OmpL1 vaccines against *Leptospira*. *Braz J Med Biol Res*. 2009;42(9):796-803. **doi:10.1590/s0100-879x2009005000013**
26. Baillie LW, et al. An anthrax subunit vaccine candidate based on protective regions of PA and LF. *Vaccine*. 2010;28(41):6740-8. **doi:10.1016/j.vaccine.2010.07.075**
27. Rhie GE, et al. Efficacy of non-toxic deletion mutants of protective antigen from *Bacillus anthracis*. *FEMS Immunol Med Microbiol*. 2005;45(2):341-7. **doi:10.1016/j.femsim.2005.05.009**
28. Schneerson R, et al. Poly(γ-D-glutamic acid)-protein conjugates induce IgG to the capsule of *Bacillus anthracis*. *PNAS*. 2003;100(15):8945-50. **doi:10.1073/pnas.1633512100**
29. Zhang Y, et al. Plasmid-based vaccination with candidate anthrax vaccine antigens (EA1/SAP/capA) induces durable Th1/Th2 responses. *Vaccine*. 2007;26(5):614-22. **doi:10.1016/j.vaccine.2007.11.072**
