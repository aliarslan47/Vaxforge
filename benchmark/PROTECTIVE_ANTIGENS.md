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
