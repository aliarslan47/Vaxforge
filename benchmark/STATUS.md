# VaxForge yeni-virüs benchmark — DURUM

- Son güncelleme: **2026-10-18 11:17:19**
- Faz: ****
- Son olay: **Faz 2: listeria bitti**
- Detay: Recall@antijen: 5/5  ·  fold: 3.83  ·  Fisher p: 0.00123 (rc=0)

## Son 10 olay

- `2026-10-10 12:53:09` [] örnek veri demosu — demo proteom (26 nihai aday): ESKİ ToxinPred2 hibrit t=0.6 sert kapı adayların %61.5'ini (16/26) elerdi; YENİ ToxinPred3 ET yumuşak kapı 0 işaret / 0 eleme. benchmark/results/demo_toxicity_gate.json
- `2026-10-10 13:34:36` [] ssuis kalibrasyon-sonrası koşu başladı — kalibre kapılarla (ToxinPred3 yumuşak + Gram-duyarlı lokalizasyon) tam proteom; eski sonuçlar ssuis_headtohead.pre_calibration.json olarak yedeklendi
- `2026-10-10 13:51:25` [] ssuis kalibrasyon-sonrası koşu — 3 düzeltme aktif: ToxinPred3 yumuşak toksisite + Gram-duyarlı lokalizasyon + otoimmünite granülaritesi (protein yumuşak/epitop sert). Eski sonuçlar ssuis_headtohead.pre_calibration.json
- `2026-10-10 16:29:48` [] ssuis kalibrasyon-sonrası koşu BİTTİ — 9194s. recall 1/5 (değişmedi), fold 1.34→0.79, huni 474→490, aday 1474→4872, mimikri 16 tam eşleşme eledi. Kaçan 4 antijenin bağlayıcı kısıtı lokalizasyon: PSORTb dördünü de sitoplazmik veriyor (9.97/9.97/7.5/9.97) — moonlighting/klasik-olmayan salgı, dizi-temelli tahminciyle erişilemez. KISITLAMA olarak raporlanacak.
- `2026-10-12 00:18:02` [] askıda 3 yayın kapatıldı + 4. kayıp mekanizması — V5 (doi:10.1093/infdis/jiz352) ve SBP2' (doi:10.1016/j.rvsc.2021.05.006) açık erişimli ikincil makalelerden ÇÖZÜLDÜ (eLife PMC12185102; AEM PMC13599487) - eşleme ada göre değil ölçümle: V5'te uzunluk-birebir (SrtA 249, SLY 497) + ScpA hizalamasıyla ScpCL/SspA ayrımı; SBP2'=ZY05719_RS11000 dört kontrolle (srtBCD komşuluk + tek YPKN + RrgB N-term %31.1 vs yayın %29 + 33 kDa). HK/Gtf-2/PsbP ÇÖZÜLEMEDİ (paywall; belirsizlik sayıldı: HK 10, Gtf 6 eşleşme). Küratörlük 9->15, recall 1/4 DEĞİŞMEDİ. YENİ: SspA'nın kaybı lokalizasyon değil TM kapısı (2 TM>1, cell_wall 9.98 izinli) = 4. mekanizma; kapsam ölçüldü (antijenlerde 1/16), kapı DEĞİŞTİRİLMEDİ, kalibrasyon kuyruğuna alındı. Atıf düzeltmesi: enolaz makalesi Bichet değil Dumesnil ve ark. Makale EN+TR güncellendi + 4 referans.
- `2026-10-18 00:15:19` [] Faz 2 başladı — 10 tür sıralı, tek kod sürümü; TM sinyal-duyarlı + toksisite + otoimmünite kalibre
- `2026-10-18 02:04:29` [] Faz 2: coxiella bitti — Recall@antijen: 3/4  ·  fold: 2.37  ·  Fisher p: 0.0973 (rc=0)
- `2026-10-18 04:28:52` [] Faz 2: ssuis bitti — Recall@antijen: 1/5  ·  fold: 0.77  ·  Fisher p: 0.777 (rc=0)
- `2026-10-18 07:54:35` [] Faz 2: saureus bitti — Recall@antijen: 5/5  ·  fold: 2.81  ·  Fisher p: 0.00578 (rc=0)
- `2026-10-18 11:17:19` [] Faz 2: listeria bitti — Recall@antijen: 5/5  ·  fold: 3.83  ·  Fisher p: 0.00123 (rc=0)
