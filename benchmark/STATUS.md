# VaxForge yeni-virüs benchmark — DURUM

- Son güncelleme: **2026-10-10 12:53:09**
- Faz: ****
- Son olay: **örnek veri demosu**
- Detay: demo proteom (26 nihai aday): ESKİ ToxinPred2 hibrit t=0.6 sert kapı adayların %61.5'ini (16/26) elerdi; YENİ ToxinPred3 ET yumuşak kapı 0 işaret / 0 eleme. benchmark/results/demo_toxicity_gate.json

## Son 10 olay

- `2026-10-01 00:21:53` [done] vaxignml TAMAM — 8/9 · vsv:#3/5 | emcv:#1/11 | rvfv:#1/4 | veev:#5/10 | reovirus:#1/11 | jev:#3/15 | vaccinia:#3/223 | cowpox:#3/233
- `2026-10-02 19:32:08` [done] new6 benchmark TAMAM — 10-bakteri head-to-head (genom sırası) — VaxForge 37/48(77%) · Vaxign-ML 47/48(98%) · NERVE2 14/48(29%). Yeni-6 recall: coxiella 2/4, ssuis 1/5, brucella 4/5, mbovis 4/5, leptospira 4/5, banthracis 4/5. Vaxign+NERVE 6/6 lokal podman.
- `2026-10-02 19:32:08` [fix] netMHCpan deadlock fix — netmhc_local.py: izole TMPDIR+stdin=DEVNULL+timeout 7200->1200; brucella kilidini cozdu, 4/5 ile yeniden kostu
- `2026-10-04 12:57:10` [analysis] analysis — kacan-antijen atfi TAMAM — 48 antijen, 11 kayip tamamen atfedildi: 6 konak-homoloji sert filtresi (insan %37-52; Fisher p=3.8e-05, yakalananlarda 0/37), 2 Gram(+) 'membrane' yapilandirma acigi (CapB/ActA), 1 PSORTb yanlis atama (LigA), 2 epitop-asamasi (ESAT-6/L7-L12; toksisite hipotezi curutuldu). OLCEK-BUG YOK. Belge: MISS_ATTRIBUTION.md, veri: results/miss_attribution.json
- `2026-10-10 01:49:37` [] makale yeniden tasarımı + Gram(+) fix — MANUSCRIPT EN+TR Nature kalıbına göre baştan yazıldı (Giriş-Sonuçlar-Tartışma-Yöntemler); kaçan-antijen atfı ana sonuç bölümü oldu (Tablo 4); 3 metin hatası düzeltildi (C.burnetii/S.suis kayıpları lokalizasyon degil konak-homoloji, L7/L12 epitop asamasi, ActA Gram+ acigi); yeni nesil konumlandirma eklendi (PLGDL Nat Commun 2025, Vaxign-DL, Podda ESM-2); Gram-duyarli izinli liste kodda uygulandi + test; 15/15 test yesil
- `2026-10-10 01:52:07` [] Gram(+) re-run başladı — 5 tür sıralı: ssuis listeria saureus mbovis banthracis; allowed_gram_positive aktif
- `2026-10-10 04:13:07` [] Gram(+) re-run: ssuis bitti — Recall@antijen: 1/5  ·  fold: 0.91  ·  Fisher p: 0.711 (rc=0)
- `2026-10-10 07:19:40` [] Gram(+) re-run: listeria bitti — Recall@antijen: 3/5  ·  fold: 3.18  ·  Fisher p: 0.0498 (rc=0)
- `2026-10-10 10:51:47` [] Filtre kalibrasyonu — Toksisite kapısı yanlış uzunluk sınıfında: ToxinPred2 t=0.6 insan self 9-mer %75.2, doğrulanmış IEDB epitop %63.1 yanlış-pozitif (GILGFVFTL dahil); tam proteinde %94.3/%4.7 saglam. Duzeltme: epitopta ToxinPred3 YUMUSAK + MEV'de ToxinPred2 sert; alerjenite konstruktta FAO/WHO 80aa/%35. Kalici kalibrasyon duzenegi + 2 regresyon testi; 15/15 yesil. Gram(+) yeniden kosumu durduruldu.
- `2026-10-10 12:53:09` [] örnek veri demosu — demo proteom (26 nihai aday): ESKİ ToxinPred2 hibrit t=0.6 sert kapı adayların %61.5'ini (16/26) elerdi; YENİ ToxinPred3 ET yumuşak kapı 0 işaret / 0 eleme. benchmark/results/demo_toxicity_gate.json
