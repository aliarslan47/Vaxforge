# VaxForge Benchmark — Otonom İlerleme Kaydı

Kullanıcı uyurken otonom yürütme (2026-08-19 gece). Her adımda buraya + belleğe kayıt.
Sıra: (1) Vaxign-ML rabies → (2) genişletilmiş 204-set MenB 3'lü flagship → (3) rakipler tam
2003 proteom → (4) figür/tablo. Soru sorulmayacak.

## TAMAMLANAN (bu oturum)
- [x] Vaxign-ML kuruldu (podman, docker.io/e4ong1031/vaxign-ml)
- [x] NERVE 2.0 kuruldu (nerve:v0.0.8 + psortb + nerve-network)
- [x] Patojen aksesyonları doğrulandı → `PATHOGEN_ACCESSIONS.md`
- [x] `stats.py` (Fisher/permütasyon/AUC-bootstrap/Jaccard) + self-test
- [x] VİRÜS: rabies × 4-konak VaxForge → `results/rabies/` (G rank #2/5; çok-konak Jaccard)
- [x] BAKTERİ: MenB 3'lü head-to-head (64-set, default) → `results/menb_headtohead.json`
      VaxForge 4/4 AUC1.00 · Vaxign-ML 4/4 AUC0.98 · NERVE 1/4 (3 antijen elendi)

## SIRADAKİ (otonom, sırayla)
- [x] ADIM 1: Vaxign-ML rabies (virus) → G rank #4/5 (pred 1.0), skorlar düz (virüste zayıf kalibre).
      VaxForge G #2/5 + çok-konak. → `results/vaxignml_rabies_summary.json`. NERVE rabies'te N/A (bacterial).
- [~] ADIM 2 BAŞLATILDI: 204-set (menb_testset_204.faa) — VaxForge (PID koşuyor, funnel 71 geçti,
      epitop fazı ~2-3h; outputs/validation_menb_204/) + Vaxign-ML (results/vaxignml_menb204/) +
      NERVE (tools/NERVE/out_menb204/). validate_menb.py env-parametrize edildi (MENB_TESTSET/OUTDIR/N_BACKGROUND).
      Bitince compare_menb.py'nin 204 varyantıyla head-to-head + stats.
- [ ] ADIM 2 (orijinal): genişletilmiş MenB seti (4 koruyucu + 200 arka plan = 204, SEED=42);
      3 aracı da bu sette koş → daha sağlam fold-enrichment; fair (özdeş girdi).
      VaxForge ~2-3h (CPU) = gecelik iş; Vaxign-ML/NERVE dakikalar.
      NOT: tam 2003-proteom VaxForge CPU'da ~17h+ → pratik değil; 204-set pragmatik+adil seçim.
- [x] ADIM 2 TAMAM: 204-set flagship head-to-head → `results/menb204_headtohead.json`.
      VaxForge 4/4 AUC0.99 fold3.12 p0.012 · Vaxign-ML 4/4 AUC0.97 fold3.64 p0.0065 ·
      NERVE 1/4 (3 antijen elendi) fold25 p0.058(anlamsız). VaxForge-204 süre ~1h55m.
      64-set kalıbı 204'te SAĞLAM tekrarlandı. AUC concordance formülü eklendi (64→1.00 doğrulandı).
- [x] ADIM 3 TAMAM: rakipler tam 2003-proteom → `results/menb_full_competitors.json`.
      Vaxign-ML tam: recall 4/4, 613 aday (gevşek), fold 3.28, ranks PorA#29/NadA#52/NHBA#59/fHbp#150.
      NERVE tam: recall 1/4 (yalnız PorA#8/9), 9 aday, fold 62.5 (aşırı-seçicilik → yüksek fold ama recall feda).
      ID formatı sp|ACC| (testset'te çıplaktı) — parse düzeltildi.
      KİLİT: NERVE HER ÖLÇEKTE (64/204/2003) 3/4 antijen kaçırıyor; yüksek fold = aşırı-seçicilik, recall değil.
- [x] ADIM 4 TAMAM: 4 figür üretildi + görsel doğrulandı → `figures/` (fig1 antijen-kurtarma ★,
      fig2 recall-bar, fig3 fold-vs-recall, fig4 çok-konak Jaccard ★). Script `make_figures.py`.
      Master özet: `BENCHMARK_SUMMARY.md`.

## ✅ TÜM ADIMLAR TAMAM (2026-08-20 ~02:25). Kullanıcı uyanınca hazır.
Kalan (yazım aşaması, kullanıcı kuralı): patojen/konak seçilim gerekçeleri DOI'li; VaxForge
tam-2003 flagship (GPU gelince). Kozmetik: fig3 alt-sağ etiket örtüşmesi.

## KARARLAR (otonom, kullanıcı uyurken alındı)
- Flagship head-to-head = 204-protein genişletilmiş set (tam 2003 VaxForge CPU'da imkânsız);
  hepsi aynı sette = adil. Tam proteom yalnız hızlı rakiplerde bonus olarak.
- Asimetri (VaxForge subsample vs rakip tam-proteom) yazıda AÇIKÇA belirtilecek.
- Patojen/konak seçilim gerekçeleri DOI'li — YAZIM aşaması (kullanıcı kuralı).

## ✅ BAKTERİYEL KOL GENİŞLETMESİ TAMAM (2026-08-27 01:09) — 3 bakteri, TAM PROTEOM
Kullanıcı isteği: 4-konağı enfekte eden 3 bakteri, genom-boyu artan, TAM proteom, tek tek, rakip-önce.
Seçilen: S. aureus (2889) · L. monocytogenes (2844) · S. Typhimurium (4533).
Yeni kod: `bacteria/*.json`, `prep_proteome.py`, `run_bacterium.py`, `run_competitors.sh`,
`compare_bacterium.py`, `run_bug_full.sh`, `rerun_vaxforge_all.sh`.
- İlk koşularda VaxForge yapay düşük çıktı (saureus 3/5, listeria 3/5, salmonella 0/5) → **iki ölçek-bug'ı**
  bulundu+düzeltildi+doğrulandı: (A) PSORTb tek-konteyner timeout → 800'lük chunk; (B) netMHCpan tek-dev-batch
  (1.37M peptit → %0) → 20k-chunk + 16-çekirdek paralel. Toksisite/alerjen suçsuzdu (kanıtlandı).
- **NİHAİ (2-fix, tam proteom, default): VaxForge 14/15** (saureus 5/5 AUC1.00 · listeria 4/5 · salmonella 5/5 AUC0.99)
  · Vaxign-ML 15/15 · NERVE 9/15. VaxForge ≈ Vaxign-ML recall + AUC/rank önde; NERVE geride. Tek kayıp ActA (meşru).
- Sonuçlar: `results/{saureus,listeria,salmonella}_headtohead.json`; figür `figures/fig5_bacteria_panel.png`;
  master `BENCHMARK_SUMMARY.md §1b + §metodoloji`; `PATHOGEN_ACCESSIONS.md` + `SELECTION_RATIONALE.md` güncel.
- KALAN (yazım): koruyucu-antijen "koruyucu" atıfları DOI'li; iki-bug fix Methods/validation'a girecek.
