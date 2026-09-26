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
| 3 | ESMFold yapı + konformasyonel B + RSA | #1 | ⬜ | — |
| 4 | Glikozilasyon + glikan-maskeli epitop | #2 | ⬜ | — |
| 5 | Kodon optimizasyonu + in silico klonlama | #3 | ⬜ | — |
| 6 | VaxiJen antijenite çapraz-kontrol | #4 | ⬜ | — |
| 7 | BLOSUM tabanlı kısmi mimikri | #5 | ⬜ | — |
| 8 | Otomatik çok-suş konservasyonu | #6 | ⬜ | — |
| 9 | Otoimmünite konak-DB genişletme | #9 | ⬜ | — |
| 10 | Türe-özgü ImmForge kalibrasyonu | #11 | ⬜ | — |
| 11 | Kapsamlı birim test paketi | #10 | ⬜ | — |

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
- **Sıradaki: Faz 3 (ESMFold yapı)** — YEREL ESMFold modeli önbellekte YOK + GPU yok → yerel
  folding pratik değil. Plan: **ESMAtlas fold API** (`api.esmatlas.com/foldSequence/v1/pdb/`,
  ≤400aa, anahtarsız) → PDB; RSA (Shrake-Rupley) + konformasyonel-B (yüzey/RSA tabanlı).
  Offline/uzun dizi → dürüst ARAÇ_YOK; BepiPred lineer korunur.
