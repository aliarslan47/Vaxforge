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
| 0 | Güvenlik ağı + kayıt (smoke test, ROADMAP) | #10-kısmi | 🔨 | — |
| 1 | PSORTb lokalizasyon (bakteri, Gram) | #7 | ⬜ | — |
| 2 | DeepLoc2 lokalizasyon (ökaryot/virüs) | #8 | ⬜ | — |
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
- **2026-09-26**: Faz 0 başladı. `tests/test_smoke.py` yazıldı (3-protein demo ile uçtan
  uca pipeline koşuyor). Smoke doğrulanınca baz commit + Faz 1'e geçilecek.
