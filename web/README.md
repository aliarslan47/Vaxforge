# VaxForge Web

Mevcut `vaxforge/` pipeline'ı üzerine profesyonel web uygulaması.

- **backend/** — FastAPI, `vaxforge.pipeline.run()`'ı SSE ile web'e açar (`:8011`)
- **frontend/** — Next.js 14 + Tailwind, landing + pipeline çalıştırıcı (`:3000`)
- **guardian/** — kendi kendini iyileştirme: süpervizör + sağlık/kod kontrolü

## Çalıştırma (önerilen — guardian süpervizör)

Tek komut hem backend'i hem frontend'i başlatır ve **çökerse otomatik ayağa kaldırır**:

```bash
cd ~/vaxforge
bash web/guardian/supervise.sh
```

Aç: **http://localhost:3000**

Süpervizör: bir servis çökerse `web/guardian/logs/`'a yazar, `web/guardian/incidents/INCIDENT.json` + `web/guardian/ALERTS.log` oluşturur, masaüstü bildirimi (varsa) atar ve süreci yeniden başlatır (5 denemeye kadar).

## Kalıcı otonomi (makine açılışında, Claude olmadan)

`systemd` kullanıcı servisi olarak:

```bash
mkdir -p ~/.config/systemd/user
cat > ~/.config/systemd/user/vaxforge.service <<'UNIT'
[Unit]
Description=VaxForge Web (guardian supervisor)
After=network.target

[Service]
Type=simple
WorkingDirectory=/home/msi-nb/vaxforge
ExecStart=/usr/bin/bash /home/msi-nb/vaxforge/web/guardian/supervise.sh
Restart=always
RestartSec=5

[Install]
WantedBy=default.target
UNIT

systemctl --user daemon-reload
systemctl --user enable --now vaxforge.service
loginctl enable-linger msi-nb   # oturum kapalıyken de çalışsın
```

## Manuel (guardian olmadan, geliştirme)

```bash
# backend
cd ~/vaxforge/web/backend && ../../.venv/bin/uvicorn main:app --port 8011
# frontend (ayrı terminal)
cd ~/vaxforge/web/frontend && npm run dev
```

## Kod-seviyesi kendini iyileştirme (Claude)

- `.claude/settings.json`'daki **PostToolUse hook**'u, `web/**` altında her düzenlemeden sonra
  `web/guardian/codecheck.sh` (TS tip + Python import) çalıştırır; hata varsa Claude'a geri besler.
- Çökme/hata olduğunda **`site-guardian`** subagent'ı görevlendirilir (`~/.claude/agents/site-guardian.md`):
  INCIDENT + logları okur, kök nedeni bulur, `web/` katmanını düzeltir, doğrular, ayağa kaldırır.
  (Not: hook değişikliği yeni Claude oturumunda aktifleşir.)

## Görsel doğrulama (dev aracı)

```bash
cd ~/vaxforge/web/frontend && node shot.mjs   # /tmp/vfshots/*.png
```
