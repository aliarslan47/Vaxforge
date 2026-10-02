# VaxForge yeni-virüs benchmark — DURUM

- Son güncelleme: **2026-10-02 19:32:08**
- Faz: **fix**
- Son olay: **netMHCpan deadlock fix**
- Detay: netmhc_local.py: izole TMPDIR+stdin=DEVNULL+timeout 7200->1200; brucella kilidini cozdu, 4/5 ile yeniden kostu

## Son 10 olay

- `2026-10-01 00:20:24` [vaxignml] vaxignml bitti — reovirus: sigma-1 rank #1/11 (23s)
- `2026-10-01 00:20:24` [vaxignml] vaxignml virüs — jev başladı
- `2026-10-01 00:20:47` [vaxignml] vaxignml bitti — jev: envelope protein E rank #3/15 (23s)
- `2026-10-01 00:20:47` [vaxignml] vaxignml virüs — vaccinia başladı
- `2026-10-01 00:21:19` [vaxignml] vaxignml bitti — vaccinia: L1R rank #3/223 (32s)
- `2026-10-01 00:21:19` [vaxignml] vaxignml virüs — cowpox başladı
- `2026-10-01 00:21:53` [vaxignml] vaxignml bitti — cowpox: L1R ortholog (CPXV099) rank #3/233 (34s)
- `2026-10-01 00:21:53` [done] vaxignml TAMAM — 8/9 · vsv:#3/5 | emcv:#1/11 | rvfv:#1/4 | veev:#5/10 | reovirus:#1/11 | jev:#3/15 | vaccinia:#3/223 | cowpox:#3/233
- `2026-10-02 19:32:08` [done] new6 benchmark TAMAM — 10-bakteri head-to-head (genom sırası) — VaxForge 37/48(77%) · Vaxign-ML 47/48(98%) · NERVE2 14/48(29%). Yeni-6 recall: coxiella 2/4, ssuis 1/5, brucella 4/5, mbovis 4/5, leptospira 4/5, banthracis 4/5. Vaxign+NERVE 6/6 lokal podman.
- `2026-10-02 19:32:08` [fix] netMHCpan deadlock fix — netmhc_local.py: izole TMPDIR+stdin=DEVNULL+timeout 7200->1200; brucella kilidini cozdu, 4/5 ile yeniden kostu
