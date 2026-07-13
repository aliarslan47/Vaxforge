#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# VaxForge Guardian — OTONOM SÜPERVİZÖR (Claude olmadan çalışır)
#
# Backend (uvicorn) + Frontend (next start) süreçlerini başlatır, sağlığı
# izler, biri çökerse: loglar → INCIDENT.json yazar → bildirir → yeniden
# başlatır. "Durduk yere hata alırsa sistemi ayağa kaldırsın" katmanı budur.
#
# Kullanım:  web/guardian/supervise.sh
# Durdur:    Ctrl-C  (tüm çocuk süreçler temizlenir)
# ═══════════════════════════════════════════════════════════════════════════
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BE="$ROOT/web/backend"
FE="$ROOT/web/frontend"
VENV="$ROOT/.venv"
GDIR="$ROOT/web/guardian"
LOGS="$GDIR/logs"
INC="$GDIR/incidents"
mkdir -p "$LOGS" "$INC"

BACKEND_PORT="${VF_BACKEND_PORT:-8011}"
FRONTEND_PORT="${VF_FRONTEND_PORT:-3000}"
export VF_BACKEND_PORT="$BACKEND_PORT" VF_FRONTEND_PORT="$FRONTEND_PORT"

PY="python3"; [ -x "$VENV/bin/python" ] && PY="$VENV/bin/python"

BACKEND_PID=""; FRONTEND_PID=""

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOGS/supervisor.log"; }

notify() {
  # Masaüstü bildirimi (varsa) + kalıcı bayrak dosyası (Claude hook/oturum görür).
  local title="$1" body="$2"
  command -v notify-send >/dev/null 2>&1 && notify-send -u critical "$title" "$body" 2>/dev/null || true
  printf '%s | %s\n' "$title" "$body" >> "$GDIR/ALERTS.log"
}

incident() {
  # service, reason -> incidents/INCIDENT.json (guardian subagent bunu okur)
  local service="$1" reason="$2" tail_log="$3"
  local ts; ts="$(date -u '+%FT%TZ')"
  cat > "$INC/INCIDENT.json" <<JSON
{
  "timestamp": "$ts",
  "service": "$service",
  "reason": "$reason",
  "backend_port": $BACKEND_PORT,
  "frontend_port": $FRONTEND_PORT,
  "log_tail": $(printf '%s' "$tail_log" | tail -20 | python3 -c 'import sys,json; print(json.dumps(sys.stdin.read()))' 2>/dev/null || echo '""'),
  "status": "open"
}
JSON
  cp "$INC/INCIDENT.json" "$INC/incident_${ts//[:]/-}.json" 2>/dev/null || true
}

start_backend() {
  log "backend başlatılıyor (uvicorn :$BACKEND_PORT)"
  (cd "$BE" && exec "$VENV/bin/uvicorn" main:app --host 127.0.0.1 --port "$BACKEND_PORT" \
      --log-level warning) >>"$LOGS/backend.log" 2>&1 &
  BACKEND_PID=$!
}

start_frontend() {
  log "frontend başlatılıyor (next start :$FRONTEND_PORT)"
  (cd "$FE" && VAXFORGE_API="http://127.0.0.1:$BACKEND_PORT" exec npm run start) \
      >>"$LOGS/frontend.log" 2>&1 &
  FRONTEND_PID=$!
}

cleanup() {
  log "kapatılıyor — çocuk süreçler durduruluyor"
  [ -n "$BACKEND_PID" ] && kill "$BACKEND_PID" 2>/dev/null || true
  [ -n "$FRONTEND_PID" ] && kill "$FRONTEND_PID" 2>/dev/null || true
  exit 0
}
trap cleanup INT TERM

# --- başlangıç ---
log "════ VaxForge Guardian başladı ════"
start_backend
start_frontend

# servislerin ayağa kalkması için ısınma
sleep 4

declare -A FAILS=([backend]=0 [frontend]=0)
MAX_RESTARTS="${VF_MAX_RESTARTS:-5}"

while true; do
  sleep 10

  # süreç ölmüş mü? (PID kontrol)
  if ! kill -0 "$BACKEND_PID" 2>/dev/null; then
    FAILS[backend]=$(( FAILS[backend] + 1 ))
    log "⚠ backend süreci öldü (deneme ${FAILS[backend]}/$MAX_RESTARTS)"
    incident "backend" "process_exited" "$(tail -30 "$LOGS/backend.log" 2>/dev/null)"
    notify "VaxForge backend çöktü" "Yeniden başlatılıyor (:$BACKEND_PORT)"
    if [ "${FAILS[backend]}" -le "$MAX_RESTARTS" ]; then start_backend; sleep 3
    else log "backend $MAX_RESTARTS kez çöktü — otomatik restart durduruldu (kod fix gerek)"; fi
  fi

  if ! kill -0 "$FRONTEND_PID" 2>/dev/null; then
    FAILS[frontend]=$(( FAILS[frontend] + 1 ))
    log "⚠ frontend süreci öldü (deneme ${FAILS[frontend]}/$MAX_RESTARTS)"
    incident "frontend" "process_exited" "$(tail -30 "$LOGS/frontend.log" 2>/dev/null)"
    notify "VaxForge frontend çöktü" "Yeniden başlatılıyor (:$FRONTEND_PORT)"
    if [ "${FAILS[frontend]}" -le "$MAX_RESTARTS" ]; then start_frontend; sleep 3
    else log "frontend $MAX_RESTARTS kez çöktü — otomatik restart durduruldu (kod fix gerek)"; fi
  fi

  # süreç canlı ama HTTP yanıt vermiyor mu? (asılı kalma)
  hc="$("$GDIR/healthcheck.sh" 2>&1)"
  if echo "$hc" | grep -q "DOWN backend"; then
    log "⚠ backend HTTP yanıt vermiyor (süreç canlı ama asılı)"
  fi
  if echo "$hc" | grep -q "DOWN frontend"; then
    log "⚠ frontend HTTP yanıt vermiyor (süreç canlı ama asılı)"
  fi
done
