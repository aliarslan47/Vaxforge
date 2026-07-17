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

# Bir portu tutan süreçleri bul ve öldür, port boşalana kadar bekle.
# KÖK NEDEN FİX: restart öncesi eski süreç portu bırakmazsa yeni uvicorn
# "Errno 98 address already in use" ile anında çöker → sonsuz restart döngüsü.
free_port() {
  local port="$1" tries=0 pids
  while :; do
    pids="$( { ss -ltnpH "sport = :$port" 2>/dev/null || ss -ltnp 2>/dev/null | grep -E "[:.]$port\b"; } \
             | grep -oE 'pid=[0-9]+' | cut -d= -f2 | sort -u )"
    [ -z "$pids" ] && return 0
    if [ "$tries" -ge 10 ]; then
      log "⚠ port :$port hâlâ dolu (pid: $pids) — zorla kapatılıyor"
      echo "$pids" | xargs -r kill -9 2>/dev/null || true
      sleep 1; return 0
    fi
    [ "$tries" -eq 0 ] && log "port :$port dolu (pid: $pids) — boşaltılıyor" && echo "$pids" | xargs -r kill 2>/dev/null || true
    tries=$((tries+1)); sleep 1
  done
}

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
  free_port "$BACKEND_PORT"
  log "backend başlatılıyor (uvicorn :$BACKEND_PORT)"
  (cd "$BE" && exec "$VENV/bin/uvicorn" main:app --host 127.0.0.1 --port "$BACKEND_PORT" \
      --log-level warning) >>"$LOGS/backend.log" 2>&1 &
  BACKEND_PID=$!
}

start_frontend() {
  free_port "$FRONTEND_PORT"
  log "frontend başlatılıyor (next start :$FRONTEND_PORT)"
  (cd "$FE" && VAXFORGE_API="http://127.0.0.1:$BACKEND_PORT" exec npm run start) \
      >>"$LOGS/frontend.log" 2>&1 &
  FRONTEND_PID=$!
}

cleanup() {
  log "kapatılıyor — çocuk süreçler durduruluyor"
  [ -n "$BACKEND_PID" ] && kill "$BACKEND_PID" 2>/dev/null || true
  [ -n "$FRONTEND_PID" ] && kill "$FRONTEND_PID" 2>/dev/null || true
  rm -f "$GDIR/supervisor.pid" 2>/dev/null || true
  exit 0
}
trap cleanup INT TERM

# --- başlangıç ---
# ÇİFT-INSTANCE KORUMASI: başka bir supervisor koşuyorsa çık (iki süpervizörün
# aynı portları restart etmesi de sonsuz "address in use" döngüsü doğuruyordu).
PIDFILE="$GDIR/supervisor.pid"
if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE" 2>/dev/null)" 2>/dev/null; then
  log "başka bir guardian zaten çalışıyor (pid $(cat "$PIDFILE")) — çıkılıyor"
  exit 1
fi
echo "$$" > "$PIDFILE"

log "════ VaxForge Guardian başladı ════"
start_backend
start_frontend

# servislerin ayağa kalkması için ısınma
sleep 4

declare -A FAILS=([backend]=0 [frontend]=0)
declare -A GAVE_UP=([backend]=0 [frontend]=0)   # MAX aşıldı → sus (incident spam yok)
MAX_RESTARTS="${VF_MAX_RESTARTS:-5}"

# Bir servisi denetle: ölmüşse restart et, MAX aşıldıysa BİR KEZ pes edip sus.
supervise_one() {
  local svc="$1" pid_var="$2" port="$3" starter="$4" logf="$5"
  local pid="${!pid_var}"
  # pes edilmiş servis için hiçbir şey yapma (sonsuz incident yazımını önler)
  [ "${GAVE_UP[$svc]}" -eq 1 ] && return 0

  if ! kill -0 "$pid" 2>/dev/null; then
    FAILS[$svc]=$(( FAILS[$svc] + 1 ))
    log "⚠ $svc süreci öldü (deneme ${FAILS[$svc]}/$MAX_RESTARTS)"
    incident "$svc" "process_exited" "$(tail -30 "$logf" 2>/dev/null)"
    notify "VaxForge $svc çöktü" "Yeniden başlatılıyor (:$port)"
    if [ "${FAILS[$svc]}" -le "$MAX_RESTARTS" ]; then
      "$starter"; sleep 3
    else
      GAVE_UP[$svc]=1
      log "$svc $MAX_RESTARTS kez çöktü — otomatik restart durduruldu (kod fix gerek). Bu servis artık sessiz."
      notify "VaxForge $svc PES EDİLDİ" "Kod düzeltmesi gerekiyor — site-guardian devreye alınmalı"
    fi
  else
    # süreç sağlıklı ayakta: geçici flapping sayacını sıfırla ki eski
    # hatalar birikip yanlışlıkla MAX'e ulaşmasın.
    [ "${FAILS[$svc]}" -gt 0 ] && FAILS[$svc]=0
  fi
}

while true; do
  sleep 10
  supervise_one backend  BACKEND_PID  "$BACKEND_PORT"  start_backend  "$LOGS/backend.log"
  supervise_one frontend FRONTEND_PID "$FRONTEND_PORT" start_frontend "$LOGS/frontend.log"

  # süreç canlı ama HTTP yanıt vermiyor mu? (asılı kalma)
  hc="$("$GDIR/healthcheck.sh" 2>&1)"
  echo "$hc" | grep -q "DOWN backend"  && log "⚠ backend HTTP yanıt vermiyor (süreç canlı ama asılı)"
  echo "$hc" | grep -q "DOWN frontend" && log "⚠ frontend HTTP yanıt vermiyor (süreç canlı ama asılı)"
done
