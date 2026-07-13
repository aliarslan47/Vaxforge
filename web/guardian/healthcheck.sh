#!/usr/bin/env bash
# VaxForge guardian — ÇALIŞMA-ZAMANI sağlık kontrolü.
# Backend ve frontend'in canlı olup olmadığını curl'ler. Sağlıksızsa exit≠0.
# Supervisor bunu döngüde çağırır; hangi servisin düştüğünü stdout'a yazar.

set -u
BACKEND_PORT="${VF_BACKEND_PORT:-8011}"
FRONTEND_PORT="${VF_FRONTEND_PORT:-3000}"
BACKEND_URL="http://127.0.0.1:${BACKEND_PORT}/api/health"
FRONTEND_URL="http://127.0.0.1:${FRONTEND_PORT}/"

rc=0

if curl -fs -m 5 "$BACKEND_URL" >/dev/null 2>&1; then
  echo "ok backend :$BACKEND_PORT"
else
  echo "DOWN backend :$BACKEND_PORT"
  rc=1
fi

if curl -fs -m 5 "$FRONTEND_URL" >/dev/null 2>&1; then
  echo "ok frontend :$FRONTEND_PORT"
else
  echo "DOWN frontend :$FRONTEND_PORT"
  rc=2
fi

exit $rc
