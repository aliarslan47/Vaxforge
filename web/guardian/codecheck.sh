#!/usr/bin/env bash
# VaxForge guardian — KOD sağlık kontrolü (hook tarafından çağrılır).
# web/** düzenlemesi sonrası: frontend TypeScript tip kontrolü + backend import kontrolü.
# Hata varsa exit≠0 ve stderr'e özet — Claude bunu görür, guardian'ı görevlendirir.

set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FE="$ROOT/web/frontend"
BE="$ROOT/web/backend"
VENV="$ROOT/.venv"
rc=0

# --- Frontend: TypeScript tip kontrolü (hızlı, derleme değil) ---
if [ -d "$FE/node_modules" ]; then
  if ! (cd "$FE" && npx --no-install tsc --noEmit >/tmp/vf_tsc.log 2>&1); then
    echo "FRONTEND TİP HATASI:" >&2
    tail -25 /tmp/vf_tsc.log >&2
    rc=1
  fi
fi

# --- Backend: Python import + söz dizimi kontrolü ---
PY="python3"
[ -x "$VENV/bin/python" ] && PY="$VENV/bin/python"
if ! (cd "$BE" && "$PY" -c "import main" >/tmp/vf_py.log 2>&1); then
  echo "BACKEND IMPORT HATASI:" >&2
  tail -25 /tmp/vf_py.log >&2
  rc=1
fi

if [ $rc -eq 0 ]; then
  echo "codecheck: temiz ✓"
fi
exit $rc
