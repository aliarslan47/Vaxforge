#!/usr/bin/env bash
# VaxForge guardian — PostToolUse hook sarmalayıcısı.
# Claude bir Edit/Write yaptığında çağrılır. stdin'den hook JSON'u gelir.
# YALNIZCA düzenlenen dosya vaxforge/web/ altındaysa codecheck çalıştırır;
# değilse sessizce çıkar (başka projeleri etkilemez).

input="$(cat)"
# tool_input.file_path'i çıkar (jq varsa onunla, yoksa python)
path=""
if command -v jq >/dev/null 2>&1; then
  path="$(printf '%s' "$input" | jq -r '.tool_input.file_path // .tool_input.path // empty' 2>/dev/null)"
else
  path="$(printf '%s' "$input" | python3 -c 'import sys,json;
try:
    d=json.load(sys.stdin); ti=d.get("tool_input",{}); print(ti.get("file_path") or ti.get("path") or "")
except Exception: print("")' 2>/dev/null)"
fi

case "$path" in
  */vaxforge/web/*)
    ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    out="$("$ROOT/web/guardian/codecheck.sh" 2>&1)"; rc=$?
    if [ $rc -ne 0 ]; then
      # exit 2 => Claude'a stderr geri beslenir (guardian'ı görevlendirir)
      echo "🛡️ VaxForge guardian: kod kontrolü BAŞARISIZ — düzeltilmeli:" >&2
      echo "$out" >&2
      exit 2
    fi
    ;;
  *) : ;;  # web dışı — sessiz
esac
exit 0
