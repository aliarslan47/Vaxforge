"""Kopma-dayanıklı kayıt sistemi (append-only JSONL + atomik STATUS.md).

Kullanıcı kuralı: her bitirişte kayıt al ("ne olur olmaz") + her 30 dk'da bir güncelleme.
[[egs-project]] kayıt sistemiyle aynı desen: append-only olay günlüğü asla üzerine yazılmaz,
STATUS.md atomik (tmp+rename) yazılır → süreç ortada kesilse bile son tutarlı durum kalır.

Kullanım:
    python3 benchmark/checkpoint.py "olay" "detay"        # elle kayıt
    from checkpoint import log_event, heartbeat           # koddan
"""
from __future__ import annotations

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

BENCH = Path(__file__).resolve().parent
LOG = BENCH / "PROGRESS.jsonl"
STATUS = BENCH / "STATUS.md"


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S")


def log_event(event: str, detail: str = "", phase: str = "") -> None:
    """Append-only olay günlüğü (asla üzerine yazılmaz) + STATUS.md atomik güncelle."""
    rec = {"ts": _now(), "event": event, "detail": detail, "phase": phase, "pid": os.getpid()}
    with open(LOG, "a") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    _write_status(rec)


def _write_status(last: dict) -> None:
    lines = ["# VaxForge yeni-virüs benchmark — DURUM", "",
             f"- Son güncelleme: **{last['ts']}**",
             f"- Faz: **{last.get('phase','-')}**",
             f"- Son olay: **{last['event']}**"]
    if last.get("detail"):
        lines.append(f"- Detay: {last['detail']}")
    lines += ["", "## Son 10 olay", ""]
    if LOG.exists():
        recs = LOG.read_text().splitlines()[-10:]
        for r in recs:
            try:
                d = json.loads(r)
                lines.append(f"- `{d['ts']}` [{d.get('phase','')}] {d['event']} — {d.get('detail','')}")
            except Exception:
                pass
    tmp = STATUS.with_suffix(".md.tmp")
    tmp.write_text("\n".join(lines) + "\n")
    os.replace(tmp, STATUS)  # atomik


def heartbeat(phase: str, note: str = "") -> None:
    """30-dk periyodik güncelleme için (uzun koşularda arka planda çağrılır)."""
    log_event("heartbeat", note or f"{phase} sürüyor", phase=phase)


if __name__ == "__main__":
    ev = sys.argv[1] if len(sys.argv) > 1 else "manual"
    det = sys.argv[2] if len(sys.argv) > 2 else ""
    ph = sys.argv[3] if len(sys.argv) > 3 else ""
    log_event(ev, det, ph)
    print(f"[checkpoint] {_now()} :: {ev} — {det}")
