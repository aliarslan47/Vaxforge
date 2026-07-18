"""Aktif koşu kayıt defteri (registry) + arka-plan yürütme.

Neden: Eskiden pipeline doğrudan SSE POST generator'ının içinde çalışıyordu;
tarayıcı sayfayı yenileyince canlı görünüm kopuyor, koşuyu iptal etmenin ya da
"hâlâ çalışıyor mu" diye sormanın yolu yoktu (run_id ancak iş bitince oluşuyor).

Çözüm: Koşu, SSE bağlantısından BAĞIMSIZ bir arka-plan thread'inde çalışır.
Başlangıçta bir `job_id` verilir. Tüm event'ler job tamponunda saklanır; istemci
(yenilese bile) `job_id` ile akışa baştan/imleçten yeniden bağlanabilir, koşuyu
iptal edebilir, aktif koşuları listeleyebilir.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import queue as _queue
import signal
import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Iterator, Optional

# Bitmiş job'ları reconnect için bir süre bellekte tut (sonra temizle).
_RETAIN_SECONDS = 900  # 15 dk


@dataclass
class Job:
    id: str
    filename: str
    profile: str
    status: str = "running"          # running | done | error | cancelled
    events: list[dict] = field(default_factory=list)   # tüm event'ler (replay için)
    run_id: Optional[str] = None     # bitince pipeline'ın ürettiği çıktı dizini
    error: Optional[str] = None
    started: float = field(default_factory=time.time)
    finished: Optional[float] = None
    tmp_path: Optional[str] = None   # bitince silinecek geçici girdi dosyası
    _proc: Optional[mp.Process] = field(default=None, repr=False)  # pipeline işlemi
    _cond: threading.Condition = field(default_factory=threading.Condition, repr=False)
    _cancel: threading.Event = field(default_factory=threading.Event, repr=False)

    def add_event(self, ev: dict) -> None:
        with self._cond:
            self.events.append(ev)
            self._cond.notify_all()

    def finish(self, status: str, run_id: Optional[str] = None,
               error: Optional[str] = None) -> None:
        with self._cond:
            self.status = status
            if run_id:
                self.run_id = run_id
            if error:
                self.error = error
            self.finished = time.time()
            self._cond.notify_all()

    @property
    def cancelled(self) -> bool:
        return self._cancel.is_set()

    def request_cancel(self) -> None:
        self._cancel.set()


_JOBS: dict[str, Job] = {}
_LOCK = threading.Lock()


def _reap() -> None:
    """Süresi dolmuş bitmiş job'ları bellekten at."""
    now = time.time()
    with _LOCK:
        for jid in list(_JOBS):
            j = _JOBS[jid]
            if j.finished and (now - j.finished) > _RETAIN_SECONDS:
                del _JOBS[jid]


def create_job(filename: str, profile: str, tmp_path: str | None = None) -> Job:
    _reap()
    job = Job(id=uuid.uuid4().hex[:12], filename=filename, profile=profile,
              tmp_path=tmp_path)
    with _LOCK:
        _JOBS[job.id] = job
    return job


def get_job(job_id: str) -> Optional[Job]:
    with _LOCK:
        return _JOBS.get(job_id)


def _summary(j: Job) -> dict:
    last = j.events[-1] if j.events else {}
    return {
        "job_id": j.id,
        "filename": j.filename,
        "profile": j.profile,
        "status": j.status,
        "run_id": j.run_id,
        "phase": last.get("phase"),
        "msg": last.get("msg"),
        "n_events": len(j.events),
        "elapsed": round((j.finished or time.time()) - j.started, 1),
    }


def active_jobs() -> list[dict]:
    """UI'nın 'çalışıyor' rozetleri için: yalnız devam eden koşular."""
    _reap()
    with _LOCK:
        return [_summary(j) for j in _JOBS.values() if j.status == "running"]


def job_status(job_id: str) -> Optional[dict]:
    j = get_job(job_id)
    return _summary(j) if j else None


# ── İptal: pipeline işleminin TÜM process-grubunu öldür ─────────────────────
def _kill_process_group(pid: int) -> None:
    """`pid`'in process-grubunu (kendisi + tüm alt araçları: netMHC, bepipred…)
    topluca sonlandırır. Çocuk `setsid()` ile kendi grubunda olduğundan bu grup
    yalnız o koşuya aittir — başka projeler (ör. ConVir) asla etkilenmez."""
    try:
        pgid = os.getpgid(pid)
    except OSError:
        pgid = pid
    try:
        os.killpg(pgid, signal.SIGTERM)
    except OSError:
        pass
    time.sleep(0.4)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except OSError:
        pass


def cancel_job(job_id: str) -> bool:
    j = get_job(job_id)
    if j is None or j.status != "running":
        return False
    j.request_cancel()
    proc = j._proc
    if proc is not None and proc.pid:
        _kill_process_group(proc.pid)   # pipeline + araçlar ANINDA ölür
    return True


# ── Arka-plan yürütme (AYRI İŞLEM) ──────────────────────────────────────────
# Pipeline netMHC/bepipred gibi ağır alt-süreçler başlatır ve tek "adım" içinde
# ara event vermez. Thread olsaydı iptal ancak sıradaki event'te işlerdi (yavaş).
# Bu yüzden pipeline ayrı bir İŞLEMDE çalışır: iptal = process-grubunu SIGKILL —
# çalışan netMHC dahil her şey anında durur.
def _pipeline_child(q: "mp.Queue", kwargs: dict) -> None:
    """Çocuk işlem: kendi oturumunu açar, run_pipeline event'lerini kuyruğa yazar."""
    try:
        os.setsid()   # kendi process-grubu → topluca öldürülebilir
    except OSError:
        pass
    try:
        from runner import run_pipeline
        for ev in run_pipeline(**kwargs):
            q.put(ev)
    except Exception as exc:  # noqa: BLE001
        try:
            q.put({"phase": "__error__", "status": "error",
                   "msg": f"Sunucu hatası: {exc}", "data": None})
        except Exception:
            pass
    finally:
        try:
            q.put(None)   # sentinel: normal bitiş
        except Exception:
            pass


_MP = mp.get_context("fork")


def start_job(job: Job, kwargs: dict) -> None:
    """Pipeline'ı ayrı bir işlemde başlatır; event'leri okuyan bir daemon thread
    job tamponuna yazar, bitiş/hata/iptal'i işaretler, tmp'yi siler."""
    q: "mp.Queue" = _MP.Queue()
    proc = _MP.Process(target=_pipeline_child, args=(q, kwargs),
                       name=f"vfpipe-{job.id}", daemon=True)
    proc.start()
    job._proc = proc

    def _worker() -> None:
        try:
            while True:
                try:
                    ev = q.get(timeout=1.0)
                except _queue.Empty:
                    if not proc.is_alive():
                        break   # işlem öldü (iptal/çökme), kuyruk boş → çık
                    continue
                if ev is None:
                    break       # normal bitiş sentinel'i
                if job.cancelled:
                    break
                phase = ev.get("phase")
                job.add_event(ev)
                if phase == "__done__":
                    job.finish("done", run_id=(ev.get("data") or {}).get("run_id"))
                    break
                if phase == "__error__":
                    job.finish("error", error=ev.get("msg"))
                    break
        finally:
            if proc.is_alive():
                _kill_process_group(proc.pid)
            try:
                proc.join(timeout=3)
            except Exception:
                pass
            if job.tmp_path:
                try:
                    os.unlink(job.tmp_path)
                except OSError:
                    pass
            if job.status == "running":
                if job.cancelled:
                    job.add_event({"phase": "__cancelled__", "status": "cancelled",
                                   "msg": "Koşu iptal edildi", "data": None})
                    job.finish("cancelled")
                else:
                    # işlem done/error vermeden öldü
                    job.add_event({"phase": "__error__", "status": "error",
                                   "msg": "Koşu beklenmedik şekilde sonlandı",
                                   "data": None})
                    job.finish("error", error="beklenmeyen sonlanma")

    threading.Thread(target=_worker, name=f"vfjob-{job.id}", daemon=True).start()


def stream_job(job: Job, cursor: int = 0) -> Iterator[dict]:
    """`cursor`'dan itibaren event'leri akıtır; devam ederken yeni event bekler,
    boşta keepalive `__ping__` gönderir, koşu bitip tampon boşalınca sonlanır."""
    idx = max(0, cursor)
    # İlk bağlantı/yeniden bağlantı için job kimliğini garanti gönder.
    yield {"phase": "__job__", "status": "", "msg": "",
           "data": {"job_id": job.id, "status": job.status, "cursor": idx}}
    while True:
        with job._cond:
            while idx >= len(job.events) and job.status == "running":
                if not job._cond.wait(timeout=12):
                    break  # keepalive
            batch = job.events[idx:]
            idx += len(batch)
            status = job.status
        for ev in batch:
            yield ev
        if status != "running" and idx >= len(job.events):
            return
        if not batch and status == "running":
            yield {"phase": "__ping__", "status": "", "msg": "", "data": None}
