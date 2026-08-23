"""VaxForge web backend — FastAPI.

Uçlar:
  GET  /api/health                 canlılık (supervisor + hook için)
  GET  /api/config                 profiller, konaklar, araç durumu
  GET  /api/runs                   geçmiş run özetleri (outputs/)
  GET  /api/runs/{id}              tek run.json (adaylar + meta)
  GET  /api/runs/{id}/file/{n}     çıktı dosyası (report.html/pdf, csv, xlsx, fasta)
  POST /api/run                    dosya yükle -> job başlat + SSE ilerleme akışı
  GET  /api/active-runs            devam eden koşular (UI 'çalışıyor' rozeti için)
  GET  /api/jobs/{jid}             tek job durumu (yenileme sonrası kontrol)
  GET  /api/jobs/{jid}/stream      job akışına yeniden bağlan (?cursor=N)
  POST /api/jobs/{jid}/cancel      devam eden koşuyu iptal et
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse

import jobs
from runner import (delete_run, get_config, get_run, get_run_file, list_runs)


def _sse(ev: dict) -> str:
    """Bir event dict'ini SSE 'data:' satırına çevirir."""
    return f"data: {json.dumps(ev, ensure_ascii=False)}\n\n"


_SSE_HEADERS = {
    "Cache-Control": "no-cache",
    "X-Accel-Buffering": "no",   # nginx/proxy tampon kapat
    "Connection": "keep-alive",
}

app = FastAPI(title="VaxForge API", version="1.0")

# Frontend (Next.js dev :3000 / prod) origin'lerine izin.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # yerel/tek-makine dağıtım; gerekirse daraltılır
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)

_MEDIA = {
    "report.html": "text/html",
    "report.pdf": "application/pdf",
    "candidates.csv": "text/csv",
    "candidates_full.xlsx":
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    "top_peptides.fasta": "text/plain",
    "run.json": "application/json",
}


@app.get("/api/health")
def health():
    return {"status": "ok", "service": "vaxforge-api"}


@app.get("/api/config")
def config():
    return get_config()


@app.get("/api/runs")
def runs():
    return {"runs": list_runs()}


@app.get("/api/runs/{run_id}")
def run_detail(run_id: str):
    data = get_run(run_id)
    if data is None:
        raise HTTPException(status_code=404, detail="run bulunamadı")
    return data


@app.delete("/api/runs/{run_id}")
def run_delete(run_id: str):
    """Bir koşu çıktısını kalıcı olarak siler."""
    if not delete_run(run_id):
        raise HTTPException(status_code=404, detail="run bulunamadı")
    return {"deleted": run_id}


@app.get("/api/runs/{run_id}/file/{name}")
def run_file(run_id: str, name: str):
    f = get_run_file(run_id, name)
    if f is None:
        raise HTTPException(status_code=404, detail="dosya bulunamadı")
    return FileResponse(
        str(f), media_type=_MEDIA.get(name, "application/octet-stream"),
        filename=f"{run_id}__{name}",
    )


@app.post("/api/run")
async def run(
    file: UploadFile = File(...),
    profile: str = Form("bacteria"),
    hosts: str = Form(""),          # virgülle ayrık konak adları
    gram: str = Form(""),           # negative | positive | ""
    lang: str = Form("tr"),
    adjuvant: str = Form("beta_defensin"),   # MEV adjuvan anahtarı
    strains: list[UploadFile] = File(default=[]),  # opsiyonel çok-suş (konservasyon)
):
    """Yüklenen dosyayı geçici diske yazar, koşuyu arka-plan job'ı olarak başlatır
    ve o job'ın akışını SSE ile döndürür.

    Koşu artık SSE bağlantısından bağımsız (thread'de) çalışır: tarayıcı kapansa/
    yenilense bile iş sürer; istemci `job_id` ile /api/jobs/{jid}/stream'e yeniden
    bağlanabilir. İlk event `__job__` job_id'yi taşır.
    """
    suffix = Path(file.filename or "input.dat").suffix or ".dat"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(await file.read())
        tmp_path = tmp.name

    host_names = [h.strip() for h in hosts.split(",") if h.strip()] or None
    gram_val = gram.strip() or None
    filename = file.filename or Path(tmp_path).name

    # Opsiyonel çok-suş proteomları (konservasyon). Verilmezse liste boş → dürüstçe
    # "hesaplanmadı". Her biri geçici diske yazılır (pipeline dosya-yolu bekler).
    strain_paths: list[str] = []
    for sf in strains or []:
        if not sf or not sf.filename:
            continue
        ssuf = Path(sf.filename).suffix or ".faa"
        with tempfile.NamedTemporaryFile(delete=False, suffix=ssuf) as stmp:
            stmp.write(await sf.read())
            strain_paths.append(stmp.name)

    job = jobs.create_job(filename, profile, tmp_path=tmp_path)
    jobs.start_job(job, {
        "input_path": tmp_path, "filename": filename, "profile": profile,
        "host_names": host_names, "gram": gram_val,
        "lang": lang, "adjuvant": adjuvant, "strain_paths": strain_paths,
    })

    def event_stream():
        for ev in jobs.stream_job(job, cursor=0):
            yield _sse(ev)

    return StreamingResponse(event_stream(), media_type="text/event-stream",
                             headers=_SSE_HEADERS)


@app.get("/api/active-runs")
def active_runs():
    """Devam eden koşular — UI 'çalışıyor' rozeti / yenileme sonrası kurtarma."""
    return {"runs": jobs.active_jobs()}


@app.get("/api/jobs/{job_id}")
def job_detail(job_id: str):
    st = jobs.job_status(job_id)
    if st is None:
        raise HTTPException(status_code=404, detail="job bulunamadı")
    return st


@app.get("/api/jobs/{job_id}/stream")
def job_stream(job_id: str, cursor: int = Query(0, ge=0)):
    """Devam eden (ya da yeni bitmiş) bir koşunun akışına yeniden bağlan."""
    job = jobs.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job bulunamadı")

    def event_stream():
        for ev in jobs.stream_job(job, cursor=cursor):
            yield _sse(ev)

    return StreamingResponse(event_stream(), media_type="text/event-stream",
                             headers=_SSE_HEADERS)


@app.post("/api/jobs/{job_id}/cancel")
def job_cancel(job_id: str):
    """Devam eden koşuyu iptal et (alt-süreçleri öldürür)."""
    if not jobs.cancel_job(job_id):
        raise HTTPException(status_code=404,
                            detail="job bulunamadı ya da zaten bitmiş")
    return {"cancelled": job_id}


@app.get("/")
def root():
    return JSONResponse({"service": "VaxForge API", "docs": "/docs",
                         "endpoints": ["/api/health", "/api/config",
                                       "/api/runs", "/api/run"]})
