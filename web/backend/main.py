"""VaxForge web backend — FastAPI.

Uçlar:
  GET  /api/health              canlılık (supervisor + hook için)
  GET  /api/config              profiller, konaklar, araç durumu
  GET  /api/runs                geçmiş run özetleri (outputs/)
  GET  /api/runs/{id}           tek run.json (adaylar + meta)
  GET  /api/runs/{id}/file/{n}  çıktı dosyası (report.html/pdf, csv, xlsx, fasta)
  POST /api/run                 dosya yükle -> SSE ilerleme akışı
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse

from runner import (get_config, get_run, get_run_file, list_runs,
                    run_pipeline)

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
):
    """Yüklenen dosyayı geçici diske yazar ve pipeline'ı SSE ile akıtır."""
    suffix = Path(file.filename or "input.dat").suffix or ".dat"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(await file.read())
        tmp_path = tmp.name

    host_names = [h.strip() for h in hosts.split(",") if h.strip()] or None
    gram_val = gram.strip() or None
    filename = file.filename or Path(tmp_path).name

    def event_stream():
        try:
            for ev in run_pipeline(tmp_path, filename, profile,
                                   host_names, gram_val, lang):
                yield f"data: {json.dumps(ev, ensure_ascii=False)}\n\n"
        except Exception as exc:  # pipeline içi beklenmeyen çökme
            err = {"phase": "__error__", "status": "error",
                   "msg": f"Sunucu hatası: {exc}", "data": None}
            yield f"data: {json.dumps(err, ensure_ascii=False)}\n\n"
        finally:
            try:
                Path(tmp_path).unlink(missing_ok=True)
            except Exception:
                pass

    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",   # nginx/proxy tampon kapat
            "Connection": "keep-alive",
        },
    )


@app.get("/")
def root():
    return JSONResponse({"service": "VaxForge API", "docs": "/docs",
                         "endpoints": ["/api/health", "/api/config",
                                       "/api/runs", "/api/run"]})
