"""Uçtan uca pipeline smoke testi — pipeline'ın BOZULMADIĞINI garanti eder.

Her fazdan sonra koşar. Küçük demo proteom (3 protein) ile `pipeline.run()`
generator'ı tüketilir; `__result__` event'i gelmeli, meta temel anahtarları
içermeli ve run.json yazılmalı. Dış araç yoksa dürüst fallback devrede kalır —
test yine GEÇMELİ (pipeline hiçbir koşulda düşmemeli).

Hız: pipeline bir kez (module-scope fixture) koşar, tüm assert'ler onu paylaşır.
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DEMO = REPO / "data" / "samples" / "proteome_demo.faa"


@pytest.fixture(scope="module")
def pipeline_run(tmp_path_factory):
    """Demo proteom üzerinde pipeline'ı BİR kez koşar; (events, result) döndürür."""
    from vaxforge import pipeline
    from vaxforge.config_loader import ThresholdConfig
    from vaxforge.detect import detect
    from vaxforge.hosts import HostRegistry

    outdir = tmp_path_factory.mktemp("smoke_out")
    det = detect(str(DEMO))
    det.filename = DEMO.name
    cfg = ThresholdConfig.load()
    reg = HostRegistry.load()

    events, result = [], None
    for ev in pipeline.run(
        str(DEMO), det, cfg, "bacteria",
        host_names=["human"], overrides={}, has_gpu=False,
        outdir=str(outdir), host_registry=reg,
        organism_taxon=None, gram="negative", lang="tr",
    ):
        events.append(ev)
        if ev["phase"] == "__result__":
            result = ev["data"]
        if ev["phase"] == "__error__":
            pytest.fail(f"Pipeline hata verdi: {ev.get('msg')}")
    return events, result


def test_pipeline_runs_end_to_end(pipeline_run):
    events, result = pipeline_run
    assert result is not None, "__result__ event'i gelmedi"
    meta = result["meta"]
    for key in ("n_input", "n_discovery", "n_funnel", "n_epitope", "timestamp"):
        assert key in meta, f"meta['{key}'] eksik"
    assert Path(result["paths"]["json"]).exists(), "run.json yazılmadı"


def test_phase_events_are_wellformed(pipeline_run):
    events, _ = pipeline_run
    assert len(events) > 5
    for ev in events:
        assert "phase" in ev and "status" in ev and "msg" in ev
