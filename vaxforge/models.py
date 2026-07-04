"""Pipeline boyunca taşınan ortak veri yapıları."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ProteinRecord:
    """Bir aday antijen proteini ve pipeline boyunca biriken açıklamaları."""

    id: str
    seq: str
    source: str = ""                       # nasıl elde edildi (proteome/translate/orf)
    annotations: dict[str, Any] = field(default_factory=dict)
    # eleme izi: her adımın kararı ve gerekçesi
    trace: list[dict] = field(default_factory=list)
    passed: bool = True                    # sert filtrelerden geçti mi

    def mark(self, step: str, ok: bool, reason: str, **data: Any) -> None:
        self.trace.append({"step": step, "ok": ok, "reason": reason, **data})
        if not ok:
            self.passed = False


@dataclass
class Peptide:
    """Aday epitop peptidi ve tüm metrikleri."""

    seq: str
    parent: str                            # kaynak protein id
    kind: str                              # B | MHC-I | MHC-II
    start: int = 0
    metrics: dict[str, Any] = field(default_factory=dict)
    methods: dict[str, str] = field(default_factory=dict)   # metrik -> 'heuristik'/'gerçek araç'
    candidacy: float = 0.0
    passed: bool = True
    notes: list[str] = field(default_factory=list)
