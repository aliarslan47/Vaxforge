"""Merkezi eşik yapılandırmasını yükler, doğrular ve organizma profiline çözer.

Pipeline'daki tüm eşikler config/thresholds.yaml içinde yaşar. Bu modül:
  - YAML'i yükler,
  - her parametrenin 'range' bilgisine göre değerleri doğrular,
  - seçilen organizma profili (bacteria/virus/parasite) için etkin eşikleri çözer,
  - kullanıcı override'larını (arayüzden gelen düzenlemeler) uygular.

Böylece kodun hiçbir yerinde gömülü sabit (magic number) bulunmaz.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "thresholds.yaml"


@dataclass
class ResolvedParam:
    """Belirli bir organizma için çözülmüş tek bir eşik parametresi."""

    tool: str
    name: str
    value: Any
    description: str = ""
    unit: str = ""
    range: tuple | None = None

    def in_range(self) -> bool:
        if self.range is None or not isinstance(self.value, (int, float)):
            return True
        lo, hi = self.range
        return lo <= self.value <= hi


@dataclass
class ResolvedTool:
    tool: str
    step: str
    engine: str
    description: str
    hard_filter: bool
    params: dict[str, ResolvedParam] = field(default_factory=dict)


class ThresholdConfig:
    """Yüklü config + belirli bir organizma profili için çözülmüş eşikler."""

    def __init__(self, raw: dict, path: Path):
        self.raw = raw
        self.path = path
        self.profiles: list[str] = raw["meta"]["profiles"]
        self.default_profile: str = raw["meta"]["default_profile"]

    # -- yükleme -------------------------------------------------------------
    @classmethod
    def load(cls, path: Path | str = DEFAULT_CONFIG_PATH) -> "ThresholdConfig":
        path = Path(path)
        with path.open("r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
        cfg = cls(raw, path)
        cfg.validate()
        return cfg

    # -- doğrulama -----------------------------------------------------------
    def validate(self) -> list[str]:
        """Yapısal tutarlılığı kontrol eder; sorun listesini döndürür (boşsa temiz)."""
        problems: list[str] = []
        profiles = set(self.profiles)
        for tool, spec in self.raw.get("tools", {}).items():
            for pname, pspec in spec.get("params", {}).items():
                default = pspec.get("default", {})
                if isinstance(default, dict):
                    missing = profiles - set(default)
                    if missing:
                        problems.append(f"{tool}.{pname}: eksik profil varsayılanı {sorted(missing)}")
        return problems

    # -- çözme ---------------------------------------------------------------
    def resolve(self, profile: str, overrides: dict | None = None) -> dict[str, ResolvedTool]:
        """Verilen organizma profili için tüm araçların etkin eşiklerini çözer.

        overrides: {"tool.param": deger} biçiminde arayüz düzenlemeleri.
        """
        if profile not in self.profiles:
            raise ValueError(f"Bilinmeyen profil: {profile!r}. Seçenekler: {self.profiles}")
        overrides = overrides or {}
        resolved: dict[str, ResolvedTool] = {}
        for tool, spec in self.raw.get("tools", {}).items():
            rtool = ResolvedTool(
                tool=tool,
                step=spec.get("step", ""),
                engine=spec.get("engine", ""),
                description=spec.get("description", ""),
                hard_filter=bool(spec.get("hard_filter", False)),
            )
            for pname, pspec in spec.get("params", {}).items():
                default = pspec.get("default", {})
                value = default.get(profile) if isinstance(default, dict) else default
                key = f"{tool}.{pname}"
                if key in overrides:
                    value = overrides[key]
                rng = pspec.get("range")
                rtool.params[pname] = ResolvedParam(
                    tool=tool,
                    name=pname,
                    value=value,
                    description=pspec.get("desc", ""),
                    unit=pspec.get("unit", ""),
                    range=tuple(rng) if rng else None,
                )
            resolved[tool] = rtool
        return resolved

    def candidacy_weights(self) -> dict[str, float]:
        w = dict(self.raw.get("candidacy_score", {}).get("weights", {}))
        total = sum(w.values()) or 1.0
        return {k: v / total for k, v in w.items()}


def flatten_for_report(resolved: dict[str, ResolvedTool]) -> list[dict]:
    """Çözülmüş eşikleri rapor/tablo için düz satırlara çevirir."""
    rows = []
    for rtool in resolved.values():
        for p in rtool.params.values():
            rows.append(
                {
                    "step": rtool.step,
                    "tool": rtool.tool,
                    "engine": rtool.engine,
                    "param": p.name,
                    "value": p.value,
                    "unit": p.unit,
                    "hard_filter": rtool.hard_filter,
                    "in_range": p.in_range(),
                    "description": p.description,
                }
            )
    return rows
