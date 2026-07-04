"""Konak (host) kayıtları — seçilebilir konak + MHC allelleri.

config/hosts.yaml'den konakları yükler. Her konak MHC-I ve MHC-II allel
listelerini (IMGT/HLA + IPD-MHC kaynaklı) ve hangi yordayıcının uygun olduğunu
taşır. Pipeline seçilen konak(lar)ın allellerini tarar; insana kilitli değildir.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

DEFAULT_HOSTS_PATH = Path(__file__).resolve().parent.parent / "config" / "hosts.yaml"


@dataclass
class Host:
    name: str
    label: str
    source: str = ""
    mhc_i: list[str] = field(default_factory=list)
    mhc_ii: list[str] = field(default_factory=list)
    predictors: dict[str, str] = field(default_factory=dict)

    def alleles(self, mhc_class: str) -> list[str]:
        return self.mhc_i if mhc_class == "mhc_i" else self.mhc_ii

    def predictor(self, mhc_class: str) -> str:
        return self.predictors.get(mhc_class, "proxy")


class HostRegistry:
    def __init__(self, raw: dict, path: Path):
        self.raw = raw
        self.path = path
        self.hosts: dict[str, Host] = {}
        for name, spec in raw.get("hosts", {}).items():
            self.hosts[name] = Host(
                name=name, label=spec.get("label", name), source=spec.get("source", ""),
                mhc_i=list(spec.get("mhc_i", [])), mhc_ii=list(spec.get("mhc_ii", [])),
                predictors=dict(spec.get("predictors", {})),
            )
        self.default_hosts: list[str] = raw.get("default_hosts", list(self.hosts)[:1])

    @classmethod
    def load(cls, path: Path | str = DEFAULT_HOSTS_PATH) -> "HostRegistry":
        with open(path, "r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
        return cls(raw, Path(path))

    def names(self) -> list[str]:
        return list(self.hosts)

    def get(self, name: str) -> Host:
        return self.hosts[name]
