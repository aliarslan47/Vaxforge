"""VaxForge → ImmForge glue: MEV/aday konstruktu → in silico immünizasyon → OLUR/OLMAZ.

VaxForge pipeline'ının SON adımı bunu çağırır: üretilen aşı konstruktunu (MEVConstruct.seq)
seçili konak(lar)ın MHC allelleri için GERÇEK NetMHCpan/NetMHCIIpan'den geçirip ImmForge
motorunda immünize eder → immünojenisite skoru + karar. ImmForge = adayları eleyen immün test.

ORTAM-BAĞIMSIZ: NetMHCpan/BepiPred nerede kuruluysa (yerel VEYA servis sunucusu) orada çalışır;
predict_fn enjekte edilir. immforge paketi IMMFORGE_PATH (varsayılan ~/immforge) üzerinden bulunur
(serviste pip-kurulu olur — netmhc_local'ın path deseniyle aynı).

DÜRÜSTLÜK: NetMHCpan yoksa uydurma yapmaz → {"verdict": "ARAÇ_YOK"}. Tür-özgü ImmForge
kalibrasyonu olmayan konak (ör. chicken) → human motoru proxy + açık not.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from vaxforge import bepipred, netmhc_local

IMMFORGE_PATH = os.environ.get("IMMFORGE_PATH", str(Path.home() / "immforge"))

# ImmForge'un kalibre (motor+K1+K3) türleri — bkz. immforge/config/calibration.*.yaml
_CALIBRATED = {"human", "mouse", "bovine", "pig", "goat", "sheep", "rabbit", "fish", "horse"}


def _ensure_immforge() -> None:
    if IMMFORGE_PATH not in sys.path and Path(IMMFORGE_PATH).is_dir():
        sys.path.insert(0, IMMFORGE_PATH)


def available() -> bool:
    """ImmForge skorlaması bu ortamda çalışabilir mi? (NetMHCpan + immforge paketi)."""
    if not Path(IMMFORGE_PATH).is_dir():
        return False
    return netmhc_local.available()


def _immforge_species(host_name: str) -> str:
    return host_name if host_name in _CALIBRATED else "human"


def _load_params(species: str):
    _ensure_immforge()
    from immforge.params import load_params
    cfg = Path(IMMFORGE_PATH) / "config" / f"calibration.{species}.yaml"
    return load_params(calibration_path=str(cfg)) if cfg.exists() else load_params()


def _bepitopes(seq: str) -> list[tuple[str, float]]:
    """BepiPred artık-skorları → lineer B-epitop peptitleri (varsa)."""
    if not bepipred.available():
        return []
    _ensure_immforge()
    from immforge.integrations.vaxforge import bepipred_residues_to_bepitopes
    res = bepipred.predict_residues(seq)
    return bepipred_residues_to_bepitopes(seq, res) if res else []


def score_construct(seq: str, hosts, *, name: str = "MEV") -> dict:
    """Aday dizisi + konaklar → her konak için OLUR/OLMAZ kararı + genel özet.

    hosts: vaxforge.hosts.Host listesi (.name, .label, .alleles("mhc_i"/"mhc_ii")).
    Dönüş (JSON): {tool, bepipred, headline, per_host:[...], note}.
    """
    if not seq or not hosts:
        return {"tool": "none", "verdict": "ARAÇ_YOK",
                "note": "Boş konstrukt ya da konak yok."}
    if not available():
        return {"tool": "none", "verdict": "ARAÇ_YOK",
                "note": "NetMHCpan/ImmForge bu ortamda yok — serviste kurulu olacak."}

    _ensure_immforge()
    from immforge.integrations.vaxforge import evaluate_candidate

    b_eps = _bepitopes(seq)
    per_host: list[dict] = []
    for h in hosts:
        species = _immforge_species(h.name)
        params = _load_params(species)
        v = evaluate_candidate(
            name, seq,
            mhci_alleles=list(h.alleles("mhc_i")),
            mhcii_alleles=list(h.alleles("mhc_ii")),
            predict_fn=netmhc_local.predict,
            species=species, params=params,
            b_epitopes=b_eps or None,
        )
        d = v.to_dict()
        d["host"] = h.name
        d["host_label"] = getattr(h, "label", h.name)
        d["species_calibrated"] = h.name in _CALIBRATED
        if not d["species_calibrated"]:
            d["reasons"] = list(d.get("reasons", [])) + [
                f"'{h.name}' için ImmForge tür-özgü kalibrasyonu yok → human motoru proxy (dürüst)."]
        per_host.append(d)

    headline = per_host[0]
    return {
        "tool": "netmhcpan",
        "bepipred": bool(b_eps),
        "n_bepitopes": len(b_eps),
        "headline": {"verdict": headline["verdict"], "score": headline["score"],
                     "host": headline["host"]},
        "per_host": per_host,
        "note": "In silico immünojenisite ÖNGÖRÜSÜ (garanti değil); eşikler operasyonel "
                "(OLUR≥60/ZAYIF≥35). Motor kinetiği peptid-aşısı held-out'unda doğrulandı.",
    }
