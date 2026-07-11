"""PSORTb — prokaryota-özel hücre-altı lokalizasyon (podman container).

DeepLoc ökaryot-yanlıdır; bakteri proteinlerini yanlış konumlandırabilir
(ör. MenB'de yüzey lipoproteini fHbp'yi 'sitoplazma' sandı). PSORTb bakteriye
özel eğitilmiştir ve Gram-tipine göre doğru bölmeleri (sitoplazma / iç membran /
periplazma / dış membran / hücre duvarı / hücre-dışı) tanır. RV'nin özü yüzey/
salgı proteinlerini bulmak olduğundan bakteri dalında PSORTb biyolojik doğrudur.

Kurulum: podman container (brinkmanlab/psortb_commandline). Kurulu değilse
available()=False → funnel DeepLoc'a düşer (dürüst 'method' etiketi).

Gram bilgisi ZORUNLU: --positive (Gram+) veya --negative (Gram−). Yanlış Gram =
yanlış tahmin, bu yüzden kullanıcı arayüzden seçer.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

IMAGE = "docker.io/brinkmanlab/psortb_commandline:1.0.2"

# PSORTb 'Final_Localization' -> config bakteri lokalizasyon sözlüğü.
# config allowed(bacteria) = [extracellular, outer_membrane, cell_wall, periplasm]
_LOC_MAP = {
    "extracellular": "extracellular",
    "outermembrane": "outer_membrane",
    "cellwall": "cell_wall",
    "periplasmic": "periplasm",
    "cytoplasmicmembrane": "membrane",   # iç membran — bakteri allowed'da yok → elenir
    "cytoplasmic": "cytoplasm",          # elenir
    "unknown": "unknown",                # PSORTb belirsiz → funnel gate'i 'atlar' (elemez)
}


def _podman() -> str | None:
    return shutil.which("podman")


def available() -> bool:
    """podman var VE image yerelde mevcut mu?"""
    pod = _podman()
    if not pod:
        return False
    try:
        r = subprocess.run([pod, "image", "exists", IMAGE],
                           capture_output=True, timeout=30)
        return r.returncode == 0
    except Exception:
        return False


def _parse_terse(text: str) -> dict[str, dict]:
    """PSORTb -o terse çıktısı: başlık satırı + SeqID<TAB>Localization<TAB>Score.

    SeqID PSORTb tarafından ilk boşluğa kadar kısaltılır; FASTA başlığını
    boşluksuz (s0, s1…) verdiğimiz için birebir eşleşir.
    """
    out: dict[str, dict] = {}
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return out
    # başlık satırını atla (SeqID ile başlar)
    for ln in lines:
        if ln.lower().startswith("seqid"):
            continue
        parts = ln.split("\t")
        if len(parts) < 2:
            continue
        sid = parts[0].strip()
        raw_loc = parts[1].strip()
        try:
            score = float(parts[2]) if len(parts) > 2 and parts[2].strip() else 0.0
        except ValueError:
            score = 0.0
        key = raw_loc.replace(" ", "").replace("/", "").lower()
        mapped = _LOC_MAP.get(key, "unknown")
        out[sid] = {"localization": mapped, "raw": raw_loc, "score": round(score, 2)}
    return out


def predict(seq_pairs: list[tuple[str, str]], gram: str) -> dict[str, dict]:
    """[(id, seq)] -> {id: {localization, raw, score}}. Hata/eksik -> {}.

    gram: 'positive' | 'negative'. Sentetik s0/s1… ID kullanılır (boşluksuz),
    gerçek ID'ye çağıran tarafta eşlenir.
    """
    pod = _podman()
    if not pod or gram not in ("positive", "negative") or not seq_pairs:
        return {}
    flag = "--positive" if gram == "positive" else "--negative"
    workdir = Path(tempfile.mkdtemp(prefix="psortb_"))
    infile = workdir / "input.faa"
    # sentetik boşluksuz ID -> gerçek ID
    id_map = {f"s{i}": pid for i, (pid, _) in enumerate(seq_pairs)}
    infile.write_text(
        "".join(f">s{i}\n{seq}\n" for i, (_, seq) in enumerate(seq_pairs)),
        encoding="utf-8")
    # PSORTb 'psort' wrapper quirk'leri:
    #  - '-r' seçeneğini KABUL ETMEZ; sonucu container-içi /tmp/results/ altına yazar.
    #  - Bayesian model yolu CWD'ye GÖRELİ → mutlaka /usr/local/src'de koşmalı.
    # Bu yüzden entrypoint'i bash'e çeviriyoruz: doğru WD'de koştur, sonucu cat et.
    script = (f"cd /usr/local/src && /usr/local/psortb/bin/psort "
              f"-i /data/input.faa {flag} -o terse >/dev/null 2>&1; "
              f"cat /tmp/results/*_psortb_*.txt 2>/dev/null")
    try:
        cmd = [pod, "run", "--rm", "-v", f"{workdir}:/data:Z",
               "--entrypoint", "/bin/bash", IMAGE, "-c", script]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except Exception:
        shutil.rmtree(workdir, ignore_errors=True)
        return {}
    shutil.rmtree(workdir, ignore_errors=True)
    parsed = _parse_terse(r.stdout)
    # sentetik ID'leri gerçek ID'ye çevir
    return {id_map[sid]: d for sid, d in parsed.items() if sid in id_map}
