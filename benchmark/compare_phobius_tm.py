"""Phobius vs TMHMM-2.0+SignalP-5.0 — AYNI kontrol setlerinde baş-başa karşılaştırma.

Gerekçe: Faz 1'de TM kapısı sinyal-duyarlı sayıma geçirildi (`R3`), çünkü TMHMM
Gram(+) yüzey proteinlerinin N-ucundaki salgı sinyalini sarmal sayıyor. O düzeltme
ELLE yazılmış bir yamadır: "sinyal peptidi(+) ise ilk sarmalı düş". Phobius ise aynı
problemi MODELLEYEREK çözer — sinyal peptidi ve TM topolojisini tek HMM'de ortak
çözümler, yani bir bölgenin sinyal mi çapa mı olduğuna kendi karar verir
(Käll, Krogh & Sonnhammer, J Mol Biol 338:1027-1036, 2004; PMID 15111065).

Soru: Phobius, elle yazdığımız yamadan DAHA İYİ mi? Yoksa yama yeterli mi?

Kontrol setleri YENİDEN SEÇİLMEZ — Faz 1'in tam aynı setleri kullanılır
(`results/tm_gate_calibration.json` içindeki pozitif/negatif listeler). Araç
değişirken setin de değişmesi karşılaştırmayı yorumlanamaz hale getirirdi.

Sınanan kurallar:
  T1  TMHMM, kalibrasyon ÖNCESİ:   tm <= 1
  T3  TMHMM, ÜRETİMDEKİ (Faz 1):   sinyal-duyarlı etkin tm <= 1
  P1  Phobius, doğal:              phobius_tm <= 1
  P2  Phobius, gevşek:             phobius_tm <= 2
  P3  Phobius + aynı yama:         sinyal-duyarlı etkin phobius_tm <= 1
      (P3 varsayım DEĞİL, ölçüm: Phobius sinyali ayrı raporluyorsa P3 == P1 çıkar;
       çıkmazsa Phobius da sinyali sarmal sayıyor demektir.)

ÖN KAYITLI karar ölçütü (sonuçlara BAKILMADAN yazıldı): Phobius, üretimdeki T3'ün
yerine ancak şu koşulla önerilir —
  (a) duyarlılık T3'ten DÜŞÜK DEĞİL, VE
  (b) özgüllük T3'ten DÜŞÜK DEĞİL, VE
  (c) en az birinde KESİN artış var.
Eşitlik çıkarsa öneri "DEĞİŞTİRME": aynı sonucu veren yeni bir dış bağımlılık,
kurulum yükü ve atıf yükü getirir, karşılığını vermez. Bu durumda Phobius yalnız
*çapraz-doğrulama kanıtı* olarak kayda geçer (bağımsız araç aynı sonucu buldu).
Betik kapıyı KENDİLİĞİNDEN DEĞİŞTİRMEZ; yalnız ölçüp önerir.

Koşum: .venv/bin/python3 benchmark/compare_phobius_tm.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

CALIB = ROOT / "benchmark" / "results" / "tm_gate_calibration.json"
OUT = ROOT / "benchmark" / "results" / "phobius_vs_tmhmm.json"
PHOBIUS = ROOT / "tools" / "phobius-1.01" / "phobius.pl"
SP_MIN = 0.5  # Faz 1 ile aynı: signalp eşiği


def read_faa(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    acc, buf = None, []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if acc:
                out[acc] = "".join(buf)
            acc, buf = ln[1:].split(None, 1)[0], []
        else:
            buf.append(ln.strip())
    if acc:
        out[acc] = "".join(buf)
    return out


def run_phobius(seqs: dict[str, str]) -> dict[str, tuple[int, bool]]:
    """{accession: (tm_sarmal_sayısı, sinyal_peptidi_var_mı)} — phobius.pl -short."""
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.faa"
        with fa.open("w", encoding="utf-8") as fh:
            for acc, s in seqs.items():
                fh.write(f">{acc}\n")
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + "\n")
        proc = subprocess.run(
            ["perl", str(PHOBIUS), "-short", str(fa)],
            capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=1800,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"phobius rc={proc.returncode}: {proc.stderr[-800:]}")
    # -short biçimi: "ID  TM  SP  PREDICTION"; SP kolonu "Y" ya da "0"
    # (sayı değil — ilk sürümde int() denenmiş, sinyal(+) proteinlerin HEPSİ
    #  sessizce düşmüştü; ölçüm tam da sormak istediğimiz sınıfı kaybediyordu).
    res: dict[str, tuple[int, bool]] = {}
    for ln in proc.stdout.splitlines():
        parts = ln.split()
        if len(parts) < 3 or parts[1] == "ID":
            continue
        try:
            tm = int(parts[1])
        except ValueError:
            continue
        sp_tok = parts[2].strip().upper()
        if sp_tok not in ("Y", "0", "N"):
            raise RuntimeError(f"phobius SP kolonu beklenmedik: {parts[2]!r} ({ln!r})")
        res[parts[0]] = (tm, sp_tok == "Y")
    return res


def rules(tm_t: int, sp_t: bool, tm_p: int, sp_p: bool) -> dict[str, bool]:
    """True = kapıyı geçer (aday olarak kalır)."""
    eff_t = max(0, tm_t - 1) if sp_t else tm_t
    eff_p = max(0, tm_p - 1) if sp_p else tm_p
    return {
        "T1_tmhmm_max1": tm_t <= 1,
        "T3_tmhmm_sp_aware_max1": eff_t <= 1,
        "P1_phobius_max1": tm_p <= 1,
        "P2_phobius_max2": tm_p <= 2,
        "P3_phobius_sp_aware_max1": eff_p <= 1,
    }


def main() -> int:
    calib = json.loads(CALIB.read_text(encoding="utf-8"))
    if not PHOBIUS.exists():
        print(f"HATA: {PHOBIUS} yok", file=sys.stderr)
        return 2

    # --- Faz 1'in TAM AYNI setleri; dizileri proteomlardan topla ---
    proteomes: dict[str, dict[str, str]] = {}
    items: list[dict] = []
    for grp, recs in (("pos", calib["positives"]), ("neg", calib["negatives"])):
        for r in recs:
            bug = r["bug"]
            if bug not in proteomes:
                p = ROOT / "data" / "validation" / f"{bug}_prepared.faa"
                proteomes[bug] = read_faa(p) if p.exists() else {}
            seq = proteomes[bug].get(r["accession"])
            if not seq:
                print(f"  ! dizi yok, atlandı: {bug}/{r['accession']}")
                continue
            items.append({
                "group": grp, "bug": bug, "accession": r["accession"],
                "name": r.get("name") or r.get("description", "")[:60],
                "length": r["length"], "seq": seq,
                "tmhmm_tm": r["tm_helices"],
                "signalp_prob": r.get("sp_prob", 0.0),
                "signalp_pos": float(r.get("sp_prob", 0.0)) >= SP_MIN,
            })

    print(f"Phobius koşuyor: {len(items)} protein "
          f"({sum(1 for i in items if i['group'] == 'pos')} poz / "
          f"{sum(1 for i in items if i['group'] == 'neg')} neg)")
    ph = run_phobius({i["accession"]: i["seq"] for i in items})
    missing = [i["accession"] for i in items if i["accession"] not in ph]
    if missing:
        print(f"  ! Phobius çıktısı eksik: {len(missing)} protein")
    items = [i for i in items if i["accession"] in ph]
    for i in items:
        i["phobius_tm"], i["phobius_sp"] = ph[i["accession"]]
        i.pop("seq")
        i["rules"] = rules(i["tmhmm_tm"], i["signalp_pos"],
                           i["phobius_tm"], i["phobius_sp"])

    rule_names = list(rules(0, False, 0, False).keys())
    pos = [i for i in items if i["group"] == "pos"]
    neg = [i for i in items if i["group"] == "neg"]
    scores = {}
    for rn in rule_names:
        kept = sum(1 for i in pos if i["rules"][rn])
        rej = sum(1 for i in neg if not i["rules"][rn])
        scores[rn] = {
            "positives_kept": kept, "n_positives": len(pos),
            "sensitivity_pct": round(100.0 * kept / len(pos), 1),
            "negatives_rejected": rej, "n_negatives": len(neg),
            "specificity_pct": round(100.0 * rej / len(neg), 1),
        }

    # --- araç uyumu (kapıdan bağımsız, aracın kendi davranışı) ---
    sp_agree = sum(1 for i in items if i["phobius_sp"] == i["signalp_pos"])
    tm_agree = sum(1 for i in items if i["phobius_tm"] == i["tmhmm_tm"])
    agreement = {
        "signal_peptide_agreement": f"{sp_agree}/{len(items)}",
        "tm_count_exact_agreement": f"{tm_agree}/{len(items)}",
        "phobius_lower_tm": sum(1 for i in items if i["phobius_tm"] < i["tmhmm_tm"]),
        "phobius_higher_tm": sum(1 for i in items if i["phobius_tm"] > i["tmhmm_tm"]),
        # Faz 1'in can alıcı sınıfı: sinyal(+) ve TMHMM'e göre >1 sarmal
        "sp_pos_and_tmhmm_gt1": [
            {"acc": i["accession"], "name": i["name"], "group": i["group"],
             "tmhmm": i["tmhmm_tm"], "phobius": i["phobius_tm"],
             "phobius_sp": i["phobius_sp"]}
            for i in items if i["signalp_pos"] and i["tmhmm_tm"] > 1
        ],
    }

    base = scores["T3_tmhmm_sp_aware_max1"]
    qualify = []
    for rn in ("P1_phobius_max1", "P2_phobius_max2", "P3_phobius_sp_aware_max1"):
        s = scores[rn]
        d_sens = s["sensitivity_pct"] - base["sensitivity_pct"]
        d_spec = s["specificity_pct"] - base["specificity_pct"]
        if d_sens >= 0 and d_spec >= 0 and (d_sens > 0 or d_spec > 0):
            qualify.append((rn, d_sens, d_spec))
    qualify.sort(key=lambda t: (-t[2], -t[1]))
    recommend = qualify[0][0] if qualify else "T3_tmhmm_sp_aware_max1 (DEĞİŞTİRME)"

    rec = {
        "question": "Phobius (ortak SP+TM modeli), elle yazılan sinyal-duyarlı yamadan daha iyi mi?",
        "control_sets": "Faz 1 ile AYNI — results/tm_gate_calibration.json",
        "tools": {
            "phobius": "Phobius 1.01 (yerel, decodeanhmm 64-bit), PMID 15111065",
            "tmhmm": "TMHMM-2.0d (yerel) — DTU tarafından 'outdated' etiketli",
            "signalp": "SignalP-5.0b (yerel), eşik 0.5",
        },
        "decision_criterion": "ön kayıtlı: duyarlılık >= T3 VE özgüllük >= T3 VE en az birinde kesin artış; yoksa DEĞİŞTİRME",
        "set_sizes": {"positives": len(pos), "negatives": len(neg), "skipped": len(missing)},
        "rules": scores,
        "agreement": agreement,
        "qualifying_rules": [{"rule": r, "d_sensitivity_pp": a, "d_specificity_pp": b}
                             for r, a, b in qualify],
        "recommendation": recommend,
        "per_protein": items,
    }
    OUT.write_text(json.dumps(rec, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"\nkontrol seti: {len(pos)} pozitif / {len(neg)} negatif")
    print(f"{'kural':<28} {'duyarlılık':>18} {'özgüllük':>18}")
    for rn in rule_names:
        s = scores[rn]
        print(f"{rn:<28} {s['positives_kept']:>3}/{s['n_positives']:<3} "
              f"{s['sensitivity_pct']:>7.1f}%  {s['negatives_rejected']:>3}/"
              f"{s['n_negatives']:<3} {s['specificity_pct']:>7.1f}%")
    print(f"\nsinyal peptidi uyumu (Phobius vs SignalP): {agreement['signal_peptide_agreement']}")
    print(f"TM sayısı birebir uyum (Phobius vs TMHMM): {agreement['tm_count_exact_agreement']}"
          f"  (Phobius daha az: {agreement['phobius_lower_tm']}, daha çok: {agreement['phobius_higher_tm']})")
    print(f"\nsinyal(+) VE TMHMM>1 olan {len(agreement['sp_pos_and_tmhmm_gt1'])} protein:")
    for e in agreement["sp_pos_and_tmhmm_gt1"]:
        print(f"  [{e['group']}] {e['acc']} TMHMM={e['tmhmm']} → Phobius={e['phobius']} "
              f"(sp={e['phobius_sp']})  {e['name'][:44]}")
    print(f"\nÖNERİ: {recommend}")
    print(f"kayıt: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
