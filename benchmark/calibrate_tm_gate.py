"""TM kapısı (transmembrane.max_helices) kalibrasyonu — uzunluk ve sınıf eşli kontrol setleri.

Gerekçe: `measure_tm_gate.py` kapının KAPSAMINI ölçtü (S. suis'te proteomun %19.1'i,
sinyal-peptidi(+) havuzun %22.2'si, deneysel antijenlerin %6.2'si = 1/16, SspA).
Tek antijen "yanlış kalibrasyon" demeye yetmez; o yüzden kapı değiştirilmedi ve
FILTER_CALIBRATION kuyruğuna alındı. Bu betik kuyruktaki deneyi koşar.

FILTER_CALIBRATION.md kuralı: her sert kapı, uzunluk-eşli kontrol setlerinde
kalibre edilmeden üretime girmez. Burada ölçülen tam olarak o:

  POZİTİF set = deneysel doğrulanmış Gram(+) koruyucu antijenler, lokalizasyon
    kapısını GEÇENLER (yani TM, geri kalan tek bağlayıcı kısıt). 5 Gram(+) türün
    hepsinden toplanır — tek türün 16 antijeni bu kararı vermeye yetmez.

KARIŞIM (mixture) BİLEŞENLERİ BURADA DAHİL — ve bu, PROSPECTIVE_GROUNDTRUTH.md'deki
ölçütten bilinçli bir SAPMA; gerekçesi kayda geçiyor. `mixture=true` kuralı *recall
atfı* için vardır: koruma beş bileşenli bir karışımdan geldiğinde tek bir proteine
atfedilemez, o yüzden tek-antijen recall sayımına girmez. Buradaki soru ise atıf
değil YAPI: TM sayıcısı, sinyal peptidi(+) ve C-ucunda LPXTG/membran çapası taşıyan
Gram(+) yüzey proteinlerini iki sarmal sayıp eliyor mu? Deneysel olarak üretilmiş,
PSORTb'nin yüzey/duvar dediği bir karışım bileşeni bu yapısal sınıfın meşru üyesidir.
İlk koşumda bu bileşenler dışlanmıştı ve sonuç şuydu: pozitif setteki HİÇBİR proteinin
TM'si 1'den büyük değildi (23/23 = %100, dört kural için de aynı) — yani ölçüm, SspA
dışlandığı için kendi sorusunu göremiyordu. Kayıt: ilk koşumun çıktısı
`tm_gate_calibration.pre_mixture.json`.
Her iki tabaka da raporlanır: BİRİNCİL = karışım dahil (yapısal soru), İKİNCİL =
karışım hariç (recall sayımıyla süreklilik). Karar ölçütü birincil tabakaya uygulanır.
  NEGATİF set = gerçek çok-geçişli integral membran proteinleri (taşıyıcılar,
    permeazlar, PTS, efluks). Seçim YALNIZ anotasyonla yapılır — TM sayısıyla
    seçmek ölçümü döngüsel yapardı. Uzunluk dağılımı pozitif sete eşlenir.

Sınanan kurallar (hiçbiri varsayılmıyor, hepsi ölçülüyor):
  R1  mevcut kapı:          TM <= max_helices (1)
  R2  kapı gevşetilmiş:     TM <= 2
  R3  sinyal-duyarlı:       etkin TM = TM - 1 (sinyal peptidi(+) ise); <= 1
  R4  sinyal-duyarlı+gevşek: etkin TM = TM - 1 (sinyal peptidi(+) ise); <= 2

ÖN KAYITLI karar ölçütü (sonuçlara BAKILMADAN önce yazıldı): mevcut kural R1
yerine bir kural ancak şu iki koşulu birlikte sağlarsa önerilir —
  (a) pozitif sette duyarlılığı artırır (daha çok gerçek antijen geçer), VE
  (b) negatif sette özgüllük kaybı <= 5 yüzde puanı.
İkisini de sağlayan birden fazla kural varsa özgüllüğü en yüksek olanı seçilir.
Hiçbiri sağlamazsa karar "KAPI DEĞİŞMEZ" olur ve SspA sınıf-dışı tek vaka olarak
kayda geçer. Betik kapıyı KENDİLİĞİNDEN DEĞİŞTİRMEZ; yalnız ölçüp önerir.

Koşum: .venv/bin/python3 benchmark/calibrate_tm_gate.py
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from vaxforge import psortb, signalp, tmhmm_local  # noqa: E402
from vaxforge.config_loader import ThresholdConfig  # noqa: E402

OUT = ROOT / "benchmark" / "results" / "tm_gate_calibration.json"
SEED = 1307  # S. suis taxon — kayıt için sabit

GRAM_POSITIVE = ["ssuis", "listeria", "saureus", "mbovis", "banthracis"]

# Gerçek çok-geçişli integral membran sınıfı. Anotasyon tabanlı seçim: TM sayıcısını
# hiç kullanmaz, dolayısıyla ölçüm döngüsel değil.
NEG_KEYWORDS = [
    "permease", "abc transporter", "transporter", "symporter", "antiporter",
    "efflux", "pts system", "mfs transporter", "channel protein",
    "multidrug", "sodium:proton", "amino acid transporter", "sugar transporter",
]
# Negatif sete sızmaması gereken sınıflar: bunlar yüzeyde/salgıda olabilir ya da
# membran proteini olmayıp sadece adında taşıyıcı geçer.
NEG_EXCLUDE = [
    "atp-binding protein", "substrate-binding", "periplasmic", "lipoprotein",
    "binding protein precursor", "solute-binding",
]


def read_faa(path: Path) -> dict[str, tuple[str, str]]:
    """{accession: (description, sequence)}"""
    out: dict[str, tuple[str, str]] = {}
    acc, desc, buf = None, "", []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if acc:
                out[acc] = (desc, "".join(buf))
            head = ln[1:]
            parts = head.split(None, 1)
            acc, desc, buf = parts[0], (parts[1] if len(parts) > 1 else ""), []
        else:
            buf.append(ln.strip())
    if acc:
        out[acc] = (desc, "".join(buf))
    return out


def collect_positives(cfg_allowed: set[str]) -> list[dict]:
    """Deneysel doğrulanmış Gram(+) antijenler; lokalizasyon kapısını geçenler."""
    cands: list[dict] = []
    for bug in GRAM_POSITIVE:
        cfg_path = ROOT / "benchmark" / "bacteria" / f"{bug}.json"
        faa = ROOT / "data" / "validation" / f"{bug}_prepared.faa"
        if not (cfg_path.exists() and faa.exists()):
            print(f"  ! {bug}: config ya da proteom yok, atlandı")
            continue
        bug_cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
        seqs = read_faa(faa)
        ag: dict[str, tuple[str, bool]] = {k: (v, False)
                                           for k, v in bug_cfg.get("protective", {}).items()}
        pros = ROOT / "benchmark" / "results" / f"prospective_groundtruth_{bug}.json"
        if pros.exists():
            for a in json.loads(pros.read_text(encoding="utf-8"))["antigens"]:
                # mixture=true DAHİL: atıf sorusu değil yapı sorusu — bkz. modül
                # başlığındaki gerekçe. Tabakalı raporlanır, körlemesine katılmaz.
                if a.get("target_accession"):
                    ag.setdefault(a["target_accession"], (a["name"], bool(a.get("mixture"))))
        for acc, (name, is_mix) in sorted(ag.items()):
            if acc not in seqs:
                continue
            desc, seq = seqs[acc]
            cands.append({"bug": bug, "accession": acc, "name": name, "mixture": is_mix,
                          "length": len(seq), "seq": seq, "description": desc})
    return cands


def collect_negatives(positives: list[dict], rng: random.Random) -> list[dict]:
    """Anotasyonla seçilmiş çok-geçişli integral membran proteinleri, uzunluk eşli."""
    pool: list[dict] = []
    for bug in GRAM_POSITIVE:
        faa = ROOT / "data" / "validation" / f"{bug}_prepared.faa"
        if not faa.exists():
            continue
        pos_acc = {p["accession"] for p in positives if p["bug"] == bug}
        for acc, (desc, seq) in read_faa(faa).items():
            if acc in pos_acc:
                continue
            d = desc.lower()
            if not any(k in d for k in NEG_KEYWORDS):
                continue
            if any(k in d for k in NEG_EXCLUDE):
                continue
            pool.append({"bug": bug, "accession": acc, "name": desc.split("OS=")[0].strip(),
                         "length": len(seq), "seq": seq, "description": desc})
    # uzunluk eşleme: pozitif setin uzunluk çeyrekleri, her kovada pozitif sayısının
    # ~4 katı negatif (özgüllük tahmini için yeterli örneklem, aşırı değil)
    if not positives:
        return []
    lens = sorted(p["length"] for p in positives)
    q = [lens[0], lens[len(lens) // 4], lens[len(lens) // 2],
         lens[3 * len(lens) // 4], lens[-1]]
    picked: list[dict] = []
    for i in range(4):
        lo, hi = q[i], q[i + 1]
        if i == 3:
            hi = max(hi, q[-1])
        n_pos = sum(1 for p in positives if lo <= p["length"] <= hi)
        bucket = [x for x in pool if lo <= x["length"] <= hi
                  and x["accession"] not in {y["accession"] for y in picked}]
        rng.shuffle(bucket)
        take = bucket[: n_pos * 4]
        for t in take:
            t["length_bucket"] = f"{lo}-{hi}"
        picked.extend(take)
    return picked


def rules(tm: int, has_sp: bool, max_tm: int) -> dict[str, bool]:
    """Her kural için 'kapıyı geçer mi' (True = aday olarak kalır)."""
    eff = max(0, tm - 1) if has_sp else tm
    return {
        "R1_current_max1": tm <= max_tm,
        "R2_max2": tm <= 2,
        "R3_sp_aware_max1": eff <= 1,
        "R4_sp_aware_max2": eff <= 2,
    }


def main() -> int:
    rng = random.Random(SEED)
    cfg = ThresholdConfig.load().resolve("bacteria")
    max_tm = int(cfg["transmembrane"].params["max_helices"].value)
    sp_min = float(cfg["signal_peptide"].params["min_probability"].value)
    loc_allowed = set(cfg["localization"].params["allowed_gram_positive"].value)

    print(f"max_helices={max_tm} · signalp_min={sp_min} · izinli loc={sorted(loc_allowed)}")
    print("\n[1/4] pozitif aday antijenler toplanıyor (5 Gram(+) tür)…")
    pos_raw = collect_positives(loc_allowed)
    print(f"      {len(pos_raw)} deneysel antijen bulundu")

    # Lokalizasyon kapısı: pozitif sete yalnız yüzey/duvar sınıfı girer. TM'nin
    # geri kalan tek bağlayıcı kısıt olduğundan emin olmak için GERÇEK PSORTb.
    print("[2/4] PSORTb (Gram+) ile lokalizasyon doğrulanıyor…")
    if not psortb.available():
        print("      ! PSORTb yok — kalibrasyon yapılamaz (dürüst durdurma)")
        return 1
    loc = psortb.predict([(p["accession"], p["seq"]) for p in pos_raw], "positive")
    positives = []
    excluded_loc = []
    for p in pos_raw:
        rec = loc.get(p["accession"]) or {}
        site = str(rec.get("localization") or rec.get("site") or "").lower().replace(" ", "_")
        site = {"cytoplasmic_membrane": "membrane", "cellwall": "cell_wall"}.get(site, site)
        p["psortb"] = site or "bilinmiyor"
        p["psortb_score"] = rec.get("score")
        if site in loc_allowed:
            positives.append(p)
        else:
            excluded_loc.append({"accession": p["accession"], "name": p["name"],
                                 "bug": p["bug"], "psortb": p["psortb"]})
    print(f"      lokalizasyon kapısını geçen: {len(positives)} · "
          f"geçmeyen (sitoplazmik/moonlighting → pozitif sete girmez): {len(excluded_loc)}")
    if len(positives) < 8:
        print(f"      ! pozitif set çok küçük (n={len(positives)}) — karar verilemez")

    print("[3/4] negatif set (anotasyonla seçilmiş integral membran, uzunluk eşli)…")
    negatives = collect_negatives(positives, rng)
    print(f"      {len(negatives)} negatif seçildi")

    print("[4/4] TM + SignalP koşuluyor (gerçek araçlar)…")
    allp = positives + negatives
    pairs = [(x["accession"] + "|" + x["bug"], x["seq"]) for x in allp]
    # accession'lar türler arası çakışabilir → bileşik anahtar
    tm_real = tmhmm_local.predict(pairs)
    sp_real = signalp.predict(pairs, "bacteria")
    if not tm_real:
        print("      ! tmhmm_local sonuç vermedi — kalibrasyon yapılamaz")
        return 1
    for x in allp:
        k = x["accession"] + "|" + x["bug"]
        x["tm_helices"] = int(tm_real.get(k, 0))
        x["sp_prob"] = round(float((sp_real.get(k) or {}).get("sp_prob", 0.0)), 3)
        x["has_sp"] = x["sp_prob"] >= sp_min

    rule_names = list(rules(0, False, max_tm).keys())
    nonmix = [x for x in positives if not x.get("mixture")]

    def score(pos_subset: list[dict]) -> dict:
        out = {}
        for r in rule_names:
            tp = sum(1 for x in pos_subset if rules(x["tm_helices"], x["has_sp"], max_tm)[r])
            tn = sum(1 for x in negatives if not rules(x["tm_helices"], x["has_sp"], max_tm)[r])
            out[r] = {
                "positives_kept": tp, "n_positives": len(pos_subset),
                "sensitivity_pct": round(100.0 * tp / len(pos_subset), 1) if pos_subset else None,
                "negatives_rejected": tn, "n_negatives": len(negatives),
                "specificity_pct": round(100.0 * tn / len(negatives), 1) if negatives else None,
            }
        return out

    res = score(positives)          # BİRİNCİL: karışım dahil (yapısal soru)
    res_nonmix = score(nonmix)      # İKİNCİL: karışım hariç (recall sayımıyla süreklilik)

    # --- ön kayıtlı karar ölçütü (birincil tabakaya uygulanır) ---
    base = res["R1_current_max1"]
    qualify = []
    for r in rule_names:
        if r == "R1_current_max1":
            continue
        d_sens = res[r]["sensitivity_pct"] - base["sensitivity_pct"]
        d_spec = res[r]["specificity_pct"] - base["specificity_pct"]
        if d_sens > 0 and d_spec >= -5.0:
            qualify.append((r, d_sens, d_spec))
    qualify.sort(key=lambda t: -res[t[0]]["specificity_pct"])
    recommend = qualify[0][0] if qualify else "R1_current_max1"
    verdict = ("KAPI DEĞİŞMEZ" if recommend == "R1_current_max1"
               else f"KAPI DEĞİŞTİRİLMESİ ÖNERİLİR → {recommend}")
    # Gücün nereden geldiği denetlenebilir olsun: kapının elediği pozitifler
    dropped_by_current = [{"bug": x["bug"], "accession": x["accession"], "name": x["name"],
                           "mixture": x.get("mixture", False), "tm_helices": x["tm_helices"],
                           "sp_prob": x["sp_prob"], "psortb": x["psortb"], "length": x["length"]}
                          for x in positives if x["tm_helices"] > max_tm]

    rep = {
        "question": "TM kapısı (max_helices) Gram(+) yüzey sınıfında yanlış yerde mi duruyor?",
        "seed": SEED, "max_helices_current": max_tm, "signalp_min_probability": sp_min,
        "localization_allowed_gram_positive": sorted(loc_allowed),
        "method": {
            "positive_set": ("5 Gram(+) türün deneysel doğrulanmış koruyucu antijenleri "
                             "(+ prospektif küratörlük, mixture hariç), GERÇEK PSORTb ile "
                             "lokalizasyon kapısını geçenler"),
            "negative_set": ("aynı proteomlardan YALNIZ anotasyonla seçilmiş çok-geçişli "
                             "integral membran proteinleri (taşıyıcı/permeaz/PTS/efluks); "
                             "TM sayısıyla seçim YAPILMADI → döngüsellik yok"),
            "length_matching": "pozitif setin uzunluk çeyrekleri, kova başına 4x negatif",
            "tools": {"tm": "tmhmm_local (GERÇEK)", "signalp": "SignalP (GERÇEK)",
                      "localization": "PSORTb Gram+ (GERÇEK)"},
            "decision_criterion": ("ön kayıtlı: duyarlılık artışı > 0 VE özgüllük kaybı "
                                   "<= 5 yüzde puanı; birden fazla uyarsa en yüksek özgüllük"),
        },
        "set_sizes": {"positives": len(positives),
                      "positives_mixture": len(positives) - len(nonmix),
                      "positives_nonmixture": len(nonmix),
                      "negatives": len(negatives),
                      "positives_excluded_by_localization": len(excluded_loc)},
        "excluded_by_localization": excluded_loc,
        "rules": res,
        "rules_nonmixture_stratum": res_nonmix,
        "positives_dropped_by_current_gate": dropped_by_current,
        "qualifying_rules": [{"rule": r, "d_sensitivity_pp": round(ds, 1),
                              "d_specificity_pp": round(dp, 1)} for r, ds, dp in qualify],
        "recommendation": recommend,
        "verdict": verdict,
        "positives": [{k: v for k, v in x.items() if k != "seq"} for x in positives],
        "negatives": [{k: v for k, v in x.items() if k != "seq"} for x in negatives],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rep, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"\nBİRİNCİL tabaka (karışım dahil, n_poz={len(positives)}):")
    print(f"{'kural':<22}{'duyarlılık':<26}{'özgüllük'}")
    for r in rule_names:
        v = res[r]
        mark = "  <<< mevcut" if r == "R1_current_max1" else ""
        print(f"{r:<22}{v['positives_kept']}/{v['n_positives']} = %{v['sensitivity_pct']:<18}"
              f"{v['negatives_rejected']}/{v['n_negatives']} = %{v['specificity_pct']}{mark}")
    print(f"\nİKİNCİL tabaka (karışım hariç, n_poz={len(nonmix)}):")
    for r in rule_names:
        v = res_nonmix[r]
        print(f"{r:<22}{v['positives_kept']}/{v['n_positives']} = %{v['sensitivity_pct']:<18}"
              f"{v['negatives_rejected']}/{v['n_negatives']} = %{v['specificity_pct']}")
    if dropped_by_current:
        print(f"\nMevcut kapının elediği pozitifler ({len(dropped_by_current)}):")
        for d in dropped_by_current:
            mx = " [karışım]" if d["mixture"] else ""
            print(f"  {d['bug']:<11}{d['accession']:<14}TM={d['tm_helices']} sp={d['sp_prob']:.2f} "
                  f"{d['psortb']:<14}{d['name'][:38]}{mx}")
    else:
        print("\nMevcut kapı hiçbir pozitifi elemiyor — bu sette bedeli yok.")
    print(f"\nKARAR: {verdict}")
    print(f"yazıldı: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
