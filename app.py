"""VaxForge — Streamlit arayüzü (prototip, hafif omurga).

Akış: dosya yükle -> otomatik tanı -> planlanan adımları göster ->
organizma profiline göre düzenlenebilir eşikleri göster.

Ağır adımlar (AlphaFold peptit-MHC, docking/MD) bu makinede GPU olmadığı için
'deferred' rozetiyle işaretlenir; ilerideki GPU makinesine/buluta gider.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from vaxforge import mhc_real
from vaxforge.config_loader import ThresholdConfig
from vaxforge.detect import detect
from vaxforge.hosts import HostRegistry
from vaxforge.plan import build_plan, plan_table

st.set_page_config(page_title="VaxForge", page_icon="🧬", layout="wide")

CFG = ThresholdConfig.load()
HOSTS = HostRegistry.load()

# --- Başlık -----------------------------------------------------------------
st.title("🧬 VaxForge")
st.caption(
    "Ajan destekli in silico reverse vaccinology pipeline — prototip. "
    "Dosyayı yükleyin; sistem tipini tanıyıp ne yapacağını planlar."
)

# --- Kenar çubuğu: organizma profili + ortam --------------------------------
with st.sidebar:
    st.header("Ayarlar")
    profile = st.selectbox(
        "Patojen profili (etken)",
        CFG.profiles,
        index=CFG.profiles.index(CFG.default_profile),
        help="Eşik varsayılanları patojene göre değişir.",
    )
    host_choices = HOSTS.names()
    selected_hosts = st.multiselect(
        "Konak(lar) — MHC taranacak",
        host_choices,
        default=HOSTS.default_hosts,
        format_func=lambda n: HOSTS.get(n).label,
        help="Seçilen her konağın MHC-I ve MHC-II allelleri taranır. "
             "İnsana kilitli değil; birden çok konak seçilebilir (tür-kapsamı).",
    )
    for n in selected_hosts:
        h = HOSTS.get(n)
        pi = "GERÇEK" if (h.predictor("mhc_i") == "mhcflurry" and mhc_real.available()) else "proxy"
        st.caption(f"• {h.label}: MHC-I {len(h.mhc_i)} allel ({pi}), "
                   f"MHC-II {len(h.mhc_ii)} allel (proxy)")
    has_gpu = st.toggle(
        "Yerel GPU var", value=False,
        help="Kapalıysa AlphaFold/MD adımları 'deferred' (uzak worker) işaretlenir.",
    )
    st.divider()
    st.caption(f"config: `{CFG.path.name}` · hosts: `{HOSTS.path.name}` · "
               f"MHCflurry: {'var' if mhc_real.available() else 'yok'}")

# --- Girdi: dosya yükleme veya örnek ----------------------------------------
st.subheader("📁 Girdi dosyası")
_SAMPLES = {
    "🧪 Patojen proteinleri (gerçek VFDB)": "data/samples/pathogen_demo.faa",
    "🧬 Proteom (protein FASTA)": "data/samples/proteome_demo.faa",
    "🔤 Genler / CDS (nükleotid)": "data/samples/genes_demo.fna",
    "📄 Ham okumalar (FASTQ)": "data/samples/reads_demo.fastq",
}
_up_col, _s_col = st.columns([2, 1])
with _up_col:
    uploaded = st.file_uploader(
        "Dosyanızı sürükleyip bırakın — FASTA / FASTQ / GenBank (.gz destekli)",
        type=["fasta", "fa", "faa", "fna", "fastq", "fq", "gz", "gb", "gbk", "genbank", "gbff"],
        help="Protein/nükleotid FASTA, ham okuma FASTQ, veya anotasyonlu GenBank (.gb/.gbk). "
             "GenBank'ta CDS'ler çeviri + gen + lokus ile doğrudan alınır. Tip otomatik tanınır.",
    )
with _s_col:
    st.caption("ya da örnek dosyayla dene:")
    for _label, _path in _SAMPLES.items():
        if st.button(_label, key=f"smp_{_path}", use_container_width=True):
            st.session_state["sample"] = _path

# Girdiyi çöz: yüklenen dosya öncelikli, yoksa seçilen örnek
input_path = input_name = None
if uploaded is not None:
    suffix = Path(uploaded.name).suffix or ".dat"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(uploaded.getbuffer())
        input_path = tmp.name
    input_name = uploaded.name
    st.session_state.pop("sample", None)
elif st.session_state.get("sample") and Path(st.session_state["sample"]).exists():
    input_path = st.session_state["sample"]
    input_name = Path(input_path).name

if input_path is None:
    st.info("👆 Bir dosya yükleyin ya da sağdaki örneklerden birini seçin. "
            "Desteklenen: protein/nükleotid FASTA veya FASTQ (.gz olabilir).")
    st.stop()

st.success(f"Girdi: **{input_name}**")
det = detect(input_path)
det.filename = input_name

# --- Tanıma sonucu ----------------------------------------------------------
st.subheader("1) Dosya tanıma")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Format", det.fmt.upper())
c2.metric("Tip", det.seq_type)
c3.metric("Molekül", det.molecule)
c4.metric("Kayıt", f"{det.num_seqs}")
c1.metric("Ort. uzunluk", f"{det.avg_len:.0f}")
c2.metric("En kısa", f"{det.min_len}")
c3.metric("En uzun", f"{det.max_len}")
c4.metric("Sıkıştırma", "gzip" if det.is_gzipped else "yok")

if not det.confident:
    st.warning("Tanıma düşük güvenle yapıldı — sonuçları elle doğrulayın.")
for note in det.notes:
    st.write("•", note)

# --- Planlanan adımlar ------------------------------------------------------
st.subheader("2) Planlanan pipeline")
steps = build_plan(det, has_gpu=has_gpu)
df_plan = pd.DataFrame(plan_table(steps))


def _style_status(val: str) -> str:
    return {
        "deferred": "color:#b26a00;font-weight:600",
        "done": "color:#2e7d32",
        "running": "color:#1565c0",
    }.get(val, "")


st.dataframe(
    df_plan.style.map(_style_status, subset=["durum"]),
    use_container_width=True, hide_index=True,
)
n_def = sum(1 for s in steps if s.status == "deferred")
if n_def:
    st.info(f"{n_def} ağır adım (GPU) bu makinede ertelendi; ayrı GPU makinesine/buluta gidecek. "
            "Yan panelden 'Yerel GPU var' açılırsa yerelde koşarlar.")

# --- Eşikler (organizma profiline göre, düzenlenebilir) ---------------------
st.subheader(f"3) Eşikler — profil: `{profile}`")
st.caption("Tüm eşikler config'ten gelir; organizmaya göre değişir ve burada düzenlenebilir. "
           "Koşulan değerler rapora yazılır (tekrarlanabilirlik).")

resolved = CFG.resolve(profile)
by_step: dict[str, list] = {}
for rtool in resolved.values():
    by_step.setdefault(rtool.step, []).append(rtool)

step_titles = {s.id: s.title for s in steps}
for step_id, tools in by_step.items():
    with st.expander(f"▸ {step_titles.get(step_id, step_id)}", expanded=False):
        for rtool in tools:
            badge = "🔒 sert filtre" if rtool.hard_filter else "◦ skorlama"
            st.markdown(f"**{rtool.tool}** — {badge}  \n<small>{rtool.engine} · {rtool.description}</small>",
                        unsafe_allow_html=True)
            for p in rtool.params.values():
                label = f"{p.name}" + (f" ({p.unit})" if p.unit else "")
                key = f"{rtool.tool}.{p.name}"
                if isinstance(p.value, bool):
                    st.checkbox(label, value=p.value, key=key, help=p.description)
                elif isinstance(p.value, (int, float)) and p.range:
                    st.number_input(label, value=float(p.value),
                                    min_value=float(p.range[0]), max_value=float(p.range[1]),
                                    key=key, help=p.description)
                else:
                    st.text_input(label, value=str(p.value), key=key, help=p.description)

# --- Çalıştır ---------------------------------------------------------------
st.subheader("4) Çalıştır")
_tool_status = []
for _mod, _lbl in [("discovery", "DIAMOND+VFDB"), ("deeploc", "DeepLoc"),
                   ("tmhmm_local", "TMHMM"), ("signalp", "SignalP"),
                   ("iapred", "IApred"), ("bepipred", "BepiPred"),
                   ("netmhc_local", "NetMHCpan"), ("toxinpred", "ToxinPred2")]:
    try:
        _m = __import__(f"vaxforge.{_mod}", fromlist=["available"])
        _tool_status.append(f"{'✅' if _m.available() else '⚠️'} {_lbl}")
    except Exception:
        _tool_status.append(f"⚠️ {_lbl}")
st.caption("Gerçek araç durumu: " + " · ".join(_tool_status)
           + "  \n(⚠️ = araç yok, o adım için dürüst-etiketli yedek yöntem kullanılır. "
             "Gerçek bir koşu birkaç dakika sürer — DeepLoc/NetMHCpan CPU'da yavaştır.)")

if st.button("▶ Pipeline'ı başlat", type="primary"):
    # arayüzden düzenlenen sayısal/bool eşikleri override olarak topla
    overrides = {}
    for rtool in resolved.values():
        for p in rtool.params.values():
            key = f"{rtool.tool}.{p.name}"
            if key in st.session_state and isinstance(p.value, (int, float, bool)):
                overrides[key] = st.session_state[key]

    from vaxforge import pipeline

    progress = st.progress(0.0, text="Başlıyor…")
    log = st.container()
    total = len(steps)
    done = 0
    result = None
    with st.status("Pipeline çalışıyor…", expanded=True) as status:
        for ev in pipeline.run(input_path, det, CFG, profile,
                               host_names=selected_hosts, overrides=overrides,
                               has_gpu=has_gpu, outdir="outputs", host_registry=HOSTS):
            ph, stt, msg = ev["phase"], ev["status"], ev["msg"]
            if ph == "__result__":
                result = ev["data"]
                continue
            if ph == "__error__":
                status.update(label="Hata", state="error")
                st.error(msg)
                st.stop()
            icon = {"running": "⏳", "done": "✅", "deferred": "⏸️"}.get(stt, "•")
            log.write(f"{icon} **{ph}** — {msg}")
            if stt in ("done", "deferred"):
                done += 1
                progress.progress(min(1.0, done / total), text=msg)
        status.update(label="Tamamlandı ✅", state="complete")

    if result:
        peptides = result["peptides"]
        construct = result["construct"]
        paths = result["paths"]
        rmeta = result["meta"]

        # --- Eleme akışı: kaç girdi → her adımda kaç kaldı ------------------
        st.subheader("📊 Eleme akışı (kaç girdi → kaç kaldı)")
        mol = rmeta.get("molecule", "dizi")
        unit = "CDS" if mol == "cds" else ("okuma" if mol == "reads" else "protein")
        n_pep_pass = sum(1 for p in peptides if p.passed)
        flow = [
            (f"1· Girdi ({unit})", rmeta.get("n_raw", "?"), "dosyadaki ham dizi sayısı"),
            ("2· Proteine çevrildi", rmeta.get("n_input", "?"), "CDS→çeviri / ORF; ≥20 aa"),
            ("3· Keşif sonrası (VFDB)", rmeta.get("n_discovery", "?"), "DIAMOND ile virülans DB eşleşmesi"),
            ("4· Huni sonrası (protein)", rmeta.get("n_funnel", "?"), "DeepLoc+TMHMM+SignalP+IApred+insan homoloji"),
            ("5· Üretilen epitop", rmeta.get("n_epitope", "?"), "SLIDING-WINDOW: B + MHC-I + MHC-II"),
            ("6· Alerjenite sonrası", rmeta.get("n_after_allergen", "?"), "FAO/WHO 6-mer — alerjenler elendi"),
            ("7· Toksisite sonrası", rmeta.get("n_after_toxicity", "?"), "ToxinPred2 — toksikler elendi"),
            ("8· Sıralanan aday", len(peptides), "adaylık puanına göre sıralandı"),
        ]
        fdf = pd.DataFrame([{"Adım": s, "Kalan": n, "Açıklama": d} for s, n, d in flow])
        cflow, cbar = st.columns([1, 1])
        cflow.dataframe(fdf, use_container_width=True, hide_index=True)
        # huni grafiği
        try:
            import matplotlib.pyplot as _plt
            steps_lbl = [s for s, n, d in flow if isinstance(n, int)]
            steps_val = [n for s, n, d in flow if isinstance(n, int)]
            fig, ax = _plt.subplots(figsize=(4.6, 3))
            ax.barh(range(len(steps_val)), steps_val, color="#0b6b4f")
            ax.set_yticks(range(len(steps_val))); ax.set_yticklabels(steps_lbl, fontsize=8)
            ax.invert_yaxis()
            for i, v in enumerate(steps_val):
                ax.text(v, i, f" {v}", va="center", fontsize=8)
            ax.set_xlabel("kalan sayı"); fig.tight_layout()
            cbar.pyplot(fig)
        except Exception:
            pass
        st.caption("Not: epitop adımında sliding-window ile her proteinden çok sayıda "
                   "peptit üretilir; sonraki adımlar bunları eler.")

        st.subheader("Sonuç — en iyi aday peptitler")
        rows = [{"sıra": i + 1, "peptit": p.seq, "tip": p.kind, "skor": p.candidacy,
                 "CDS / kaynak protein": p.parent,
                 "gen": p.metrics.get("gene") or "—",
                 "lokus_tag": p.metrics.get("locus_tag") or "—",
                 "konum": p.metrics.get("location") or "—",
                 "konak sunumu": ", ".join(p.metrics.get("hosts_presented", []) or []) or "—",
                 "en iyi allel": p.metrics.get("best_allele", "—"),
                 "alerjen": p.metrics.get("allergen"), "toksik": p.metrics.get("toxicity")}
                for i, p in enumerate(peptides[:20])]
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        # --- Her aday için yapılan TÜM analizler --------------------------
        st.subheader("🔬 Her aday için yapılan analizler")
        st.caption("Her peptidin geçtiği analizler ve bunları üreten gerçek araç.")
        for i, p in enumerate(peptides[:10], 1):
            with st.expander(f"{i}. {p.seq} — {p.kind} — adaylık {p.candidacy}"):
                a = [("Antijenite — kaynak protein (IApred)", p.metrics.get("parent_antigenicity"))]
                if p.kind in ("MHC-I", "MHC-II"):
                    tool = "NetMHCpan" if p.kind == "MHC-I" else "NetMHCIIpan"
                    a += [(f"{p.kind} bağlanma %rank ({tool})", p.metrics.get("pseudo_rank")),
                          ("En iyi allel", p.metrics.get("best_allele")),
                          ("Sunan konak(lar)", ", ".join(p.metrics.get("hosts_presented", []) or []) or "—"),
                          ("Konak/allel kapsamı", p.metrics.get("host_coverage"))]
                if p.kind == "B":
                    a += [("B-hücre skoru (BepiPred-1.0)", p.metrics.get("bepipred")),
                          ("Kolaskar-Tongaonkar antijenite", p.metrics.get("kolaskar")),
                          ("Parker hidrofilisite", p.metrics.get("parker"))]
                a += [("Alerjenite (FAO/WHO 6-mer)", "ALERJEN" if p.metrics.get("allergen") else "temiz"),
                      ("Toksisite skoru (ToxinPred2)", p.metrics.get("toxicity")),
                      ("Kaynak CDS / protein", p.parent)]
                if p.metrics.get("gene"):
                    a += [("Gen / lokus / konum",
                           f"{p.metrics.get('gene')} / {p.metrics.get('locus_tag')} / {p.metrics.get('location')}")]
                st.table(pd.DataFrame(a, columns=["Analiz (araç)", "Sonuç"]))

        st.subheader("mRNA konstrüktü")
        cm = construct.metrics
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("mRNA uzunluk", f"{cm['mrna_len']} nt")
        m2.metric("GC %", cm["gc_percent"])
        m3.metric("CAI (insan)", cm["cai_human"])
        m4.metric("CTL/HTL/B", f"{cm['n_ctl']}/{cm['n_htl']}/{cm['n_bcell']}")
        st.code(construct.mrna, language="text")

        st.subheader("İndirilebilir çıktılar")
        cols = st.columns(len(paths))
        labels = {"csv": "Adaylar (CSV)", "fasta": "Peptitler (FASTA)",
                  "genbank": "mRNA (GenBank)", "json": "Tam koşum (JSON)",
                  "html": "Rapor (HTML)", "pdf": "Rapor (PDF)"}
        for col, (kind, path) in zip(cols, paths.items()):
            data = Path(path).read_bytes()
            col.download_button(labels.get(kind, kind), data=data,
                                file_name=Path(path).name, key=f"dl_{kind}")

        # Kaynaklar / atıflar (kullanılan tüm araçlar)
        refs = result["meta"].get("citations", [])
        if refs:
            with st.expander(f"📚 Kaynaklar / atıflar ({len(refs)} araç)"):
                for i, r in enumerate(refs, 1):
                    doi = r["doi"]
                    link = doi if doi.startswith("http") else f"https://doi.org/{doi}"
                    st.markdown(f"**[{i}] {r['tool']}** — _{r['step']}_  \n{r['citation']} [{link}]({link})")
