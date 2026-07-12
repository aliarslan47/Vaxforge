"""VaxForge — Streamlit arayüzü (prototip, hafif omurga).

Akış: dosya yükle -> otomatik tanı -> planlanan adımları göster ->
organizma profiline göre düzenlenebilir eşikleri göster.

Yapısal doğrulama (AlphaFold peptit-MHC) + docking/MD adımları şimdilik
çıkarıldı; odak aday-belirleme. GPU gelince geri eklenebilir.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from vaxforge import netmhc_local
from vaxforge.config_loader import ThresholdConfig
from vaxforge.detect import detect
from vaxforge.hosts import HostRegistry
from vaxforge.plan import build_plan, plan_table

st.set_page_config(page_title="VaxForge — Reverse Vaccinology", page_icon="🧬",
                   layout="wide", initial_sidebar_state="expanded")

CFG = ThresholdConfig.load()
HOSTS = HostRegistry.load()

# --- Kurumsal tema (CSS) ----------------------------------------------------
st.markdown("""
<style>
:root{
  --vf-bg:#0b1220; --vf-panel:#111a2e; --vf-line:#1e2a44;
  --vf-ink:#e6ecf7; --vf-mut:#93a2c0;
  --vf-brand:#2dd4bf; --vf-brand2:#6366f1; --vf-accent:#38bdf8;
}
html,body,[class*="css"]{ font-feature-settings:"cv02","cv03","cv04"; }
.stApp{ background:radial-gradient(1200px 600px at 12% -8%, #14213d 0%, var(--vf-bg) 55%); }
#MainMenu, footer, [data-testid="stToolbar"]{ visibility:hidden; height:0; }
.block-container{ padding-top:1.4rem; padding-bottom:3rem; max-width:1280px; }
h1,h2,h3{ letter-spacing:-.01em; }
[data-testid="stSidebar"]{ background:linear-gradient(180deg,#0d1526,#0b1220);
  border-right:1px solid var(--vf-line); }
[data-testid="stSidebar"] .stMarkdown, [data-testid="stSidebar"] label{ color:var(--vf-ink); }

/* Hero */
.vf-hero{ display:flex; align-items:center; gap:18px; padding:22px 26px; margin-bottom:6px;
  background:linear-gradient(110deg, rgba(45,212,191,.14), rgba(99,102,241,.14) 60%, rgba(56,189,248,.10));
  border:1px solid var(--vf-line); border-radius:18px;
  box-shadow:0 10px 40px -20px rgba(56,189,248,.5); }
.vf-hero .vf-logo{ font-size:44px; line-height:1;
  filter:drop-shadow(0 4px 14px rgba(45,212,191,.5)); }
.vf-hero h1{ margin:0; font-size:30px; font-weight:800; color:var(--vf-ink); }
.vf-hero .vf-sub{ margin:4px 0 0; color:var(--vf-mut); font-size:14px; max-width:720px; }
.vf-chip{ display:inline-block; padding:3px 10px; border-radius:999px; font-size:11px;
  font-weight:700; letter-spacing:.03em; text-transform:uppercase;
  background:rgba(45,212,191,.16); color:var(--vf-brand); border:1px solid rgba(45,212,191,.35); }
.vf-chips{ display:flex; gap:8px; flex-wrap:wrap; margin:14px 0 26px; }
.vf-chips .c{ padding:6px 12px; border-radius:10px; font-size:12px; color:var(--vf-mut);
  background:var(--vf-panel); border:1px solid var(--vf-line); }
.vf-chips .c b{ color:var(--vf-ink); font-weight:700; }

/* Bileşenler */
.stButton>button{ border-radius:11px; border:1px solid var(--vf-line);
  background:var(--vf-panel); color:var(--vf-ink); font-weight:600; transition:.15s; }
.stButton>button:hover{ border-color:var(--vf-brand); color:#fff;
  box-shadow:0 6px 20px -10px rgba(45,212,191,.7); transform:translateY(-1px); }
.stButton>button[kind="primary"]{ border:0;
  background:linear-gradient(90deg,var(--vf-brand),var(--vf-brand2)); color:#04121a; font-weight:800; }
[data-testid="stMetric"]{ background:var(--vf-panel); border:1px solid var(--vf-line);
  border-radius:14px; padding:12px 16px; }
[data-testid="stMetricValue"]{ color:var(--vf-ink); font-weight:800; }
[data-testid="stExpander"]{ border:1px solid var(--vf-line); border-radius:12px;
  background:rgba(17,26,46,.5); }
[data-testid="stExpander"] summary:hover{ color:var(--vf-brand); }
hr{ border-color:var(--vf-line); }
.vf-sec{ display:flex; align-items:center; gap:10px; margin:26px 0 4px; }
.vf-sec .n{ width:26px; height:26px; border-radius:8px; display:grid; place-items:center;
  font-size:13px; font-weight:800; color:#04121a;
  background:linear-gradient(135deg,var(--vf-brand),var(--vf-accent)); }
.vf-sec h3{ margin:0; font-size:19px; color:var(--vf-ink); }
</style>
""", unsafe_allow_html=True)

# --- Hero başlık ------------------------------------------------------------
_n_tools = len(CFG.resolve(CFG.default_profile))
_n_hosts = len(HOSTS.names())
st.markdown(f"""
<div class="vf-hero">
  <div class="vf-logo">🧬</div>
  <div>
    <h1>VaxForge <span class="vf-chip">prototip</span></h1>
    <p class="vf-sub"><b><i>in&nbsp;silico</i> Reverse Vaccinology Pipeline</b></p>
  </div>
</div>
<div class="vf-chips">
  <div class="c">🔬 <b>{_n_tools}</b> analiz aracı</div>
  <div class="c">🧫 <b>{_n_hosts}</b> konak / MHC paneli</div>
  <div class="c">🧠 B-hücre · MHC-I/II · IFN-γ · işleme</div>
  <div class="c">🛡️ alerjenite · toksisite · insan-homoloji</div>
  <div class="c">📄 PDF · HTML · Excel · atıflı rapor</div>
</div>
""", unsafe_allow_html=True)


def _sec(num: str, title: str) -> None:
    """Numaralı profesyonel bölüm başlığı."""
    st.markdown(f'<div class="vf-sec"><div class="n">{num}</div>'
                f'<h3>{title}</h3></div>', unsafe_allow_html=True)

# --- Kenar çubuğu: organizma profili + ortam --------------------------------
with st.sidebar:
    st.markdown('<div class="vf-sec" style="margin:2px 0 10px">'
                '<div class="n">⚙</div><h3>Çalışma ayarları</h3></div>',
                unsafe_allow_html=True)
    profile = st.selectbox(
        "Patojen profili (etken)",
        CFG.profiles,
        index=CFG.profiles.index(CFG.default_profile),
        help="Eşik varsayılanları patojene göre değişir.",
    )
    # Gram-tipi: yalnız bakteride anlamlı (PSORTb doğru bölmeleri Gram'a göre tanır).
    gram = None
    if profile == "bacteria":
        _gram_lbl = st.radio(
            "Gram boyaması",
            ["Gram-negatif", "Gram-pozitif"],
            horizontal=True,
            help="PSORTb lokalizasyonu Gram-tipine göre çalışır. Yanlış seçim = "
                 "yanlış tahmin. Gram−: dış membran/periplazma var; Gram+: hücre duvarı.",
        )
        gram = "negative" if _gram_lbl == "Gram-negatif" else "positive"
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
        pi = "GERÇEK" if netmhc_local.runnable("mhc_i") else "IEDB/proxy"
        pii = "GERÇEK" if netmhc_local.runnable("mhc_ii") else "IEDB/proxy"
        st.caption(f"• {h.label}: MHC-I {len(h.mhc_i)} allel ({pi}), "
                   f"MHC-II {len(h.mhc_ii)} allel ({pii})")
    st.divider()
    st.caption(f"config: `{CFG.path.name}` · hosts: `{HOSTS.path.name}` · "
               f"NetMHCpan: {'yerel' if netmhc_local.available() else 'IEDB/proxy'}")

# IEDB literatür/bilinen-epitop eşleştirmesi otomatik çalışır (organizma seçimi
# olmadan tüm IEDB'ye karşı canlı tarama). Yapısal/docking-MD adımları çıkarıldı
# (odak aday-belirleme). Sidebar sade tutuldu (kullanıcı isteği, 2026-07-10).
organism_taxon = None
has_gpu = False

# --- Girdi: dosya yükleme veya örnek ----------------------------------------
_sec("1", "Girdi dosyası")
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
_sec("2", "Dosya tanıma")
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
_sec("3", "Planlanan pipeline")
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
_sec("4", f"Eşikler — profil: {profile}")
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
_sec("5", "Çalıştır")
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
                               has_gpu=has_gpu, outdir="outputs", host_registry=HOSTS,
                               organism_taxon=organism_taxon, gram=gram):
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
            ("3· Virülans skoru (VFDB)", rmeta.get("n_discovery", "?"), "DIAMOND→VFDB skorlama (eleme DEĞİL; NERVE2 gibi)"),
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

        # --- Aday epitoplar TİPE GÖRE sıralı tablolar (literatür stili) -----
        from vaxforge import evaluate
        tt = evaluate.type_tables(peptides, rmeta, top_n=15)
        total = {k: sum(1 for p in peptides if p.kind == k) for k in ("MHC-I", "MHC-II", "B")}
        st.subheader("🔬 Aday epitoplar — tipe göre sıralı (CTL / HTL / B-hücre)")
        st.caption("Her tip kendi içinde bağlanma gücüne göre sıralı (T-hücre: %rank artan; "
                   "B-hücre: BepiPred azalan). ★ + yeşil satır = tüm zorunlu ölçütleri geçen "
                   "**final seçilen** epitop. ✔=geçti, 📖=IEDB literatürde, Kons.=suş verisi yok. "
                   "Tam araç dökümü **candidates_full.xlsx** indirmesindedir.")

        def _yn(b):
            return "✔" if b else "✗"

        def _n(v, nd=2):
            return "—" if v is None else f"{v:.{nd}f}"

        titles = {"MHC-I": "🟦 CTL — MHC-I / CD8⁺ T-hücre  ·  %rank artan",
                  "MHC-II": "🟨 HTL — MHC-II / CD4⁺ T-hücre  ·  %rank artan",
                  "B": "🟩 B-hücre — antikor  ·  BepiPred azalan"}
        for kind in ("MHC-I", "MHC-II", "B"):
            rows = tt.get(kind, [])
            if not rows:
                continue
            st.markdown(f"**{titles[kind]}**  — {len(rows)}/{total[kind]} gösteriliyor")
            data = []
            for r in rows:
                base = {"★": "★" if r["star"] else "", "Epitop": r["epitope"],
                        "Kaynak": r["source"], "Poz.": r["pos"]}
                if kind == "B":
                    base.update({"Uzun.": r["length"], "BepiPred": _n(r["bcell"]),
                                 "Antij.": _n(r["antigenicity"]), "Alerjen": _yn(r["allergen_ok"]),
                                 "Toksik": _yn(r["toxic_ok"]), "📖": "📖" if r["iedb"] else "–"})
                elif kind == "MHC-I":
                    base.update({"%rank": _n(r["rank"]), "Allel": r["allele"],
                                 "Antij.": _n(r["antigenicity"]), "Alerjen": _yn(r["allergen_ok"]),
                                 "Toksik": _yn(r["toxic_ok"]), "İmmüno.": _n(r["immunogenicity"]),
                                 "İşleme": _n(r["processing"]), "Kons.": "–",
                                 "📖": "📖" if r["iedb"] else "–"})
                else:
                    base.update({"%rank": _n(r["rank"]), "Allel": r["allele"],
                                 "Antij.": _n(r["antigenicity"]), "Alerjen": _yn(r["allergen_ok"]),
                                 "Toksik": _yn(r["toxic_ok"]), "IFN-γ": _yn(r["ifn_ok"]),
                                 "Kons.": "–", "📖": "📖" if r["iedb"] else "–"})
                data.append(base)
            df_t = pd.DataFrame(data)
            sty = df_t.style.apply(
                lambda row: ['background-color:#e7f7ee' if row["★"] == "★" else '' for _ in row],
                axis=1)
            st.dataframe(sty, use_container_width=True, hide_index=True)

        # --- Popülasyon kapsamı (IEDB) ------------------------------------
        popcov = rmeta.get("population_coverage")
        if popcov:
            st.subheader("🌍 Popülasyon kapsamı (IEDB HLA frekansları)")
            if not popcov.get("available"):
                st.info(popcov.get("note", "IEDB Population Coverage aracı kurulu değil."))
            else:
                st.caption("Aday epitop setinin, bir bireyin en az bir epitop-bağlayan "
                           "allele sahip olma olasılığı (%). Yalnız insan HLA için gerçek "
                           "frekans verisi vardır.")
                areas = popcov.get("areas", [])
                for hname, hentry in popcov.get("hosts", {}).items():
                    st.markdown(f"**{hentry.get('label', hname)}**")
                    if hentry.get("note"):
                        st.caption("ℹ️ " + hentry["note"])
                        continue
                    rows = []
                    for cls, lbl in (("mhc_i", "MHC-I"), ("mhc_ii", "MHC-II")):
                        cov = hentry.get(cls)
                        if not cov:
                            continue
                        rows.append({"Sınıf": lbl,
                                     **{a: f"{cov.get(a, {}).get('coverage', '—')}%"
                                        if a in cov else "—" for a in areas}})
                    if rows:
                        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        # --- IEDB literatür/bilinen-epitop taraması + validasyon ----------
        im = rmeta.get("iedb_match")
        if im:
            st.subheader("📖 IEDB literatür / bilinen-epitop taraması")
            if not im.get("available"):
                st.info(im.get("note", "IEDB taraması yapılamadı."))
            else:
                st.caption(f"Kaynak: `{im.get('source','')}` · eşleşen aday: "
                           f"**{im.get('n_matched',0)}/{im.get('n_candidates','?')}**. "
                           "Bir adayın deneysel doğrulanmış epitopla örtüşmesi güçlü pozitif "
                           "kontrol/literatür sinyalidir. Bu adım adaylık puanını **değiştirmez**.")
                irows = []
                for i, p in enumerate(peptides, 1):
                    ie = p.metrics.get("iedb") or {}
                    if ie.get("matched") is not True:
                        continue
                    pmids = ie.get("pmids") or []
                    irows.append({
                        "sıra": i, "peptit": p.seq, "tip": p.kind,
                        "eşleşme": ie.get("match_type", ""),
                        "IEDB epitobu": ie.get("epitope_seq", ""),
                        "organizma": (ie.get("organisms") or ["—"])[0],
                        "antijen": (ie.get("antigens") or ["—"])[0],
                        "PMID": ", ".join(pmids[:4]) or (f"IEDB {ie.get('iedb_id')}" if ie.get("iedb_id") else "—"),
                    })
                if irows:
                    st.dataframe(pd.DataFrame(irows), use_container_width=True, hide_index=True)
                else:
                    st.write("Hiçbir aday bilinen IEDB epitobuyla eşleşmedi (yeni olabilir).")
                bm = im.get("benchmark") or {}
                if bm and bm.get("recall") is not None:
                    st.markdown("**Validasyon — bilinen epitop recall'ü**")
                    st.caption("IEDB'deki deneysel doğrulanmış lineer epitoplar 'ground truth'; "
                               "pipeline tahminleriyle örtüşme üzerinden ölçülür.")
                    b1, b2, b3 = st.columns(3)
                    b1.metric("Bilinen epitop (benzersiz)", bm["n_known"])
                    b2.metric("Yakalanan (recall)", f"{round(bm['recall']*100,1)}%",
                              f"{bm['n_known_hit']}/{bm['n_known']}")
                    b3.metric("Bilinene eşleşen aday",
                              f"{round((bm['precision_like'] or 0)*100,1)}%",
                              f"{bm['n_pred_matched']}/{bm['n_pred']}")

        st.subheader("İndirilebilir çıktılar")
        cols = st.columns(len(paths))
        labels = {"csv": "Adaylar (CSV)", "xlsx": "Tam liste (Excel)",
                  "fasta": "Peptitler (FASTA)", "json": "Tam koşum (JSON)",
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
