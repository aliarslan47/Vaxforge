"""Basit iki-dilli (TR/EN) metin sözlüğü + t() yardımcısı.

Kullanım: t(lang, "anahtar")  → lang ∈ {"tr","en"}. Anahtar yoksa TR'ye,
o da yoksa anahtarın kendisine düşer (dürüst, çökmez). Bilimsel veri (peptit,
allel, sayılar, 'GERÇEK'/'proxy' method etiketleri) çevrilmez — yalnız arayüz/
rapor iskeleti metinleri.
"""

from __future__ import annotations

STRINGS: dict[str, dict[str, str]] = {
    # ---- Ortak / arayüz ----
    "app_tagline":        {"tr": "in silico Reverse Vaccinology Pipeline",
                           "en": "in silico Reverse Vaccinology Pipeline"},
    "chip_tools":         {"tr": "analiz aracı", "en": "analysis tools"},
    "chip_hosts":         {"tr": "konak / MHC paneli", "en": "hosts / MHC panel"},
    "chip_epi":           {"tr": "B-hücre · MHC-I/II · IFN-γ · işleme",
                           "en": "B-cell · MHC-I/II · IFN-γ · processing"},
    "chip_filters":       {"tr": "alerjenite · toksisite · insan-homoloji",
                           "en": "allergenicity · toxicity · self-homology"},
    "chip_report":        {"tr": "PDF · HTML · Excel · atıflı rapor",
                           "en": "PDF · HTML · Excel · cited report"},
    "sec_input":          {"tr": "Girdi dosyası", "en": "Input file"},
    "sec_thresholds":     {"tr": "Eşikler — profil", "en": "Thresholds — profile"},
    "settings":           {"tr": "Çalışma ayarları", "en": "Run settings"},
    "language":           {"tr": "Dil / Language", "en": "Language / Dil"},
    "pathogen_profile":   {"tr": "Patojen profili (etken)", "en": "Pathogen profile"},
    "gram_stain":         {"tr": "Gram boyaması", "en": "Gram stain"},
    "gram_neg":           {"tr": "Gram-negatif", "en": "Gram-negative"},
    "gram_pos":           {"tr": "Gram-pozitif", "en": "Gram-positive"},
    "hosts_label":        {"tr": "Konak(lar) — MHC taranacak", "en": "Host(s) — MHC scanned"},
    "upload_hint":        {"tr": "Dosyanızı sürükleyip bırakın — FASTA / FASTQ / GenBank (.gz destekli)",
                           "en": "Drag & drop your file — FASTA / FASTQ / GenBank (.gz supported)"},
    "try_sample":         {"tr": "ya da örnek dosyayla dene:", "en": "or try a sample file:"},
    "sample_proteome_vfdb": {"tr": "Patojen proteinleri (gerçek VFDB)",
                             "en": "Pathogen proteins (real VFDB)"},
    "sample_proteome":    {"tr": "Proteom (protein FASTA)", "en": "Proteome (protein FASTA)"},
    "sample_cds":         {"tr": "Genler / CDS (nükleotid)", "en": "Genes / CDS (nucleotide)"},
    "sample_reads":       {"tr": "Ham okumalar (FASTQ)", "en": "Raw reads (FASTQ)"},
    "pick_file":          {"tr": "Bir dosya yükleyin ya da sağdaki örneklerden birini seçin. "
                                 "Desteklenen: protein/nükleotid FASTA veya FASTQ (.gz olabilir).",
                           "en": "Upload a file or pick one of the samples on the right. "
                                 "Supported: protein/nucleotide FASTA or FASTQ (.gz allowed)."},
    "run_btn":            {"tr": "▶ Pipeline'ı çalıştır", "en": "▶ Run pipeline"},
    "running":            {"tr": "Pipeline çalışıyor…", "en": "Pipeline running…"},
    "result_header":      {"tr": "Sonuç — en iyi aday peptitler", "en": "Result — top candidate peptides"},
    "downloads":          {"tr": "İndirmeler", "en": "Downloads"},
    "sec_detect":         {"tr": "Dosya tanıma", "en": "File detection"},
    "sec_plan":           {"tr": "Planlanan pipeline", "en": "Planned pipeline"},
    "sec_run":            {"tr": "Çalıştır", "en": "Run"},
    "run_start":          {"tr": "▶ Pipeline'ı başlat", "en": "▶ Start pipeline"},
    "res_flow":           {"tr": "📊 Eleme akışı (kaç girdi → kaç kaldı)",
                           "en": "📊 Filtering flow (input → survivors)"},
    "res_bytype":         {"tr": "🔬 Aday epitoplar — tipe göre sıralı (CTL / HTL / B-hücre)",
                           "en": "🔬 Candidate epitopes — ranked by type (CTL / HTL / B-cell)"},
    "res_pop":            {"tr": "🌍 Popülasyon kapsamı (IEDB HLA frekansları)",
                           "en": "🌍 Population coverage (IEDB HLA frequencies)"},
    "res_iedb":           {"tr": "📖 IEDB literatür / bilinen-epitop taraması",
                           "en": "📖 IEDB literature / known-epitope screening"},
    "res_downloads":      {"tr": "İndirilebilir çıktılar", "en": "Downloadable outputs"},
    "pick_info":          {"tr": "👆 Bir dosya yükleyin ya da sağdaki örneklerden birini seçin. "
                                 "Desteklenen: protein/nükleotid FASTA veya FASTQ (.gz olabilir).",
                           "en": "👆 Upload a file or pick a sample on the right. "
                                 "Supported: protein/nucleotide FASTA or FASTQ (.gz allowed)."},

    # ---- Rapor / tip tabloları (HTML + PDF ortak) ----
    "rep_title":          {"tr": "VaxForge — Aşı Adayı Raporu", "en": "VaxForge — Vaccine Candidate Report"},
    "rep_generated":      {"tr": "Oluşturma", "en": "Generated"},
    "rep_input":          {"tr": "Girdi", "en": "Input"},
    "rep_profile":        {"tr": "Patojen profili", "en": "Pathogen profile"},
    "rep_hosts":          {"tr": "Konak(lar)", "en": "Host(s)"},
    "rep_summary":        {"tr": "Özet", "en": "Summary"},
    "rep_sum_input":      {"tr": "Girdi proteini", "en": "Input proteins"},
    "rep_sum_discovery":  {"tr": "Keşif sonrası", "en": "After discovery"},
    "rep_sum_funnel":     {"tr": "Huni sonrası", "en": "After funnel"},
    "rep_sum_epitope":    {"tr": "Epitop", "en": "Epitopes"},
    "rep_sum_survivors":  {"tr": "Sağ kalan aday", "en": "Surviving candidates"},
    "rep_top15":          {"tr": "En iyi 15 aday peptit (özet)", "en": "Top 15 candidate peptides (summary)"},
    "rep_bytype":         {"tr": "Aday epitoplar — tipe göre sıralı (CTL / HTL / B-hücre)",
                           "en": "Candidate epitopes — ranked by type (CTL / HTL / B-cell)"},
    "rep_thresholds":     {"tr": "Kullanılan eşikler (tekrarlanabilirlik)",
                           "en": "Thresholds used (reproducibility)"},
    "rep_plan":           {"tr": "Çalıştırılan plan", "en": "Executed plan"},
    "rep_methods":        {"tr": "Kullanılan yöntemler ve araçlar", "en": "Methods and tools used"},
    "rep_references":     {"tr": "Referanslar", "en": "References"},
    "tt_ctl":             {"tr": "CTL — MHC-I / CD8⁺ T-hücre", "en": "CTL — MHC-I / CD8⁺ T-cell"},
    "tt_htl":             {"tr": "HTL — MHC-II / CD4⁺ T-hücre", "en": "HTL — MHC-II / CD4⁺ T-cell"},
    "tt_bcell":           {"tr": "B-hücre — antikor", "en": "B-cell — antibody"},
    "tt_rank_asc":        {"tr": "bağlanma %rank artan", "en": "binding %rank ascending"},
    "tt_bepi_desc":       {"tr": "BepiPred azalan", "en": "BepiPred descending"},
    "tt_showing":         {"tr": "gösteriliyor", "en": "shown"},
    "tt_intro":           {"tr": "Epitoplar tipe göre ayrı tablolarda, her tip kendi içinde bağlanma "
                                 "gücüne göre sıralı. ★ + yeşil satır = tüm zorunlu ölçütleri geçen final "
                                 "seçilen epitop (antijenik + güçlü bağlanma [+ HTL'de IFN-γ]). ✔=geçti, "
                                 "📖=IEDB literatürde, Kons.=suş verisi yok. Tam araç dökümü "
                                 "candidates_full.xlsx dosyasındadır.",
                           "en": "Epitopes are shown in separate tables per type, each ranked by binding "
                                 "strength. ★ + green row = final selected epitope passing all mandatory "
                                 "criteria (antigenic + strong binding [+ IFN-γ for HTL]). ✔=pass, "
                                 "📖=in IEDB literature, Cons.=no strain data. Full per-tool dump is in "
                                 "candidates_full.xlsx."},
    # tablo kolon başlıkları
    "col_epitope":        {"tr": "Epitop", "en": "Epitope"},
    "col_source":         {"tr": "Kaynak antijen", "en": "Source antigen"},
    "col_pos":            {"tr": "Poz.", "en": "Pos."},
    "col_len":            {"tr": "Uzun.", "en": "Len."},
    "col_allele":         {"tr": "Allel", "en": "Allele"},
    "col_antig":          {"tr": "Antij.", "en": "Antig."},
    "col_allergen":       {"tr": "Alerjen", "en": "Allergen"},
    "col_toxic":          {"tr": "Toksik", "en": "Toxic"},
    "col_immuno":         {"tr": "İmmüno.", "en": "Immuno."},
    "col_proc":           {"tr": "İşleme", "en": "Processing"},
    "col_cons":           {"tr": "Kons.", "en": "Cons."},
    "col_step":           {"tr": "Adım", "en": "Step"},
    "col_tool":           {"tr": "Araç", "en": "Tool"},
    "col_param":          {"tr": "Parametre", "en": "Parameter"},
    "col_value":          {"tr": "Değer", "en": "Value"},
    "col_type":           {"tr": "Tip", "en": "Type"},
    "col_status":         {"tr": "Durum", "en": "Status"},
    "col_note":           {"tr": "Not", "en": "Note"},
    "type_hard":          {"tr": "sert", "en": "hard"},
    "type_score":         {"tr": "skor", "en": "score"},
    "rep_disclaimer":     {"tr": "Bu rapor VaxForge in silico reverse vaccinology hattı ile üretilmiştir. "
                                 "Kullanılan araçlar ve eşikler aşağıda listelenir; yöntem etiketleri "
                                 "(GERÇEK/proxy) her tabloda gösterilir.",
                           "en": "This report was produced by the VaxForge in silico reverse vaccinology "
                                 "pipeline. Tools and thresholds are listed below; method labels "
                                 "(real/proxy) are shown in each table."},
    # ---- IEDB bölümü ----
    "iedb_title":         {"tr": "IEDB literatür/bilinen-epitop taraması",
                           "en": "IEDB literature / known-epitope screening"},
    "iedb_none":          {"tr": "Hiçbir aday bilinen IEDB epitobuyla eşleşmedi.",
                           "en": "No candidate matched a known IEDB epitope."},
    "iedb_intro":         {"tr": "Kaynak: {src} · eşleşen aday: <b>{n}/{tot}</b>. Eşleşme, adayın deneysel "
                                 "doğrulanmış bir epitopla (exact/içerme/ortak çekirdek) örtüştüğünü gösterir "
                                 "— güçlü pozitif kontrol sinyali. Bu adım adaylık puanını <b>değiştirmez</b> "
                                 "(salt yorumlama).",
                           "en": "Source: {src} · matched candidates: <b>{n}/{tot}</b>. A match means the "
                                 "candidate overlaps an experimentally validated epitope (exact/containment/"
                                 "shared core) — a strong positive-control signal. This step does <b>not</b> "
                                 "change the candidacy score (interpretation only)."},
    "iedb_val_title":     {"tr": "Validasyon — bilinen epitop recall'ü",
                           "en": "Validation — known-epitope recall"},
    "iedb_val_text":      {"tr": "Bu organizma için IEDB'de deneysel doğrulanmış lineer epitoplar 'ground "
                                 "truth' alınır; pipeline'ın tahminleri bunlarla örtüşme (≥{k} aa ortak "
                                 "çekirdek / içerme / exact) üzerinden değerlendirilir. NOT: pipeline seçici "
                                 "olarak KISA bir öncelik listesi üretir; bu yüzden tüm epitop kataloğuna "
                                 "karşı recall doğası gereği düşüktür — asıl anlamlı ölçüt, adayların ne "
                                 "kadarının deneysel doğrulanmış olduğudur (eşleşme oranı).",
                           "en": "For this organism, experimentally validated linear epitopes in IEDB are "
                                 "taken as ground truth; the pipeline's predictions are evaluated by overlap "
                                 "(≥{k} aa shared core / containment / exact). NOTE: the pipeline "
                                 "deliberately produces a SHORT priority list; thus recall against the full "
                                 "epitope catalog is inherently low — the meaningful metric is what fraction "
                                 "of candidates are experimentally validated (match rate)."},
    "col_order":          {"tr": "Sıra", "en": "Rank"},
    "col_candidate":      {"tr": "Aday peptit", "en": "Candidate peptide"},
    "col_match":          {"tr": "Eşleşme", "en": "Match"},
    "col_known":          {"tr": "Bilinen epitop (benzersiz)", "en": "Known epitopes (unique)"},
    "col_captured":       {"tr": "Yakalanan", "en": "Captured"},
    "col_recall":         {"tr": "Recall", "en": "Recall"},
    "col_pred":           {"tr": "Tahmin (aday)", "en": "Predicted (candidates)"},
    "col_pred_match":     {"tr": "Bilinene eşleşen aday", "en": "Candidates matching known"},
    "col_matchrate":      {"tr": "Eşleşme oranı", "en": "Match rate"},
    # ---- Popülasyon bölümü ----
    "pop_title":          {"tr": "Popülasyon kapsamı (IEDB HLA frekansları)",
                           "en": "Population coverage (IEDB HLA frequencies)"},
    "pop_text":           {"tr": "Bir bireyin en az bir epitop-bağlayan allele sahip olma olasılığı (%). "
                                 "Gerçek frekans verisi yalnız insan HLA için mevcuttur.",
                           "en": "Probability (%) that an individual carries at least one epitope-binding "
                                 "allele. Real frequency data is available only for human HLA."},
    "col_hostclass":      {"tr": "Konak / sınıf", "en": "Host / class"},
    "col_region":         {"tr": "Bölge", "en": "Region"},
    "col_coverage":       {"tr": "Kapsam", "en": "Coverage"},
    "col_iedb_epi":       {"tr": "IEDB epitobu", "en": "IEDB epitope"},
    "col_organism":       {"tr": "Organizma", "en": "Organism"},
    "col_antigen":        {"tr": "Antijen", "en": "Antigen"},
    "col_reference":      {"tr": "Referans", "en": "Reference"},
    "iedb_unavail":       {"tr": "IEDB taraması yapılamadı", "en": "IEDB screening unavailable"},
    # ---- PDF'e özgü (numaralı başlıklar + paragraflar) ----
    "pdf_top":            {"tr": "En iyi aday peptitler (CDS kaynağı + lokus ile)",
                           "en": "Top candidate peptides (with CDS source + locus)"},
    "pdf_anchor_title":   {"tr": "MHC yarığı anchor/cep motifi (yorum — sıralamayı etkilemez)",
                           "en": "MHC groove anchor/pocket motif (interpretive — does not affect ranking)"},
    "pdf_anchor_text":    {"tr": "Peptidin anchor kalıntıları (P2, C-terminal PΩ) ve o allelin NetMHCpan "
                                 "taramasından AMPİRİK çıkarılan cep tercihi. Bağlanma uyumu %rank'ta zaten "
                                 "puanlanır; bu tablo yalnız yorum içindir.",
                           "en": "The peptide's anchor residues (P2, C-terminal PΩ) and the allele's pocket "
                                 "preference derived EMPIRICALLY from a NetMHCpan scan. Binding fit is "
                                 "already scored via %rank; this table is interpretation only."},
    "pdf_methods_note":   {"tr": "GPU gerektiren yapısal adımlar (AlphaFold peptit-MHC, moleküler dinamik) "
                                 "bu prototipte çıkarılmıştır; odak aday-belirlemedir.",
                           "en": "GPU-dependent structural steps (AlphaFold peptide-MHC, molecular dynamics) "
                                 "are removed in this prototype; the focus is candidate identification."},
    "pop_text_pdf":       {"tr": "Aday epitop setinin, bir bireyin en az bir epitop-bağlayan allele sahip "
                                 "olma olasılığını (%) IEDB HLA frekanslarından hesaplar. Gerçek frekans "
                                 "verisi yalnız insan HLA için mevcuttur; diğer konaklar dürüstçe atlanır.",
                           "en": "Computes the probability (%) that an individual carries at least one "
                                 "epitope-binding allele, from IEDB HLA frequencies. Real frequency data "
                                 "exists only for human HLA; other hosts are honestly skipped."},
    "tt_intro_pdf":       {"tr": "Epitoplar tipe göre ayrı tablolarda, her tip kendi içinde bağlanma gücüne "
                                 "göre sıralı (T-hücre: %rank artan; B-hücre: BepiPred azalan). ★ + yeşil "
                                 "satır = tüm zorunlu ölçütleri geçen final seçilen epitop. ✓=geçti, "
                                 "Lit ✓=IEDB literatürde, Kons.=suş verisi yok.",
                           "en": "Epitopes shown per type, each ranked by binding strength (T-cell: %rank "
                                 "ascending; B-cell: BepiPred descending). ★ + green row = final selected "
                                 "epitope passing all mandatory criteria. ✓=pass, Lit ✓=in IEDB literature, "
                                 "Cons.=no strain data."},
    "col_anchors":        {"tr": "Anchorlar", "en": "Anchors"},
    "col_pocket":         {"tr": "Allel cep tercihi", "en": "Allele pocket pref."},
    "col_match2":         {"tr": "Uyum", "en": "Match"},
    "disc_short":         {"tr": "gösteriliyor", "en": "shown"},
}

# Dinamik 'method' etiketlerini (pipeline TR üretir) EN'e çeviren token haritası.
_METHOD_TOKENS = [
    ("GERÇEK", "REAL"), ("heuristik proxy", "heuristic proxy"), ("heuristik", "heuristic"),
    ("ÇALIŞTIRILMADI", "NOT RUN"), ("çalıştırılmadı", "not run"),
    ("yerel", "local"), ("çevrimdışı", "offline"), ("yedek", "fallback"),
    ("yok", "n/a"), ("klasik yöntemler", "classical methods"), ("klasik yöntem", "classical method"),
    ("İnsan", "Human"), ("Fare", "Mouse"), ("Sığır", "Cattle"), ("Domuz", "Pig"), ("Tavuk", "Chicken"),
    ("insan", "human"), ("Homo sapiens", "Homo sapiens"),
    ("VaxiJen DEĞİL", "not VaxiJen"), ("DEĞİL", "NOT"), ("resmi model değil", "not the official model"),
    ("insan Swiss-Prot", "human Swiss-Prot"), ("proteom", "proteome"),
]


def method_label(s: str, lang: str) -> str:
    """TR üretilmiş 'method' etiketini EN'e çevirir (token değişimi). lang!=en → aynen."""
    if lang != "en" or not s:
        return s
    for tr, en in _METHOD_TOKENS:
        s = s.replace(tr, en)
    return s


def t(lang: str, key: str, default: str | None = None) -> str:
    entry = STRINGS.get(key)
    if not entry:
        return default if default is not None else key
    return entry.get(lang) or entry.get("tr") or key
