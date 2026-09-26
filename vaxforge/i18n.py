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
    "iedb_intro":         {"tr": "Kaynak: {src} · doğrulanmış eşleşme: <b>{n_val}/{tot}</b> "
                                 "(exact veya kısa tarafın tam-içerilmesi), ayrıca zayıf benzerlik: "
                                 "<b>{n_weak}</b> (yalnız ortak çekirdek). Yalnız <b>doğrulanmış</b> eşleşmeler "
                                 "pozitif kontrol sayılır; zayıf benzerlik tam bir bilinen epitobu yeniden "
                                 "üretmez. Bu adım adaylık puanını <b>değiştirmez</b> (salt yorumlama).",
                           "en": "Source: {src} · validated matches: <b>{n_val}/{tot}</b> "
                                 "(exact or full containment of the shorter side), plus weak similarity: "
                                 "<b>{n_weak}</b> (shared core only). Only <b>validated</b> matches count as "
                                 "positive controls; a weak similarity does not reproduce a full known epitope. "
                                 "This step does <b>not</b> change the candidacy score (interpretation only)."},
    "iedb_val_title":     {"tr": "Validasyon — bilinen epitop recall'ü",
                           "en": "Validation — known-epitope recall"},
    "iedb_val_text":      {"tr": "Bu organizma için IEDB'de deneysel doğrulanmış lineer epitoplar 'ground "
                                 "truth' alınır. İki recall gösterilir: gevşek (≥{k} aa ortak çekirdek / "
                                 "içerme / exact) ve <b>tam</b> (yalnız exact/içerme). NOT: pipeline seçici "
                                 "olarak KISA bir öncelik listesi üretir; bu yüzden tüm epitop kataloğuna "
                                 "karşı recall doğası gereği düşüktür. Ayrıca zemin tek bir referansa "
                                 "yığılabilir (aşağıdaki uyarıya bakın), bu durumda ölçüt döngüselliğe eğilimlidir.",
                           "en": "For this organism, experimentally validated linear epitopes in IEDB are "
                                 "taken as ground truth. Two recalls are shown: loose (≥{k} aa shared core / "
                                 "containment / exact) and <b>exact</b> (exact/containment only). NOTE: the "
                                 "pipeline deliberately produces a SHORT priority list; thus recall against "
                                 "the full catalog is inherently low. The ground truth may also be dominated "
                                 "by a single reference (see caveat below), which biases the metric toward "
                                 "circularity."},
    "col_order":          {"tr": "Sıra", "en": "Rank"},
    "col_candidate":      {"tr": "Aday peptit", "en": "Candidate peptide"},
    "col_match":          {"tr": "Eşleşme", "en": "Match"},
    "col_known":          {"tr": "Bilinen epitop (benzersiz)", "en": "Known epitopes (unique)"},
    "col_captured":       {"tr": "Yakalanan", "en": "Captured"},
    "col_recall":         {"tr": "Recall", "en": "Recall"},
    "col_pred":           {"tr": "Tahmin (aday)", "en": "Predicted (candidates)"},
    "col_pred_match":     {"tr": "Bilinene eşleşen aday", "en": "Candidates matching known"},
    "col_matchrate":      {"tr": "Eşleşme oranı", "en": "Match rate"},
    "col_strength":       {"tr": "Kanıt gücü", "en": "Evidence"},
    "col_overlap":        {"tr": "Örtüşme (aa)", "en": "Overlap (aa)"},
    "col_assay":          {"tr": "Assay", "en": "Assay"},
    "col_recall_exact":   {"tr": "Recall (tam)", "en": "Recall (exact)"},
    "strength_validated": {"tr": "Tam", "en": "Validated"},
    "strength_weak":      {"tr": "Kısmi", "en": "Weak"},
    "iedb_caveat":        {"tr": "Not: Eşleşme, IEDB'nin ilgili referans için yatırdığı kayda dayanır; "
                                 "dizi makalenin basılı tablosunda görünmeyebilir — IEDB kayıt kimliğinden "
                                 "doğrulayın.",
                           "en": "Note: A match is based on IEDB's deposited record for the cited reference; "
                                 "the sequence may not appear in the article's printed table — verify via the "
                                 "IEDB entry ID."},
    "iedb_ref_dominance": {"tr": "Zemin uyarısı: bu organizmanın IEDB epitoplarının {share}%'i tek referanstan "
                                 "(PMID {pmid}) — recall/eşleşme oranı bu kaynağa göre yorumlanmalı.",
                           "en": "Ground-truth caveat: {share}% of this organism's IEDB epitopes come from a "
                                 "single reference (PMID {pmid}) — recall/match rate should be read accordingly."},
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
    "strain_hint":        {"tr": "Konservasyon (opsiyonel): ek suş proteom(lar)ı — çok-suş FASTA",
                           "en": "Conservation (optional): extra strain proteome(s) — multi-strain FASTA"},
    "pdf_immuno_curve":   {"tr": "Simüle antikor yanıtı — izotip (bağıl, tepe = 1)",
                           "en": "Simulated antibody response — isotype (relative, peak = 1)"},
    "pdf_immuno_day":     {"tr": "gün", "en": "day"},
    "pdf_immuno_relab":   {"tr": "bağıl antikor düzeyi", "en": "relative antibody level"},
    "pdf_bio_limits_title": {"tr": "Biyolojik kapsam ve sınırlamalar", "en": "Biological scope and limitations"},
    "pdf_bcell_caveat":   {"tr": "B-hücre epitop tahmini yalnız LİNEER'dir (BepiPred-3.0). B-hücre "
                                 "epitoplarının ~%90'ı süreksiz/konformasyoneldir (Michalik ve ark. 2016); "
                                 "yapı-tabanlı konformasyonel tahmin (DiscoTope-3), AlphaFold/GPU adımı "
                                 "canlandığında eklenecektir. Final ürün lineer-peptit MEV olduğundan bu "
                                 "sınır konstruktu değil, tam-protein antijen sıralamasını etkiler.",
                           "en": "B-cell epitope prediction is LINEAR only (BepiPred-3.0). ~90% of B-cell "
                                 "epitopes are discontinuous/conformational (Michalik et al. 2016); "
                                 "structure-based conformational prediction (DiscoTope-3) will be added when "
                                 "the AlphaFold/GPU step is live. As the final product is a linear-peptide MEV, "
                                 "this limitation affects whole-protein antigen ranking, not the construct."},
    "pdf_cons_computed":  {"tr": "Suş konservasyonu: {n} suş üzerinden hesaplandı; {c} aday antijen ≥%{m} "
                                 "korunmuş (geniş-suş koruması için yüksek istenir).",
                           "en": "Strain conservation: computed over {n} strains; {c} candidate antigens "
                                 "≥{m}% conserved (high is preferred for broad-strain protection)."},
    "pdf_cons_skipped":   {"tr": "Suş konservasyonu: ek suş proteomu verilmediği için hesaplanmadı "
                                 "(uydurma yok). Çok-suş girdisi verilirse korunmuşluk skoru devreye girer.",
                           "en": "Strain conservation: not computed (no extra strain proteomes provided). "
                                 "Providing multi-strain input enables the conservation score."},
    "pdf_mimicry_note":   {"tr": "Moleküler mimikri: {n} aday epitop konağın kendi proteomundaki bir 9-mer'e "
                                 "birebir uydu (otoimmünite/tolerans riski) — YUMUŞAK işaret (elenmedi, skor "
                                 "×0.9). Protein-seviyesi konak homolojisi funnel'da zaten sert filtredir.",
                           "en": "Molecular mimicry: {n} candidate epitopes exactly match a 9-mer in the "
                                 "host's own proteome (autoimmunity/tolerance risk) — SOFT flag (not removed, "
                                 "score ×0.9). Protein-level host homology is already a hard filter in the funnel."},
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
    # ---- Kullanılan yöntemler (funnel method özeti — Faz 1) ----
    "rep_methods":        {"tr": "Kullanılan yöntemler", "en": "Methods used"},
    "rep_methods_sub":    {"tr": "Her adımda GERÇEKTE çalışan araç. Araç kuruluysa "
                                 "\"GERÇEK\", değilse dürüst vekil/heuristik gösterilir "
                                 "(uydurma yok).",
                           "en": "The tool that ACTUALLY ran at each step. \"REAL\" if the "
                                 "tool is installed, otherwise an honest proxy/heuristic "
                                 "(no fabrication)."},
    "meth_step":          {"tr": "Adım", "en": "Step"},
    "meth_method":        {"tr": "Yöntem", "en": "Method"},
    "meth_localization":  {"tr": "Hücre-altı lokalizasyon", "en": "Subcellular localization"},
    "meth_tm":            {"tr": "Transmembran topoloji", "en": "Transmembrane topology"},
    "meth_signal":        {"tr": "Sinyal peptidi", "en": "Signal peptide"},
    "meth_antigenicity":  {"tr": "Antijenite", "en": "Antigenicity"},
    "meth_homology":      {"tr": "Öz-homoloji (otoimmünite)", "en": "Self-homology (autoimmunity)"},
}

# Dinamik 'method' etiketlerini (pipeline TR üretir) EN'e çeviren token haritası.
_METHOD_TOKENS = [
    ("GERÇEK", "REAL"), ("heuristik proxy", "heuristic proxy"), ("heuristik", "heuristic"),
    ("ÇALIŞTIRILMADI", "NOT RUN"), ("çalıştırılmadı", "not run"),
    ("yerel", "local"), ("çevrimdışı", "offline"), ("yedek", "fallback"),
    ("sert:", "hard:"), ("yumuşak:", "soft:"),
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
