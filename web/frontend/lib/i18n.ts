// VaxForge iki dilli metinler (TR varsayılan, EN). Backend'in i18n'iyle uyumlu ton.
export type Lang = "tr" | "en";

export const LANGS: Lang[] = ["tr", "en"];

type Dict = Record<string, { tr: string; en: string }>;

export const STR: Dict = {
  nav_run: { tr: "Pipeline'ı Çalıştır", en: "Run Pipeline" },
  nav_runs: { tr: "Koşular", en: "Runs" },
  nav_method: { tr: "Metodoloji", en: "Methodology" },
  nav_tools: { tr: "Araçlar", en: "Tools" },

  hero_badge: { tr: "Ajan destekli · in silico", en: "Agent-assisted · in silico" },
  hero_title_1: { tr: "Patojenden aşı adayına,", en: "From pathogen to vaccine candidate," },
  hero_title_2: { tr: "tek bir pipeline'da.", en: "in a single pipeline." },
  hero_sub: {
    tr: "Bir patojen dosyası yükleyin (FASTA / FASTQ / GenBank). VaxForge girdiyi otomatik tanır, virülans faktörlerini ve yüzey antijenlerini tarar, B- ve T-hücre epitoplarını tahmin eder, en güçlü peptit adaylarını adaylık puanıyla sıralar ve atıflı bir rapor üretir.",
    en: "Upload a pathogen file (FASTA / FASTQ / GenBank). VaxForge auto-detects the input, mines virulence factors and surface antigens, predicts B- and T-cell epitopes, ranks the strongest peptide candidates by candidacy score, and produces a fully-cited report.",
  },
  hero_cta: { tr: "Analizi Başlat", en: "Start Analysis" },
  hero_cta2: { tr: "Örnek koşuları gör", en: "See example runs" },

  stat_tools: { tr: "gerçek analiz aracı", en: "real analysis tools" },
  stat_hosts: { tr: "konak / MHC paneli", en: "hosts / MHC panels" },
  stat_stages: { tr: "pipeline aşaması", en: "pipeline stages" },
  stat_runs: { tr: "tamamlanmış koşu", en: "completed runs" },

  how_title: { tr: "Nasıl çalışır", en: "How it works" },
  how_sub: {
    tr: "Her adım gerçek bir araç ya da dürüst-etiketli yedek yöntemle koşar; koşulan tüm eşikler rapora yazılır.",
    en: "Each stage runs a real tool or an honestly-labeled fallback; every threshold used is written into the report.",
  },

  tools_title: { tr: "Entegre araçlar", en: "Integrated tools" },
  tools_sub: {
    tr: "Kurulu gerçek araçlar canlı olarak algılanır. ⚠ olanlar için o adımda dürüst-etiketli yedek yöntem kullanılır.",
    en: "Installed real tools are detected live. For those marked ⚠, an honestly-labeled fallback is used for that step.",
  },
  tools_available: { tr: "kurulu", en: "installed" },
  tools_fallback: { tr: "yedek", en: "fallback" },

  cases_title: { tr: "Örnek koşular", en: "Example runs" },
  cases_sub: {
    tr: "Gerçek patojenlerle tamamlanmış analizler — açın, adayları ve raporu inceleyin.",
    en: "Completed analyses on real pathogens — open to inspect candidates and the report.",
  },
  cases_open: { tr: "Koşuyu aç", en: "Open run" },
  cases_candidates: { tr: "aday", en: "candidates" },

  method_title: { tr: "Metodoloji & dürüstlük", en: "Methodology & honesty" },

  footer_note: {
    tr: "VaxForge yalnızca araştırma amaçlıdır — in silico tahminler laboratuvar doğrulaması yerine geçmez.",
    en: "VaxForge is for research use only — in silico predictions do not replace laboratory validation.",
  },

  // run page
  run_title: { tr: "Pipeline'ı çalıştır", en: "Run the pipeline" },
  run_drop: { tr: "Dosyayı buraya sürükleyin ya da seçin", en: "Drop a file here or browse" },
  run_formats: {
    tr: "FASTA · FASTQ · GenBank (.gb/.gbk) · .gz — tip otomatik tanınır",
    en: "FASTA · FASTQ · GenBank (.gb/.gbk) · .gz — type auto-detected",
  },
  run_profile: { tr: "Patojen profili", en: "Pathogen profile" },
  run_gram: { tr: "Gram tipi", en: "Gram type" },
  run_hosts: { tr: "Konaklar / MHC paneli", en: "Hosts / MHC panel" },
  run_start: { tr: "Analizi başlat", en: "Start analysis" },
  run_running: { tr: "Pipeline çalışıyor…", en: "Pipeline running…" },
  run_console: { tr: "Canlı ilerleme", en: "Live progress" },
  run_detect: { tr: "Girdi tanıma", en: "Input detection" },
  run_done: { tr: "Tamamlandı", en: "Completed" },
  run_error: { tr: "Hata", en: "Error" },
  run_pick: { tr: "Başlamak için bir dosya seçin.", en: "Pick a file to begin." },

  res_flow: { tr: "Eleme akışı", en: "Elimination funnel" },
  res_candidates: { tr: "Sıralı aday peptitler", en: "Ranked candidate peptides" },
  res_downloads: { tr: "İndir", en: "Downloads" },
  res_report: { tr: "Tam rapor", en: "Full report" },
  col_rank: { tr: "#", en: "#" },
  col_seq: { tr: "Peptit", en: "Peptide" },
  col_type: { tr: "Tip", en: "Type" },
  col_parent: { tr: "Kaynak", en: "Parent" },
  col_score: { tr: "Adaylık", en: "Candidacy" },
};

export function t(lang: Lang, key: string): string {
  const e = STR[key];
  if (!e) return key;
  return e[lang] ?? e.tr;
}
