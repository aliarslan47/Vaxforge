"use client";

import { motion, useReducedMotion } from "framer-motion";
import {
  FileInput,
  Radar,
  FilterX,
  Grid2x2,
  ListOrdered,
  FileText,
} from "lucide-react";
import { useLang } from "./lang-provider";

const STAGES = [
  {
    key: "ingest",
    icon: FileInput,
    tr: { t: "Girdi & Tanıma", d: "FASTA/FASTQ/GenBank otomatik tanınır, proteine çevrilir" },
    en: { t: "Ingest & Detect", d: "FASTA/FASTQ/GenBank auto-detected, translated to protein" },
  },
  {
    key: "discovery",
    icon: Radar,
    tr: { t: "Virülans Keşfi", d: "VFDB / anahtar-kelime ile virülans skoru (NERVE2 gibi)" },
    en: { t: "Virulence Discovery", d: "Virulence scoring via VFDB / keywords (NERVE2-style)" },
  },
  {
    key: "funnel",
    icon: FilterX,
    tr: { t: "Antijen Hunisi", d: "DeepLoc + TMHMM + SignalP + IApred + insan-homoloji" },
    en: { t: "Antigen Funnel", d: "DeepLoc + TMHMM + SignalP + IApred + self-homology" },
  },
  {
    key: "epitope",
    icon: Grid2x2,
    tr: { t: "Epitop Tahmini", d: "B-hücre · MHC-I/II · IFN-γ · işleme, seçili konaklarda" },
    en: { t: "Epitope Prediction", d: "B-cell · MHC-I/II · IFN-γ · processing across hosts" },
  },
  {
    key: "scoring",
    icon: ListOrdered,
    tr: { t: "Eleme & Sıralama", d: "Alerjenite · toksisite süzgeci, adaylık puanıyla sıralama" },
    en: { t: "Filter & Rank", d: "Allergenicity · toxicity filters, ranked by candidacy" },
  },
  {
    key: "report",
    icon: FileText,
    tr: { t: "Atıflı Rapor", d: "HTML · PDF · Excel — koşulan tüm eşikler ve atıflar" },
    en: { t: "Cited Report", d: "HTML · PDF · Excel — every threshold and citation" },
  },
];

export function PipelineDiagram() {
  const { lang } = useLang();
  const reduce = useReducedMotion();

  return (
    <div className="relative">
      {/* akış çizgisi (yatay, geniş ekran) */}
      <div className="pointer-events-none absolute left-0 right-0 top-[34px] hidden h-px bg-gradient-to-r from-transparent via-line to-transparent lg:block" />
      <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-6 lg:gap-3">
        {STAGES.map((s, i) => {
          const txt = s[lang];
          const Icon = s.icon;
          return (
            <motion.div
              key={s.key}
              initial={reduce ? false : { opacity: 0, y: 16 }}
              whileInView={{ opacity: 1, y: 0 }}
              viewport={{ once: true, margin: "-60px" }}
              transition={{ duration: 0.4, delay: i * 0.06, ease: [0.16, 1, 0.3, 1] }}
              className="group relative"
            >
              <div className="flex h-full flex-col rounded-2xl border border-line bg-surface/50 p-4 backdrop-blur-sm transition-colors hover:border-primary/40">
                <div className="mb-3 flex items-center gap-2">
                  <span className="grid h-9 w-9 shrink-0 place-items-center rounded-xl bg-primary/10 text-primary ring-1 ring-primary/20 transition group-hover:bg-primary/20">
                    <Icon className="h-[18px] w-[18px]" strokeWidth={2} />
                  </span>
                  <span className="font-mono text-[11px] text-fg-faint">
                    {String(i + 1).padStart(2, "0")}
                  </span>
                </div>
                <h3 className="font-display text-sm font-semibold text-fg">{txt.t}</h3>
                <p className="mt-1.5 text-[12.5px] leading-relaxed text-fg-muted">{txt.d}</p>
              </div>
            </motion.div>
          );
        })}
      </div>
    </div>
  );
}
