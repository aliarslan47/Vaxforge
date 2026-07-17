"use client";

import { useEffect, useState } from "react";
import { motion, useReducedMotion } from "framer-motion";
import {
  ArrowRight,
  ShieldCheck,
  FlaskConical,
  Braces,
  BookOpenCheck,
} from "lucide-react";
import { Navbar } from "@/components/navbar";
import { PipelineDiagram } from "@/components/pipeline-diagram";
import { ToolsGrid } from "@/components/tools-grid";
import { CaseStudies } from "@/components/case-studies";
import { Button, Card, Badge, SectionHeading } from "@/components/ui";
import { useLang } from "@/components/lang-provider";
import { getConfig, getRuns, AppConfig, RunSummary } from "@/lib/api";

export default function Home() {
  const { t } = useLang();
  const reduce = useReducedMotion();
  const [config, setConfig] = useState<AppConfig | null>(null);
  const [runs, setRuns] = useState<RunSummary[]>([]);

  useEffect(() => {
    getConfig().then(setConfig).catch(() => {});
    getRuns().then(setRuns).catch(() => {});
  }, []);

  const nTools = config?.steps?.filter((x) => x.installed).length ?? config?.tools.filter((x) => x.available).length ?? 0;
  const nHosts = config?.hosts.length ?? 0;
  const stats = [
    { v: nTools ? `${nTools}` : "—", l: t("stat_tools") },
    { v: nHosts ? `${nHosts}` : "—", l: t("stat_hosts") },
    { v: "6", l: t("stat_stages") },
    { v: runs.length ? `${runs.length}` : "—", l: t("stat_runs") },
  ];

  return (
    <main className="relative min-h-dvh overflow-x-hidden">
      <Navbar />

      {/* ---------- HERO ---------- */}
      <section className="relative px-5 pt-32 pb-20 sm:pt-40">
        <div className="grid-texture pointer-events-none absolute inset-0 -z-10" />
        <div className="mx-auto max-w-4xl text-center">
          <motion.div
            initial={reduce ? false : { opacity: 0, y: 12 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
          >
            <Badge tone="primary" className="mb-6">
              <span className="h-1.5 w-1.5 animate-pulse-dot rounded-full bg-primary" />
              {t("hero_badge")}
            </Badge>
          </motion.div>

          <motion.h1
            initial={reduce ? false : { opacity: 0, y: 16 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.55, delay: 0.05 }}
            className="font-display text-4xl font-semibold leading-[1.08] tracking-tight text-fg sm:text-6xl"
          >
            {t("hero_title_1")}
            <br />
            <span className="text-gradient">{t("hero_title_2")}</span>
          </motion.h1>

          <motion.p
            initial={reduce ? false : { opacity: 0, y: 16 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.55, delay: 0.12 }}
            className="mx-auto mt-6 max-w-2xl text-[15px] leading-relaxed text-fg-muted sm:text-base"
          >
            {t("hero_sub")}
          </motion.p>

          <motion.div
            initial={reduce ? false : { opacity: 0, y: 16 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.55, delay: 0.19 }}
            className="mt-9 flex flex-col items-center justify-center gap-3 sm:flex-row"
          >
            <Button href="/run" size="lg">
              {t("hero_cta")}
              <ArrowRight className="h-4 w-4" />
            </Button>
            <Button href="/#cases" variant="outline" size="lg">
              {t("hero_cta2")}
            </Button>
          </motion.div>
        </div>

        {/* stats */}
        <motion.div
          initial={reduce ? false : { opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.6, delay: 0.28 }}
          className="mx-auto mt-16 grid max-w-3xl grid-cols-2 gap-3 sm:grid-cols-4"
        >
          {stats.map((s) => (
            <Card key={s.l} className="px-4 py-5 text-center">
              <div className="tabular font-display text-3xl font-semibold text-primary">
                {s.v}
              </div>
              <div className="mt-1 text-xs text-fg-muted">{s.l}</div>
            </Card>
          ))}
        </motion.div>
      </section>

      {/* ---------- HOW IT WORKS ---------- */}
      <section id="how" className="scroll-mt-20 px-5 py-20">
        <div className="mx-auto max-w-6xl">
          <SectionHeading
            eyebrow="Pipeline"
            title={t("how_title")}
            sub={t("how_sub")}
          />
          <div className="mt-12">
            <PipelineDiagram />
          </div>
        </div>
      </section>

      {/* ---------- TOOLS ---------- */}
      <section id="tools" className="scroll-mt-20 px-5 py-20">
        <div className="mx-auto max-w-6xl">
          <SectionHeading eyebrow="Stack" title={t("tools_title")} sub={t("tools_sub")} />
          <div className="mt-12">
            {config ? (
              <ToolsGrid steps={config.steps ?? []} />
            ) : (
              <div className="grid grid-cols-2 gap-3 sm:grid-cols-4">
                {Array.from({ length: 8 }).map((_, i) => (
                  <div key={i} className="h-16 animate-pulse rounded-xl bg-surface/40" />
                ))}
              </div>
            )}
          </div>
        </div>
      </section>

      {/* ---------- CASES ---------- */}
      <section id="cases" className="scroll-mt-20 px-5 py-20">
        <div className="mx-auto max-w-6xl">
          <SectionHeading eyebrow="Evidence" title={t("cases_title")} sub={t("cases_sub")} />
          <div className="mt-12">
            {runs.length ? (
              <CaseStudies runs={runs} onDeleted={(id) => setRuns((prev) => prev.filter((r) => r.id !== id))} />
            ) : (
              <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                {Array.from({ length: 3 }).map((_, i) => (
                  <div key={i} className="h-48 animate-pulse rounded-2xl bg-surface/40" />
                ))}
              </div>
            )}
          </div>
        </div>
      </section>

      {/* ---------- METHODOLOGY ---------- */}
      <section className="px-5 py-20">
        <div className="mx-auto max-w-6xl">
          <SectionHeading eyebrow="Trust" title={t("method_title")} />
          <div className="mt-12 grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">
            {[
              { icon: ShieldCheck, tr: ["Dürüst etiketleme", "Gerçek araç mı, yedek yöntem mi — her metrik açıkça işaretlenir."], en: ["Honest labeling", "Real tool vs. fallback — every metric is explicitly marked."] },
              { icon: Braces, tr: ["Tekrarlanabilir", "Koşulan tüm eşikler ve sürümler rapora yazılır."], en: ["Reproducible", "Every threshold and version used is written into the report."] },
              { icon: FlaskConical, tr: ["Deneysel bağ", "Adaylar IEDB'deki bilinen epitoplarla eşleştirilir."], en: ["Experimental link", "Candidates are matched against known IEDB epitopes."] },
              { icon: BookOpenCheck, tr: ["Tam atıf", "Her araç ve yöntem literatür atfıyla listelenir."], en: ["Fully cited", "Every tool and method is listed with its citation."] },
            ].map((f, i) => {
              const Icon = f.icon;
              const txt = t("nav_run") === "Run Pipeline" ? f.en : f.tr;
              return (
                <Card key={i} className="p-5">
                  <span className="grid h-10 w-10 place-items-center rounded-xl bg-bio/10 text-bio-soft ring-1 ring-bio/20">
                    <Icon className="h-5 w-5" />
                  </span>
                  <h3 className="mt-4 font-display text-base font-semibold text-fg">{txt[0]}</h3>
                  <p className="mt-1.5 text-[13px] leading-relaxed text-fg-muted">{txt[1]}</p>
                </Card>
              );
            })}
          </div>
        </div>
      </section>

      {/* ---------- CTA ---------- */}
      <section className="px-5 py-20">
        <div className="mx-auto max-w-4xl">
          <Card glow className="relative overflow-hidden px-8 py-14 text-center">
            <div className="grid-texture pointer-events-none absolute inset-0 -z-10 opacity-60" />
            <h2 className="font-display text-3xl font-semibold tracking-tight text-fg sm:text-4xl">
              {t("hero_title_1")} <span className="text-gradient">{t("hero_title_2")}</span>
            </h2>
            <div className="mt-8">
              <Button href="/run" size="lg">
                {t("hero_cta")}
                <ArrowRight className="h-4 w-4" />
              </Button>
            </div>
          </Card>
        </div>
      </section>

      {/* ---------- FOOTER ---------- */}
      <footer className="border-t border-line px-5 py-10">
        <div className="mx-auto flex max-w-6xl flex-col items-center justify-between gap-4 sm:flex-row">
          <div className="flex items-center gap-2 font-display font-semibold text-fg">
            Vax<span className="text-primary">Forge</span>
            <span className="ml-2 font-sans text-xs font-normal text-fg-faint">
              in silico reverse vaccinology
            </span>
          </div>
          <p className="max-w-md text-center text-xs text-fg-faint sm:text-right">
            {t("footer_note")}
          </p>
        </div>
      </footer>
    </main>
  );
}
