"use client";

import Link from "next/link";
import { ArrowUpRight, FileText, Dna, Bug, Microscope } from "lucide-react";
import { RunSummary } from "@/lib/api";
import { useLang } from "./lang-provider";
import { Badge } from "./ui";

// Organizma adına göre ikon/etiket (kabaca).
function organismIcon(input: string) {
  const s = input.toLowerCase();
  if (s.includes("ebola") || s.includes("wuhan") || s.includes("sarscov") || s.includes("virus"))
    return Bug;
  if (s.includes("menb") || s.includes("neisseria") || s.includes("br1")) return Microscope;
  return Dna;
}

export function CaseStudies({ runs }: { runs: RunSummary[] }) {
  const { t } = useLang();
  const top = runs.slice(0, 6);

  return (
    <div className="grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-3">
      {top.map((r) => {
        const Icon = organismIcon(r.input);
        return (
          <Link
            key={r.id}
            href={`/results/${r.id}`}
            className="group flex flex-col rounded-2xl border border-line bg-surface/50 p-5 backdrop-blur-sm transition-all hover:border-primary/40 hover:bg-surface-hover/40"
          >
            <div className="flex items-start justify-between">
              <span className="grid h-10 w-10 place-items-center rounded-xl bg-primary/10 text-primary ring-1 ring-primary/20">
                <Icon className="h-5 w-5" />
              </span>
              <ArrowUpRight className="h-5 w-5 text-fg-faint transition-colors group-hover:text-primary" />
            </div>
            <h3 className="mt-4 truncate font-display text-lg font-semibold text-fg">
              {r.input}
            </h3>
            <div className="mt-1 flex flex-wrap items-center gap-2">
              <Badge tone="muted">{r.profile}</Badge>
              {r.hosts.slice(0, 2).map((h) => (
                <Badge key={h} tone="primary">
                  {h}
                </Badge>
              ))}
            </div>
            <div className="mt-4 flex items-end justify-between border-t border-line pt-4">
              <div>
                <div className="tabular font-display text-2xl font-semibold text-fg">
                  {r.n_candidates.toLocaleString()}
                </div>
                <div className="text-xs text-fg-muted">{t("cases_candidates")}</div>
              </div>
              {r.has_pdf && (
                <span className="inline-flex items-center gap-1 text-xs text-fg-faint">
                  <FileText className="h-3.5 w-3.5" /> PDF
                </span>
              )}
            </div>
          </Link>
        );
      })}
    </div>
  );
}
