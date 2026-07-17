"use client";

import { useState } from "react";
import Link from "next/link";
import { ArrowUpRight, FileText, Dna, Bug, Microscope, Trash2, Loader2 } from "lucide-react";
import { RunSummary, deleteRun } from "@/lib/api";
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

export function CaseStudies({ runs, onDeleted }: { runs: RunSummary[]; onDeleted?: (id: string) => void }) {
  const { t } = useLang();
  const [deleting, setDeleting] = useState<string | null>(null);
  const top = runs.slice(0, 6);

  const onDelete = async (e: React.MouseEvent, id: string) => {
    e.preventDefault();
    e.stopPropagation();
    if (!window.confirm(t("runs_delete_confirm"))) return;
    setDeleting(id);
    try {
      await deleteRun(id);
      onDeleted?.(id);
    } catch {
      /* sessizce yut */
    } finally {
      setDeleting(null);
    }
  };

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
              <div className="flex items-center gap-1">
                <button
                  onClick={(e) => onDelete(e, r.id)}
                  disabled={deleting === r.id}
                  aria-label={t("runs_delete")}
                  title={t("runs_delete")}
                  className="grid h-8 w-8 place-items-center rounded-lg border border-line text-fg-faint transition hover:border-danger/50 hover:text-danger disabled:opacity-40 cursor-pointer"
                >
                  {deleting === r.id ? (
                    <Loader2 className="h-4 w-4 animate-spin" />
                  ) : (
                    <Trash2 className="h-4 w-4" />
                  )}
                </button>
                <ArrowUpRight className="h-5 w-5 text-fg-faint transition-colors group-hover:text-primary" />
              </div>
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
