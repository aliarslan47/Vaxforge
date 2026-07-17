"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { ArrowLeft, Loader2, Calendar, FileText, FlaskConical, Plus, ChevronRight, Trash2 } from "lucide-react";
import { Navbar } from "@/components/navbar";
import { Card, Badge, Button } from "@/components/ui";
import { useLang } from "@/components/lang-provider";
import { getRuns, deleteRun, RunSummary } from "@/lib/api";

export default function RunsPage() {
  const { t, lang } = useLang();
  const [runs, setRuns] = useState<RunSummary[] | null>(null);
  const [error, setError] = useState(false);
  const [deleting, setDeleting] = useState<string | null>(null);

  useEffect(() => {
    getRuns().then(setRuns).catch(() => setError(true));
  }, []);

  const onDelete = async (e: React.MouseEvent, id: string) => {
    e.preventDefault();
    e.stopPropagation();
    if (!window.confirm(t("runs_delete_confirm"))) return;
    setDeleting(id);
    try {
      await deleteRun(id);
      setRuns((prev) => (prev ? prev.filter((r) => r.id !== id) : prev));
    } catch {
      /* sessizce yut — kart kalır */
    } finally {
      setDeleting(null);
    }
  };

  return (
    <main className="relative min-h-dvh">
      <Navbar />
      <div className="mx-auto max-w-5xl px-5 pt-28 pb-20">
        <Link
          href="/"
          className="mb-6 inline-flex items-center gap-1.5 text-sm text-fg-muted transition hover:text-fg"
        >
          <ArrowLeft className="h-4 w-4" /> {lang === "tr" ? "Ana sayfa" : "Home"}
        </Link>

        <div className="mb-8 flex flex-col gap-4 sm:flex-row sm:items-end sm:justify-between">
          <div>
            <h1 className="font-display text-3xl font-semibold text-fg">{t("runs_title")}</h1>
            <p className="mt-2 max-w-xl text-[15px] leading-relaxed text-fg-muted">{t("runs_sub")}</p>
          </div>
          <Button href="/run" className="shrink-0">
            <Plus className="h-4 w-4" /> {t("runs_new")}
          </Button>
        </div>

        {error && (
          <Card className="p-10 text-center text-fg-muted">
            {lang === "tr" ? "Koşular yüklenemedi." : "Could not load runs."}
          </Card>
        )}

        {!runs && !error && (
          <div className="flex min-h-[240px] items-center justify-center">
            <Loader2 className="h-6 w-6 animate-spin text-primary" />
          </div>
        )}

        {runs && runs.length === 0 && (
          <Card className="flex flex-col items-center justify-center gap-4 p-14 text-center">
            <FlaskConical className="h-10 w-10 text-fg-faint" />
            <p className="text-sm text-fg-muted">{t("runs_empty")}</p>
            <Button href="/run"><Plus className="h-4 w-4" /> {t("runs_new")}</Button>
          </Card>
        )}

        {runs && runs.length > 0 && (
          <div className="grid grid-cols-1 gap-3 sm:grid-cols-2">
            {runs.map((r) => (
              <Link key={r.id} href={`/results/${r.id}`} className="group">
                <Card className="h-full p-5 transition-colors group-hover:border-primary/40">
                  <div className="flex items-start justify-between gap-3">
                    <div className="min-w-0">
                      <div className="flex flex-wrap items-center gap-1.5">
                        <Badge tone="primary">{r.profile}</Badge>
                        {r.hosts?.slice(0, 3).map((h) => (
                          <Badge key={h} tone="muted">{h}</Badge>
                        ))}
                      </div>
                      <h3 className="mt-2.5 truncate font-display text-lg font-semibold text-fg">
                        {r.input || r.id}
                      </h3>
                    </div>
                    <div className="flex shrink-0 items-center gap-1">
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
                      <ChevronRight className="mt-0.5 h-5 w-5 text-fg-faint transition group-hover:translate-x-0.5 group-hover:text-primary" />
                    </div>
                  </div>
                  <div className="mt-3 flex flex-wrap items-center gap-x-4 gap-y-1 text-xs text-fg-faint">
                    <span className="inline-flex items-center gap-1.5">
                      <Calendar className="h-3.5 w-3.5" />
                      {String(r.timestamp ?? "").replace("T", " ").slice(0, 16)}
                    </span>
                    <span className="inline-flex items-center gap-1.5">
                      <FileText className="h-3.5 w-3.5" />
                      {r.n_candidates} {t("runs_cands").toLowerCase()}
                    </span>
                    <span className="inline-flex items-center gap-1.5 font-mono uppercase">
                      {r.molecule}
                    </span>
                  </div>
                  <div className="mt-3 flex gap-1.5">
                    {r.has_pdf && <span className="rounded border border-line px-1.5 py-0.5 text-[10px] font-medium text-fg-faint">PDF</span>}
                    {r.has_xlsx && <span className="rounded border border-line px-1.5 py-0.5 text-[10px] font-medium text-fg-faint">XLSX</span>}
                    {r.has_html && <span className="rounded border border-line px-1.5 py-0.5 text-[10px] font-medium text-fg-faint">HTML</span>}
                  </div>
                </Card>
              </Link>
            ))}
          </div>
        )}
      </div>
    </main>
  );
}
