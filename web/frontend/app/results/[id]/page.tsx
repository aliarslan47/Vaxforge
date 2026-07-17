"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { useParams } from "next/navigation";
import { ArrowLeft, Loader2, Calendar, FileText } from "lucide-react";
import { Navbar } from "@/components/navbar";
import { Card, Badge, SectionHeading } from "@/components/ui";
import { FunnelChart, CandidatesTable, DownloadBar, MevConstruct, PopulationCoverage, IedbValidation } from "@/components/results";
import { useLang } from "@/components/lang-provider";
import { getRun, RunDetail } from "@/lib/api";

export default function ResultsPage() {
  const { t, lang } = useLang();
  const params = useParams();
  const id = String(params.id);
  const [detail, setDetail] = useState<RunDetail | null>(null);
  const [error, setError] = useState(false);

  useEffect(() => {
    getRun(id).then(setDetail).catch(() => setError(true));
  }, [id]);

  return (
    <main className="relative min-h-dvh">
      <Navbar />
      <div className="mx-auto max-w-6xl px-5 pt-28 pb-20">
        <Link
          href="/#cases"
          className="mb-6 inline-flex items-center gap-1.5 text-sm text-fg-muted transition hover:text-fg"
        >
          <ArrowLeft className="h-4 w-4" /> {t("cases_title")}
        </Link>

        {error && (
          <Card className="p-10 text-center text-fg-muted">
            {lang === "tr" ? "Koşu bulunamadı." : "Run not found."}
          </Card>
        )}

        {!detail && !error && (
          <div className="flex min-h-[300px] items-center justify-center">
            <Loader2 className="h-6 w-6 animate-spin text-primary" />
          </div>
        )}

        {detail && (
          <div className="space-y-8">
            <div>
              <div className="flex flex-wrap items-center gap-2">
                <Badge tone="primary">{String(detail.profile ?? "")}</Badge>
                {(detail as any).hosts?.map((h: any) => (
                  <Badge key={h.name ?? h} tone="muted">
                    {h.label ?? h.name ?? String(h)}
                  </Badge>
                ))}
              </div>
              <h1 className="mt-3 font-display text-3xl font-semibold text-fg">
                {String(detail.input ?? id)}
              </h1>
              <div className="mt-2 flex items-center gap-4 text-xs text-fg-faint">
                <span className="inline-flex items-center gap-1.5">
                  <Calendar className="h-3.5 w-3.5" />
                  {String(detail.timestamp ?? "").replace("T", " ")}
                </span>
                <span className="inline-flex items-center gap-1.5">
                  <FileText className="h-3.5 w-3.5" />
                  {detail.candidates?.length ?? 0} {lang === "tr" ? "aday" : "candidates"}
                </span>
              </div>
            </div>

            <DownloadBar runId={id} meta={detail} />

            <div>
              <SectionHeading center={false} title={t("res_flow")} />
              <div className="mt-5">
                <FunnelChart meta={detail} />
              </div>
            </div>

            {detail.candidates && detail.candidates.length > 0 && (
              <div>
                <h3 className="mb-4 font-display text-xl font-semibold text-fg">
                  {t("res_candidates")}
                </h3>
                <CandidatesTable candidates={detail.candidates} />
              </div>
            )}

            {detail.mev && <MevConstruct mev={detail.mev} />}
            {(detail as any).population_coverage && <PopulationCoverage pc={(detail as any).population_coverage} />}
            {(detail as any).iedb_match && <IedbValidation iedb={(detail as any).iedb_match} />}
          </div>
        )}
      </div>
    </main>
  );
}
