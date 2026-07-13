"use client";

import { useCallback, useEffect, useRef, useState } from "react";
import Link from "next/link";
import {
  UploadCloud,
  Play,
  Loader2,
  CheckCircle2,
  XCircle,
  ChevronRight,
  FileSearch,
  ArrowLeft,
} from "lucide-react";
import { Navbar } from "@/components/navbar";
import { Button, Card, Badge, SectionHeading } from "@/components/ui";
import { FunnelChart, CandidatesTable, DownloadBar } from "@/components/results";
import { useLang } from "@/components/lang-provider";
import {
  getConfig,
  getRun,
  runPipeline,
  AppConfig,
  SSEEvent,
  RunDetail,
  fileUrl,
} from "@/lib/api";
import { cn } from "@/lib/utils";

type Phase = "idle" | "running" | "done" | "error";

export default function RunPage() {
  const { t, lang } = useLang();
  const [config, setConfig] = useState<AppConfig | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [profile, setProfile] = useState("bacteria");
  const [gram, setGram] = useState("negative");
  const [hosts, setHosts] = useState<string[]>([]);
  const [dragOver, setDragOver] = useState(false);

  const [phase, setPhase] = useState<Phase>("idle");
  const [events, setEvents] = useState<SSEEvent[]>([]);
  const [progress, setProgress] = useState(0);
  const [errorMsg, setErrorMsg] = useState("");
  const [runId, setRunId] = useState<string | null>(null);
  const [detail, setDetail] = useState<RunDetail | null>(null);
  const consoleRef = useRef<HTMLDivElement>(null);
  const totalRef = useRef(7);
  const abortRef = useRef<AbortController | null>(null);

  useEffect(() => {
    getConfig()
      .then((c) => {
        setConfig(c);
        setProfile(c.default_profile);
        setHosts(c.default_hosts);
      })
      .catch(() => {});
  }, []);

  useEffect(() => {
    consoleRef.current?.scrollTo({ top: consoleRef.current.scrollHeight, behavior: "smooth" });
  }, [events]);

  const onFile = useCallback((f: File | null) => {
    if (!f) return;
    setFile(f);
    setPhase("idle");
    setEvents([]);
    setDetail(null);
    setRunId(null);
  }, []);

  const toggleHost = (name: string) =>
    setHosts((prev) => (prev.includes(name) ? prev.filter((h) => h !== name) : [...prev, name]));

  const start = async () => {
    if (!file) return;
    setPhase("running");
    setEvents([]);
    setProgress(0);
    setErrorMsg("");
    setDetail(null);
    setRunId(null);
    let done = 0;

    const form = new FormData();
    form.append("file", file);
    form.append("profile", profile);
    form.append("hosts", hosts.join(","));
    form.append("gram", profile === "bacteria" ? gram : "");
    form.append("lang", lang);

    const ac = new AbortController();
    abortRef.current = ac;

    try {
      await runPipeline(
        form,
        (ev) => {
          if (ev.phase === "__plan__") {
            totalRef.current = ev.data?.total ?? 7;
          } else if (ev.phase === "__done__") {
            setProgress(1);
            const id = ev.data?.run_id;
            if (id) {
              setRunId(id);
              getRun(id).then(setDetail).catch(() => {});
            }
            setPhase("done");
            return;
          } else if (ev.phase === "__error__") {
            setErrorMsg(ev.msg);
            setPhase("error");
            return;
          }
          setEvents((prev) => [...prev, ev]);
          if (ev.status === "done" || ev.status === "deferred") {
            done += 1;
            setProgress(Math.min(0.98, done / totalRef.current));
          }
        },
        ac.signal,
      );
    } catch (e: any) {
      if (e?.name !== "AbortError") {
        setErrorMsg(e?.message ?? "bağlantı hatası");
        setPhase("error");
      }
    }
  };

  return (
    <main className="relative min-h-dvh">
      <Navbar />
      <div className="mx-auto max-w-6xl px-5 pt-28 pb-20">
        <Link
          href="/"
          className="mb-6 inline-flex items-center gap-1.5 text-sm text-fg-muted transition hover:text-fg"
        >
          <ArrowLeft className="h-4 w-4" /> {lang === "tr" ? "Ana sayfa" : "Home"}
        </Link>

        <div className="grid grid-cols-1 gap-8 lg:grid-cols-[380px_1fr]">
          {/* ---- Sol: kontroller ---- */}
          <div className="space-y-5">
            <div>
              <h1 className="font-display text-2xl font-semibold text-fg">{t("run_title")}</h1>
            </div>

            {/* dosya yükleme */}
            <label
              onDragOver={(e) => {
                e.preventDefault();
                setDragOver(true);
              }}
              onDragLeave={() => setDragOver(false)}
              onDrop={(e) => {
                e.preventDefault();
                setDragOver(false);
                onFile(e.dataTransfer.files?.[0] ?? null);
              }}
              className={cn(
                "flex cursor-pointer flex-col items-center justify-center rounded-2xl border border-dashed p-8 text-center transition-colors",
                dragOver ? "border-primary bg-primary/5" : "border-line bg-surface/40 hover:border-primary/40",
              )}
            >
              <input
                type="file"
                className="hidden"
                accept=".fasta,.fa,.faa,.fna,.fastq,.fq,.gz,.gb,.gbk,.genbank,.gbff"
                onChange={(e) => onFile(e.target.files?.[0] ?? null)}
              />
              <UploadCloud className={cn("h-9 w-9", file ? "text-bio-soft" : "text-primary")} />
              <div className="mt-3 text-sm font-medium text-fg">
                {file ? file.name : t("run_drop")}
              </div>
              <div className="mt-1 text-[11px] leading-relaxed text-fg-faint">{t("run_formats")}</div>
            </label>

            {/* profil */}
            <Card className="p-4">
              <div className="text-xs font-medium uppercase tracking-wide text-fg-faint">
                {t("run_profile")}
              </div>
              <div className="mt-2 flex flex-wrap gap-1.5">
                {config?.profiles.map((p) => (
                  <button
                    key={p}
                    onClick={() => setProfile(p)}
                    className={cn(
                      "rounded-lg border px-3 py-1.5 text-xs font-medium capitalize transition cursor-pointer",
                      profile === p ? "border-primary/50 bg-primary/10 text-primary" : "border-line text-fg-muted hover:text-fg",
                    )}
                  >
                    {p}
                  </button>
                ))}
              </div>

              {profile === "bacteria" && (
                <div className="mt-4">
                  <div className="text-xs font-medium uppercase tracking-wide text-fg-faint">
                    {t("run_gram")}
                  </div>
                  <div className="mt-2 flex gap-1.5">
                    {[
                      { v: "negative", tr: "Gram−", en: "Gram−" },
                      { v: "positive", tr: "Gram+", en: "Gram+" },
                    ].map((g) => (
                      <button
                        key={g.v}
                        onClick={() => setGram(g.v)}
                        className={cn(
                          "rounded-lg border px-3 py-1.5 text-xs font-medium transition cursor-pointer",
                          gram === g.v ? "border-primary/50 bg-primary/10 text-primary" : "border-line text-fg-muted hover:text-fg",
                        )}
                      >
                        {lang === "tr" ? g.tr : g.en}
                      </button>
                    ))}
                  </div>
                </div>
              )}
            </Card>

            {/* konaklar */}
            <Card className="p-4">
              <div className="text-xs font-medium uppercase tracking-wide text-fg-faint">
                {t("run_hosts")}
              </div>
              <div className="mt-2 space-y-1.5">
                {config?.hosts.map((h) => (
                  <button
                    key={h.name}
                    onClick={() => toggleHost(h.name)}
                    className={cn(
                      "flex w-full items-center justify-between rounded-lg border px-3 py-2 text-left transition cursor-pointer",
                      hosts.includes(h.name) ? "border-primary/40 bg-primary/[0.07]" : "border-line hover:border-line",
                    )}
                  >
                    <span className="text-[13px] text-fg">{h.label}</span>
                    <span className="font-mono text-[11px] text-fg-faint">
                      I:{h.n_mhc_i} · II:{h.n_mhc_ii}
                    </span>
                  </button>
                ))}
              </div>
            </Card>

            <Button
              onClick={start}
              disabled={!file || phase === "running"}
              size="lg"
              className="w-full"
            >
              {phase === "running" ? (
                <>
                  <Loader2 className="h-4 w-4 animate-spin" /> {t("run_running")}
                </>
              ) : (
                <>
                  <Play className="h-4 w-4" /> {t("run_start")}
                </>
              )}
            </Button>
          </div>

          {/* ---- Sağ: konsol + sonuçlar ---- */}
          <div className="space-y-6">
            {phase === "idle" && !detail && (
              <Card className="flex min-h-[300px] flex-col items-center justify-center p-10 text-center">
                <FileSearch className="h-10 w-10 text-fg-faint" />
                <p className="mt-4 max-w-sm text-sm text-fg-muted">{t("run_pick")}</p>
              </Card>
            )}

            {(phase === "running" || phase === "error" || (phase === "done" && events.length > 0)) && (
              <Card className="overflow-hidden">
                <div className="flex items-center justify-between border-b border-line px-4 py-3">
                  <div className="flex items-center gap-2 text-sm font-medium text-fg">
                    {phase === "running" && <Loader2 className="h-4 w-4 animate-spin text-primary" />}
                    {phase === "done" && <CheckCircle2 className="h-4 w-4 text-bio-soft" />}
                    {phase === "error" && <XCircle className="h-4 w-4 text-danger" />}
                    {t("run_console")}
                  </div>
                  <span className="tabular font-mono text-xs text-fg-faint">
                    {Math.round(progress * 100)}%
                  </span>
                </div>
                <div className="h-1 w-full bg-white/[0.04]">
                  <div
                    className="h-full bg-gradient-to-r from-primary-deep to-primary transition-all duration-500"
                    style={{ width: `${progress * 100}%` }}
                  />
                </div>
                <div ref={consoleRef} className="max-h-80 space-y-1 overflow-y-auto p-4">
                  {events.map((ev, i) => (
                    <div key={i} className="flex items-start gap-2 text-[13px]">
                      <span className="mt-0.5 shrink-0">
                        {ev.status === "done" ? (
                          <CheckCircle2 className="h-3.5 w-3.5 text-bio-soft" />
                        ) : ev.status === "running" ? (
                          <Loader2 className="h-3.5 w-3.5 animate-spin text-primary" />
                        ) : (
                          <ChevronRight className="h-3.5 w-3.5 text-fg-faint" />
                        )}
                      </span>
                      <span className="font-mono text-[11px] text-fg-faint">{ev.phase}</span>
                      <span className="text-fg-muted">{ev.msg}</span>
                    </div>
                  ))}
                  {phase === "error" && (
                    <div className="mt-2 rounded-lg border border-danger/30 bg-danger/10 px-3 py-2 text-[13px] text-danger">
                      {errorMsg}
                    </div>
                  )}
                </div>
              </Card>
            )}

            {/* sonuçlar */}
            {detail && runId && (
              <>
                <div>
                  <SectionHeading center={false} title={t("res_flow")} />
                  <div className="mt-5">
                    <FunnelChart meta={detail} />
                  </div>
                </div>

                <div>
                  <div className="mb-4 flex items-center justify-between">
                    <h3 className="font-display text-xl font-semibold text-fg">{t("res_candidates")}</h3>
                    <DownloadBar runId={runId} meta={detail} />
                  </div>
                  {detail.candidates && detail.candidates.length > 0 && (
                    <CandidatesTable candidates={detail.candidates} />
                  )}
                </div>

                {detail.candidates && (
                  <div>
                    <h3 className="mb-4 font-display text-xl font-semibold text-fg">{t("res_report")}</h3>
                    <Card className="overflow-hidden">
                      <iframe
                        src={fileUrl(runId, "report.html")}
                        className="h-[600px] w-full bg-white"
                        title="report"
                      />
                    </Card>
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      </div>
    </main>
  );
}
