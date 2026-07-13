"use client";

import { useMemo, useState } from "react";
import { Download, FileText, FileSpreadsheet, FileCode2, FileDown, Search } from "lucide-react";
import { RunDetail, Candidate, fileUrl } from "@/lib/api";
import { useLang } from "./lang-provider";
import { Card, Badge } from "./ui";
import { cn } from "@/lib/utils";

const KIND_STYLE: Record<string, { label: string; cls: string; dot: string }> = {
  B: { label: "B", cls: "text-epi-b bg-epi-b/10 border-epi-b/25", dot: "bg-epi-b" },
  "MHC-I": { label: "MHC-I", cls: "text-epi-i bg-epi-i/10 border-epi-i/25", dot: "bg-epi-i" },
  "MHC-II": { label: "MHC-II", cls: "text-epi-ii bg-epi-ii/10 border-epi-ii/25", dot: "bg-epi-ii" },
};

function kindStyle(kind: string) {
  return KIND_STYLE[kind] ?? { label: kind, cls: "text-fg-muted bg-white/5 border-line", dot: "bg-fg-muted" };
}

// ---------- Eleme hunisi ----------
export function FunnelChart({ meta }: { meta: RunDetail }) {
  const { lang } = useLang();
  const unit =
    meta.molecule === "cds" ? "CDS" : meta.molecule === "reads" ? (lang === "tr" ? "okuma" : "reads") : (lang === "tr" ? "protein" : "protein");

  const steps: { label: string; n: number | undefined }[] = [
    { label: lang === "tr" ? `Girdi (${unit})` : `Input (${unit})`, n: num(meta.n_raw) },
    { label: lang === "tr" ? "Proteine çevrildi" : "Translated" , n: num(meta.n_input) },
    { label: lang === "tr" ? "Virülans skoru" : "Virulence scored", n: num(meta.n_discovery) },
    { label: lang === "tr" ? "Huni sonrası" : "After funnel", n: num(meta.n_funnel) },
    { label: lang === "tr" ? "Alerjenite sonrası" : "After allergenicity", n: num(meta.n_after_allergen) },
    { label: lang === "tr" ? "Toksisite sonrası" : "After toxicity", n: num(meta.n_after_toxicity) },
    { label: lang === "tr" ? "Sıralanan aday" : "Ranked candidates", n: (meta.candidates?.length ?? num(meta.n_survivors)) },
  ].filter((s) => s.n !== undefined);

  const max = Math.max(...steps.map((s) => s.n ?? 0), 1);

  return (
    <div className="space-y-2.5">
      {steps.map((s, i) => (
        <div key={i} className="flex items-center gap-3">
          <div className="w-40 shrink-0 text-right text-[13px] text-fg-muted sm:w-52">{s.label}</div>
          <div className="relative h-8 flex-1 overflow-hidden rounded-lg bg-white/[0.03]">
            <div
              className="flex h-full items-center rounded-lg bg-gradient-to-r from-primary-deep/70 to-primary/70 px-3"
              style={{ width: `${Math.max(6, ((s.n ?? 0) / max) * 100)}%` }}
            >
              <span className="tabular font-mono text-[12px] font-semibold text-ink-900">
                {(s.n ?? 0).toLocaleString()}
              </span>
            </div>
          </div>
        </div>
      ))}
    </div>
  );
}

// ---------- Aday tablosu ----------
export function CandidatesTable({ candidates }: { candidates: Candidate[] }) {
  const { t, lang } = useLang();
  const [q, setQ] = useState("");
  const [kindFilter, setKindFilter] = useState<string>("");
  const [limit, setLimit] = useState(50);

  const filtered = useMemo(() => {
    let rows = candidates;
    if (kindFilter) rows = rows.filter((r) => String(r.kind) === kindFilter);
    if (q.trim()) {
      const s = q.trim().toLowerCase();
      rows = rows.filter(
        (r) =>
          String(r.peptide).toLowerCase().includes(s) ||
          String(r.gene ?? "").toLowerCase().includes(s) ||
          String(r.cds_source ?? "").toLowerCase().includes(s),
      );
    }
    return rows;
  }, [candidates, q, kindFilter]);

  const kinds = ["B", "MHC-I", "MHC-II"];

  return (
    <div>
      <div className="mb-4 flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
        <div className="relative w-full sm:max-w-xs">
          <Search className="pointer-events-none absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-fg-faint" />
          <input
            value={q}
            onChange={(e) => setQ(e.target.value)}
            placeholder={lang === "tr" ? "Peptit / gen ara…" : "Search peptide / gene…"}
            className="h-10 w-full rounded-lg border border-line bg-surface/60 pl-9 pr-3 text-sm text-fg placeholder:text-fg-faint focus:border-primary/50 focus:outline-none"
          />
        </div>
        <div className="flex items-center gap-1.5">
          <button
            onClick={() => setKindFilter("")}
            className={cn(
              "rounded-lg border px-3 py-1.5 text-xs font-medium transition cursor-pointer",
              kindFilter === "" ? "border-primary/50 bg-primary/10 text-primary" : "border-line text-fg-muted hover:text-fg",
            )}
          >
            {lang === "tr" ? "Tümü" : "All"}
          </button>
          {kinds.map((k) => {
            const st = kindStyle(k);
            return (
              <button
                key={k}
                onClick={() => setKindFilter(kindFilter === k ? "" : k)}
                className={cn(
                  "inline-flex items-center gap-1.5 rounded-lg border px-3 py-1.5 text-xs font-medium transition cursor-pointer",
                  kindFilter === k ? st.cls : "border-line text-fg-muted hover:text-fg",
                )}
              >
                <span className={cn("h-2 w-2 rounded-full", st.dot)} />
                {k}
              </button>
            );
          })}
        </div>
      </div>

      <Card className="overflow-hidden">
        <div className="overflow-x-auto">
          <table className="w-full border-collapse text-sm">
            <thead>
              <tr className="border-b border-line text-left text-xs uppercase tracking-wide text-fg-faint">
                <th className="px-4 py-3 font-medium">{t("col_rank")}</th>
                <th className="px-4 py-3 font-medium">{t("col_seq")}</th>
                <th className="px-4 py-3 font-medium">{t("col_type")}</th>
                <th className="px-4 py-3 font-medium">{t("col_parent")}</th>
                <th className="px-4 py-3 text-right font-medium">{t("col_score")}</th>
              </tr>
            </thead>
            <tbody>
              {filtered.slice(0, limit).map((r, i) => {
                const st = kindStyle(String(r.kind));
                return (
                  <tr
                    key={i}
                    className="border-b border-line/60 transition-colors last:border-0 hover:bg-white/[0.02]"
                  >
                    <td className="px-4 py-2.5 tabular font-mono text-xs text-fg-faint">{String(r.rank ?? i + 1)}</td>
                    <td className="px-4 py-2.5 font-mono text-[13px] text-fg">{String(r.peptide)}</td>
                    <td className="px-4 py-2.5">
                      <span className={cn("inline-flex items-center gap-1.5 rounded-md border px-2 py-0.5 text-[11px] font-medium", st.cls)}>
                        <span className={cn("h-1.5 w-1.5 rounded-full", st.dot)} />
                        {st.label}
                      </span>
                    </td>
                    <td className="px-4 py-2.5 text-[13px] text-fg-muted">
                      {String(r.gene || r.locus_tag || r.cds_source || "—")}
                    </td>
                    <td className="px-4 py-2.5 text-right tabular font-mono text-[13px] font-semibold text-primary">
                      {typeof r.candidacy === "number" ? r.candidacy.toFixed(3) : String(r.candidacy ?? "")}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
        {filtered.length > limit && (
          <div className="border-t border-line p-3 text-center">
            <button
              onClick={() => setLimit((l) => l + 50)}
              className="rounded-lg border border-line px-4 py-1.5 text-xs text-fg-muted transition hover:text-fg cursor-pointer"
            >
              {lang === "tr" ? `Daha fazla (${filtered.length - limit} kaldı)` : `Show more (${filtered.length - limit} left)`}
            </button>
          </div>
        )}
      </Card>
      <div className="mt-2 text-xs text-fg-faint">
        {filtered.length.toLocaleString()} / {candidates.length.toLocaleString()} {lang === "tr" ? "aday" : "candidates"}
      </div>
    </div>
  );
}

// ---------- İndirme çubuğu ----------
export function DownloadBar({ runId, meta }: { runId: string; meta: RunDetail }) {
  const { lang } = useLang();
  const items = [
    { name: "report.pdf", label: "PDF", icon: FileText, show: true },
    { name: "candidates.csv", label: "CSV", icon: FileDown, show: true },
    { name: "candidates_full.xlsx", label: "Excel", icon: FileSpreadsheet, show: true },
    { name: "top_peptides.fasta", label: "FASTA", icon: FileCode2, show: true },
  ];
  return (
    <div className="flex flex-wrap gap-2">
      {items.map((it) => {
        const Icon = it.icon;
        return (
          <a
            key={it.name}
            href={fileUrl(runId, it.name)}
            target="_blank"
            rel="noopener"
            className="inline-flex items-center gap-2 rounded-lg border border-line bg-surface/50 px-3.5 py-2 text-sm text-fg-muted transition hover:border-primary/40 hover:text-fg cursor-pointer"
          >
            <Icon className="h-4 w-4" />
            {it.label}
            <Download className="h-3.5 w-3.5 opacity-50" />
          </a>
        );
      })}
    </div>
  );
}

function num(v: unknown): number | undefined {
  return typeof v === "number" ? v : undefined;
}
