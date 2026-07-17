"use client";

import { useMemo, useState } from "react";
import { Download, FileText, FileSpreadsheet, FileCode2, FileDown, Search, Syringe, Globe2, ShieldCheck } from "lucide-react";
import { RunDetail, Candidate, MevData, fileUrl } from "@/lib/api";
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
                <th className="hidden px-4 py-3 font-medium lg:table-cell">{t("col_allele")}</th>
                <th className="hidden px-4 py-3 text-right font-medium md:table-cell">{t("col_immuno")}</th>
                <th className="hidden px-4 py-3 text-right font-medium lg:table-cell">{t("col_cover")}</th>
                <th className="hidden px-4 py-3 text-center font-medium md:table-cell">{t("col_tox")}</th>
                <th className="hidden px-4 py-3 text-center font-medium md:table-cell">{t("col_allergen")}</th>
                <th className="px-4 py-3 text-right font-medium">{t("col_score")}</th>
              </tr>
            </thead>
            <tbody>
              {filtered.slice(0, limit).map((r, i) => {
                const st = kindStyle(String(r.kind));
                const tox = typeof r.toxicity === "number" ? r.toxicity : undefined;
                const allergen = String(r.allergen) === "true" || r.allergen === 1;
                const cover = typeof r.host_coverage === "number" ? r.host_coverage : undefined;
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
                    <td className="hidden px-4 py-2.5 font-mono text-[12px] text-fg-muted lg:table-cell">
                      {String(r.best_allele ?? "—")}
                    </td>
                    <td className="hidden px-4 py-2.5 text-right tabular font-mono text-[13px] text-fg-muted md:table-cell">
                      {typeof r.immunogenicity === "number" ? r.immunogenicity.toFixed(2) : "—"}
                    </td>
                    <td className="hidden px-4 py-2.5 text-right tabular font-mono text-[13px] text-fg-muted lg:table-cell">
                      {cover !== undefined ? `${(cover * 100).toFixed(0)}%` : "—"}
                    </td>
                    <td className="hidden px-4 py-2.5 text-center md:table-cell">
                      {tox !== undefined ? (
                        <span className={cn(
                          "inline-block rounded-md px-1.5 py-0.5 font-mono text-[11px] font-medium",
                          tox >= 0.5 ? "bg-danger/10 text-danger" : "bg-bio/10 text-bio-soft",
                        )}>
                          {tox.toFixed(2)}
                        </span>
                      ) : "—"}
                    </td>
                    <td className="hidden px-4 py-2.5 text-center md:table-cell">
                      <span className={cn(
                        "inline-block rounded-md px-1.5 py-0.5 text-[11px] font-medium",
                        allergen ? "bg-warn/10 text-warn" : "bg-bio/10 text-bio-soft",
                      )}>
                        {allergen ? (lang === "tr" ? "evet" : "yes") : (lang === "tr" ? "hayır" : "no")}
                      </span>
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

// ---------- MEV konstrukt görselleştirme ----------
const ROLE_STYLE: Record<string, { cls: string; ring: string }> = {
  adjuvant: { cls: "bg-primary/15 text-primary", ring: "ring-primary/40" },
  linker: { cls: "bg-white/[0.06] text-fg-faint", ring: "ring-line" },
  "epitope:B": { cls: "bg-epi-b/15 text-epi-b", ring: "ring-epi-b/40" },
  "epitope:MHC-I": { cls: "bg-epi-i/15 text-epi-i", ring: "ring-epi-i/40" },
  "epitope:MHC-II": { cls: "bg-epi-ii/15 text-epi-ii", ring: "ring-epi-ii/40" },
};
function compStyle(c: { role: string; kind?: string }) {
  if (c.role === "epitope") return ROLE_STYLE[`epitope:${c.kind}`] ?? ROLE_STYLE["epitope:MHC-I"];
  return ROLE_STYLE[c.role] ?? { cls: "bg-white/5 text-fg-muted", ring: "ring-line" };
}

function PropStat({ label, value, tone }: { label: string; value: React.ReactNode; tone?: "ok" | "warn" }) {
  return (
    <div className="rounded-xl border border-line bg-surface/40 px-3.5 py-3">
      <div className="text-[11px] uppercase tracking-wide text-fg-faint">{label}</div>
      <div className={cn(
        "mt-1 font-mono text-[15px] font-semibold",
        tone === "ok" ? "text-bio-soft" : tone === "warn" ? "text-warn" : "text-fg",
      )}>
        {value}
      </div>
    </div>
  );
}

export function MevConstruct({ mev }: { mev: MevData }) {
  const { t } = useLang();
  const p = mev.properties ?? {};
  const comps = mev.components ?? [];
  const nk = mev.n_by_kind ?? {};

  const antigen = p.antigenicity ?? {};
  const solub = p.solubility ?? {};
  const disorder = p.disorder ?? {};
  const sec = p.secondary_structure ?? {};

  return (
    <Card className="p-5 sm:p-6">
      <div className="flex items-center gap-2.5">
        <span className="grid h-9 w-9 place-items-center rounded-xl bg-primary/10 ring-1 ring-primary/25">
          <Syringe className="h-5 w-5 text-primary" />
        </span>
        <div>
          <h3 className="font-display text-lg font-semibold text-fg">{t("mev_title")}</h3>
          <p className="text-xs text-fg-muted">{t("mev_sub")}</p>
        </div>
      </div>

      {/* epitop sayıları + adjuvan */}
      <div className="mt-4 flex flex-wrap items-center gap-2">
        {mev.adjuvant && (
          <Badge tone="primary">{t("mev_adjuvant")}: {mev.adjuvant_label ?? mev.adjuvant}</Badge>
        )}
        {Object.entries(nk).map(([k, v]) => {
          const st = kindStyle(k);
          return (
            <span key={k} className={cn("inline-flex items-center gap-1.5 rounded-full border px-2.5 py-1 text-[11px] font-medium", st.cls)}>
              <span className={cn("h-1.5 w-1.5 rounded-full", st.dot)} />
              {k}: {v}
            </span>
          );
        })}
        {typeof p.length === "number" && (
          <Badge tone="muted">{p.length} aa</Badge>
        )}
      </div>

      {/* bileşen şeması (renkli parçalar) */}
      {comps.length > 0 && (
        <div className="mt-4">
          <div className="flex flex-wrap gap-1 rounded-xl border border-line bg-ink/40 p-2.5">
            {comps.map((c, i) => {
              const st = compStyle(c);
              const label = c.role === "adjuvant" ? (c.source ?? "adjuvant")
                : c.role === "linker" ? c.seq
                : c.kind ?? "epi";
              return (
                <span
                  key={i}
                  title={`${c.role}${c.kind ? " · " + c.kind : ""}${c.source ? " · " + c.source : ""}\n${c.seq}`}
                  className={cn(
                    "cursor-default rounded-md px-2 py-1 font-mono text-[10px] font-medium ring-1",
                    st.cls, st.ring,
                  )}
                >
                  {c.role === "linker" ? label : `${label}`}
                </span>
              );
            })}
          </div>
          {/* lejant */}
          <div className="mt-2 flex flex-wrap gap-x-4 gap-y-1 text-[11px] text-fg-faint">
            <span className="inline-flex items-center gap-1"><span className="h-2 w-2 rounded-full bg-primary" />{t("mev_adjuvant")}</span>
            <span className="inline-flex items-center gap-1"><span className="h-2 w-2 rounded-full bg-epi-b" />B</span>
            <span className="inline-flex items-center gap-1"><span className="h-2 w-2 rounded-full bg-epi-i" />MHC-I</span>
            <span className="inline-flex items-center gap-1"><span className="h-2 w-2 rounded-full bg-epi-ii" />MHC-II</span>
            <span className="inline-flex items-center gap-1"><span className="h-2 w-2 rounded-full bg-fg-faint" />{t("mev_linker")} (AAY/GPGPG)</span>
          </div>
        </div>
      )}

      {/* tam dizi */}
      <div className="mt-4">
        <div className="mb-1 text-[11px] uppercase tracking-wide text-fg-faint">{t("mev_seq")}</div>
        <div className="max-h-32 overflow-y-auto rounded-xl border border-line bg-ink/40 p-3 font-mono text-[12px] leading-relaxed tracking-wide text-fg-muted break-all">
          {mev.sequence}
        </div>
      </div>

      {/* fizikokimyasal özellikler */}
      <div className="mt-4 grid grid-cols-2 gap-2.5 sm:grid-cols-3 lg:grid-cols-4">
        {typeof p.length === "number" && <PropStat label={t("mev_len")} value={p.length} />}
        {typeof p.mw_kda === "number" && <PropStat label={t("mev_mw")} value={`${p.mw_kda} kDa`} />}
        {typeof p.pI === "number" && <PropStat label={t("mev_pi")} value={p.pI} />}
        {typeof p.instability === "number" && (
          <PropStat label={t("mev_instab")} value={`${p.instability} · ${p.stable ? t("mev_stable") : t("mev_unstable")}`} tone={p.stable ? "ok" : "warn"} />
        )}
        {typeof p.aliphatic_index === "number" && <PropStat label={t("mev_aliphatic")} value={p.aliphatic_index} />}
        {typeof p.gravy === "number" && <PropStat label={t("mev_gravy")} value={p.gravy} />}
        {typeof antigen.score === "number" && (
          <PropStat label={t("mev_antigen")} value={`${antigen.score} · ${antigen.antigenic ? t("mev_yes_antigen") : t("mev_no_antigen")}`} tone={antigen.antigenic ? "ok" : "warn"} />
        )}
        {typeof solub.percent_sol === "number" && (
          <PropStat label={t("mev_solub")} value={`${solub.percent_sol}% · ${solub.soluble ? t("mev_soluble") : t("mev_insoluble")}`} tone={solub.soluble ? "ok" : "warn"} />
        )}
        {typeof disorder.percent_disordered === "number" && (
          <PropStat label={t("mev_disorder")} value={`${disorder.percent_disordered}%`} />
        )}
        {typeof sec.helix_pct === "number" && (
          <PropStat label={t("mev_secstruct")} value={`${sec.helix_pct}/${sec.strand_pct}/${sec.coil_pct}%`} />
        )}
        {p.allergen && <PropStat label={t("col_allergen")} value={p.allergen.allergenic ? "⚠" : t("mev_nonallergen")} tone={p.allergen.allergenic ? "warn" : "ok"} />}
        {p.toxicity && <PropStat label={t("col_tox")} value={p.toxicity.toxic ? "⚠" : t("mev_nontoxic")} tone={p.toxicity.toxic ? "warn" : "ok"} />}
      </div>
    </Card>
  );
}

// ---------- Popülasyon kapsamı ----------
export function PopulationCoverage({ pc }: { pc: any }) {
  const { t } = useLang();
  if (!pc || pc.available === false) return null;
  const hosts = pc.hosts ?? {};
  const hostKeys = Object.keys(hosts);
  if (hostKeys.length === 0) return null;
  const areas: string[] = pc.areas ?? [];

  return (
    <Card className="p-5 sm:p-6">
      <div className="flex items-center gap-2.5">
        <span className="grid h-9 w-9 place-items-center rounded-xl bg-bio/10 ring-1 ring-bio/25">
          <Globe2 className="h-5 w-5 text-bio-soft" />
        </span>
        <div>
          <h3 className="font-display text-lg font-semibold text-fg">{t("pc_title")}</h3>
          <p className="text-xs text-fg-muted">{t("pc_sub")}</p>
        </div>
      </div>

      {hostKeys.map((hk) => {
        const h = hosts[hk];
        const cls: Record<string, any> = h.mhc_i ?? h.combined ?? h.mhc_ii ?? {};
        const rows = (areas.length ? areas : Object.keys(cls)).filter((a) => cls[a]);
        return (
          <div key={hk} className="mt-4">
            <div className="mb-2 text-[13px] font-medium text-fg">{h.label ?? hk}</div>
            <div className="space-y-1.5">
              {rows.map((area) => {
                const cov = cls[area]?.coverage ?? 0;
                return (
                  <div key={area} className="flex items-center gap-3">
                    <div className="w-28 shrink-0 text-right text-[12px] text-fg-muted sm:w-36">{area}</div>
                    <div className="relative h-6 flex-1 overflow-hidden rounded-md bg-white/[0.03]">
                      <div
                        className="flex h-full items-center rounded-md bg-gradient-to-r from-bio/60 to-bio-soft/70 px-2"
                        style={{ width: `${Math.max(4, cov)}%` }}
                      >
                        <span className="font-mono text-[11px] font-semibold text-ink-900">{cov.toFixed(1)}%</span>
                      </div>
                    </div>
                  </div>
                );
              })}
            </div>
          </div>
        );
      })}
    </Card>
  );
}

// ---------- IEDB doğrulama ----------
export function IedbValidation({ iedb }: { iedb: any }) {
  const { t } = useLang();
  if (!iedb || iedb.available === false) return null;
  const b = iedb.benchmark ?? {};
  const items = [
    { label: t("iedb_matched"), value: `${iedb.n_matched ?? "—"} / ${iedb.n_candidates ?? "—"}` },
    { label: t("iedb_recall"), value: typeof b.recall === "number" ? `${(b.recall * 100).toFixed(1)}%` : "—" },
    { label: t("iedb_precision"), value: typeof b.precision_like === "number" ? `${(b.precision_like * 100).toFixed(1)}%` : "—" },
    { label: t("iedb_known"), value: b.n_known ?? iedb.n_records ?? "—" },
  ];
  return (
    <Card className="p-5 sm:p-6">
      <div className="flex items-center gap-2.5">
        <span className="grid h-9 w-9 place-items-center rounded-xl bg-epi-i/10 ring-1 ring-epi-i/25">
          <ShieldCheck className="h-5 w-5 text-epi-i" />
        </span>
        <div>
          <h3 className="font-display text-lg font-semibold text-fg">{t("iedb_title")}</h3>
          <p className="text-xs text-fg-muted">{t("iedb_sub")}</p>
        </div>
      </div>
      <div className="mt-4 grid grid-cols-2 gap-2.5 sm:grid-cols-4">
        {items.map((it) => (
          <PropStat key={it.label} label={it.label} value={it.value} />
        ))}
      </div>
      {iedb.source && <p className="mt-3 text-[11px] text-fg-faint">{iedb.source}</p>}
    </Card>
  );
}
