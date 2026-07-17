"use client";

import { CheckCircle2, AlertTriangle } from "lucide-react";
import { PipelineStep } from "@/lib/api";
import { useLang } from "./lang-provider";
import { cn } from "@/lib/utils";

// Grid artık "kurulu/yedek" değil, KODUN GERÇEĞİNİ gösterir:
// her kart bir pipeline adımı → o adımda gerçekte çalışan yöntem.
export function ToolsGrid({ steps }: { steps: PipelineStep[] }) {
  const { lang } = useLang();
  return (
    <div className="grid grid-cols-2 gap-3 sm:grid-cols-3 lg:grid-cols-4">
      {steps.map((s, i) => {
        const stepLabel = lang === "tr" ? s.step_tr : s.step_en;
        // installed=false → o adımda yerel/heuristik yöntem çalışıyor (yine gerçek).
        const external = s.installed;
        return (
          <div
            key={s.id}
            className={cn(
              "flex items-start gap-3 rounded-xl border p-3.5 transition-colors",
              external ? "border-bio/25 bg-bio/[0.06]" : "border-line bg-white/[0.02]",
            )}
          >
            {external ? (
              <CheckCircle2 className="mt-0.5 h-5 w-5 shrink-0 text-bio-soft" />
            ) : (
              <AlertTriangle className="mt-0.5 h-5 w-5 shrink-0 text-warn" />
            )}
            <div className="min-w-0">
              <div className="flex items-center gap-1.5 text-[10px] uppercase tracking-wide text-fg-faint">
                <span className="tabular font-mono">{String(i + 1).padStart(2, "0")}</span>
                <span className="truncate">{stepLabel}</span>
              </div>
              <div className="mt-0.5 truncate font-mono text-[13px] font-medium text-fg" title={s.method}>
                {s.method}
              </div>
            </div>
          </div>
        );
      })}
    </div>
  );
}
