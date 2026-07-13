"use client";

import { CheckCircle2, AlertTriangle } from "lucide-react";
import { ToolStatus } from "@/lib/api";
import { useLang } from "./lang-provider";
import { cn } from "@/lib/utils";

export function ToolsGrid({ tools }: { tools: ToolStatus[] }) {
  const { t } = useLang();
  return (
    <div className="grid grid-cols-2 gap-3 sm:grid-cols-3 lg:grid-cols-4">
      {tools.map((tool) => (
        <div
          key={tool.module}
          className={cn(
            "flex items-center gap-3 rounded-xl border p-3.5 transition-colors",
            tool.available
              ? "border-bio/25 bg-bio/[0.06]"
              : "border-line bg-white/[0.02]",
          )}
        >
          {tool.available ? (
            <CheckCircle2 className="h-5 w-5 shrink-0 text-bio-soft" />
          ) : (
            <AlertTriangle className="h-5 w-5 shrink-0 text-warn" />
          )}
          <div className="min-w-0">
            <div className="truncate font-mono text-[13px] font-medium text-fg">
              {tool.label}
            </div>
            <div
              className={cn(
                "text-[11px] uppercase tracking-wide",
                tool.available ? "text-bio-soft/80" : "text-warn/80",
              )}
            >
              {tool.available ? t("tools_available") : t("tools_fallback")}
            </div>
          </div>
        </div>
      ))}
    </div>
  );
}
