"use client";

import Link from "next/link";
import { cn } from "@/lib/utils";

// Hafif, on-brand UI primitifleri (shadcn yerine — Node 18 kurulum sadeliği).

export function Button({
  children,
  variant = "primary",
  size = "md",
  className,
  href,
  ...props
}: {
  children: React.ReactNode;
  variant?: "primary" | "ghost" | "outline";
  size?: "md" | "lg";
  className?: string;
  href?: string;
} & React.ButtonHTMLAttributes<HTMLButtonElement>) {
  const base =
    "inline-flex items-center justify-center gap-2 rounded-xl font-medium transition-all duration-200 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary/70 focus-visible:ring-offset-2 focus-visible:ring-offset-ink disabled:opacity-50 disabled:pointer-events-none cursor-pointer active:scale-[0.98]";
  const sizes = { md: "h-10 px-4 text-sm", lg: "h-12 px-6 text-[15px]" };
  const variants = {
    primary:
      "bg-primary text-ink-900 font-semibold hover:bg-primary-soft shadow-[0_6px_30px_-8px_rgba(34,211,238,0.6)]",
    outline:
      "border border-line text-fg hover:border-primary/60 hover:text-primary bg-white/[0.02]",
    ghost: "text-fg-muted hover:text-fg hover:bg-white/[0.04]",
  };
  const cls = cn(base, sizes[size], variants[variant], className);
  if (href)
    return (
      <Link href={href} className={cls}>
        {children}
      </Link>
    );
  return (
    <button className={cls} {...props}>
      {children}
    </button>
  );
}

export function Card({
  children,
  className,
  glow = false,
}: {
  children: React.ReactNode;
  className?: string;
  glow?: boolean;
}) {
  return (
    <div
      className={cn(
        "rounded-2xl border border-line bg-surface/60 shadow-card backdrop-blur-sm",
        glow && "shadow-glow",
        className,
      )}
    >
      {children}
    </div>
  );
}

export function Badge({
  children,
  className,
  tone = "primary",
}: {
  children: React.ReactNode;
  className?: string;
  tone?: "primary" | "bio" | "muted" | "warn";
}) {
  const tones = {
    primary: "text-primary bg-primary/10 border-primary/20",
    bio: "text-bio-soft bg-bio/10 border-bio/20",
    muted: "text-fg-muted bg-white/[0.03] border-line",
    warn: "text-warn bg-warn/10 border-warn/20",
  };
  return (
    <span
      className={cn(
        "inline-flex items-center gap-1.5 rounded-full border px-3 py-1 text-xs font-medium",
        tones[tone],
        className,
      )}
    >
      {children}
    </span>
  );
}

export function SectionHeading({
  eyebrow,
  title,
  sub,
  center = true,
}: {
  eyebrow?: string;
  title: string;
  sub?: string;
  center?: boolean;
}) {
  return (
    <div className={cn("max-w-2xl", center && "mx-auto text-center")}>
      {eyebrow && (
        <div className="mb-3 font-mono text-xs uppercase tracking-[0.2em] text-primary">
          {eyebrow}
        </div>
      )}
      <h2 className="font-display text-3xl font-semibold tracking-tight text-fg sm:text-4xl">
        {title}
      </h2>
      {sub && <p className="mt-4 text-[15px] leading-relaxed text-fg-muted">{sub}</p>}
    </div>
  );
}
