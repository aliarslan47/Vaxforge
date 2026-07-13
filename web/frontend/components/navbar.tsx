"use client";

import Link from "next/link";
import { useEffect, useState } from "react";
import { Dna, Menu, X } from "lucide-react";
import { useLang } from "./lang-provider";
import { LANGS } from "@/lib/i18n";
import { cn } from "@/lib/utils";

export function Navbar() {
  const { lang, setLang, t } = useLang();
  const [scrolled, setScrolled] = useState(false);
  const [open, setOpen] = useState(false);

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 12);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  const links = [
    { href: "/#how", label: t("nav_method") },
    { href: "/#tools", label: t("nav_tools") },
    { href: "/#cases", label: t("nav_runs") },
  ];

  return (
    <header
      className={cn(
        "fixed inset-x-0 top-0 z-50 transition-all duration-300",
        scrolled ? "border-b border-line bg-ink/80 backdrop-blur-xl" : "border-b border-transparent",
      )}
    >
      <nav className="mx-auto flex h-16 max-w-6xl items-center justify-between px-5">
        <Link href="/" className="group flex items-center gap-2.5">
          <span className="grid h-9 w-9 place-items-center rounded-xl bg-primary/10 ring-1 ring-primary/25 transition group-hover:ring-primary/50">
            <Dna className="h-5 w-5 text-primary" strokeWidth={2} />
          </span>
          <span className="font-display text-lg font-semibold tracking-tight text-fg">
            Vax<span className="text-primary">Forge</span>
          </span>
        </Link>

        <div className="hidden items-center gap-7 md:flex">
          {links.map((l) => (
            <Link
              key={l.href}
              href={l.href}
              className="text-sm text-fg-muted transition-colors hover:text-fg"
            >
              {l.label}
            </Link>
          ))}
        </div>

        <div className="flex items-center gap-3">
          <div className="hidden items-center rounded-lg border border-line p-0.5 sm:flex">
            {LANGS.map((l) => (
              <button
                key={l}
                onClick={() => setLang(l)}
                className={cn(
                  "rounded-md px-2.5 py-1 text-xs font-medium uppercase transition-colors cursor-pointer",
                  lang === l ? "bg-primary/15 text-primary" : "text-fg-faint hover:text-fg",
                )}
                aria-pressed={lang === l}
              >
                {l}
              </button>
            ))}
          </div>
          <Link
            href="/run"
            className="hidden h-9 items-center rounded-lg bg-primary px-4 text-sm font-semibold text-ink-900 transition hover:bg-primary-soft sm:inline-flex cursor-pointer"
          >
            {t("nav_run")}
          </Link>
          <button
            className="grid h-9 w-9 place-items-center rounded-lg border border-line text-fg md:hidden cursor-pointer"
            onClick={() => setOpen((v) => !v)}
            aria-label="menu"
          >
            {open ? <X className="h-5 w-5" /> : <Menu className="h-5 w-5" />}
          </button>
        </div>
      </nav>

      {open && (
        <div className="border-t border-line bg-ink/95 px-5 py-4 md:hidden">
          <div className="flex flex-col gap-3">
            {links.map((l) => (
              <Link
                key={l.href}
                href={l.href}
                onClick={() => setOpen(false)}
                className="text-sm text-fg-muted"
              >
                {l.label}
              </Link>
            ))}
            <Link
              href="/run"
              onClick={() => setOpen(false)}
              className="mt-1 inline-flex h-10 items-center justify-center rounded-lg bg-primary text-sm font-semibold text-ink-900"
            >
              {t("nav_run")}
            </Link>
          </div>
        </div>
      )}
    </header>
  );
}
