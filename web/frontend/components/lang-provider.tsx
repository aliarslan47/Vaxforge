"use client";

import { createContext, useContext, useEffect, useState } from "react";
import { Lang, t as translate } from "@/lib/i18n";

interface LangCtx {
  lang: Lang;
  setLang: (l: Lang) => void;
  t: (key: string) => string;
}

const Ctx = createContext<LangCtx>({
  lang: "tr",
  setLang: () => {},
  t: (k) => k,
});

export function LangProvider({ children }: { children: React.ReactNode }) {
  const [lang, setLangState] = useState<Lang>("tr");

  useEffect(() => {
    const saved = (typeof window !== "undefined" &&
      window.localStorage.getItem("vf_lang")) as Lang | null;
    if (saved === "tr" || saved === "en") setLangState(saved);
  }, []);

  const setLang = (l: Lang) => {
    setLangState(l);
    try {
      window.localStorage.setItem("vf_lang", l);
      document.documentElement.lang = l;
    } catch {
      /* ignore */
    }
  };

  return (
    <Ctx.Provider value={{ lang, setLang, t: (k) => translate(lang, k) }}>
      {children}
    </Ctx.Provider>
  );
}

export const useLang = () => useContext(Ctx);
