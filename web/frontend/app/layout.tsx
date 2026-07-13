import type { Metadata } from "next";
import { Exo_2, Inter, Roboto_Mono } from "next/font/google";
import "./globals.css";
import { LangProvider } from "@/components/lang-provider";

const exo = Exo_2({
  subsets: ["latin"],
  weight: ["400", "500", "600", "700"],
  variable: "--font-exo",
  display: "swap",
});
const inter = Inter({
  subsets: ["latin"],
  variable: "--font-inter",
  display: "swap",
});
const mono = Roboto_Mono({
  subsets: ["latin"],
  weight: ["400", "500", "700"],
  variable: "--font-mono",
  display: "swap",
});

export const metadata: Metadata = {
  title: "VaxForge — in silico Reverse Vaccinology",
  description:
    "Agent-assisted in silico reverse vaccinology pipeline. From pathogen to ranked vaccine candidates with a fully-cited report.",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="tr" className={`${exo.variable} ${inter.variable} ${mono.variable}`}>
      <body>
        <LangProvider>{children}</LangProvider>
      </body>
    </html>
  );
}
