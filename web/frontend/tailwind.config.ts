import type { Config } from "tailwindcss";

// VaxForge tasarım token'ları — premium biotech dark (ui-ux-pro-max: Modern Dark).
// Renkler CSS değişkenlerine bağlı (globals.css); tema tutarlılığı için semantik.
const config: Config = {
  darkMode: "class",
  content: [
    "./app/**/*.{ts,tsx}",
    "./components/**/*.{ts,tsx}",
    "./lib/**/*.{ts,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        ink: {
          DEFAULT: "#0A0F1A", // derin mürekkep zemin (pure black değil — OLED smear)
          900: "#0A0F1A",
          800: "#0E1524",
          700: "#131C2E",
          600: "#1A2438",
        },
        surface: {
          DEFAULT: "#111A2B",
          raised: "#16223A",
          hover: "#1B2842",
        },
        primary: {
          DEFAULT: "#22D3EE", // electron cyan — bilim/veri
          deep: "#0891B2",
          soft: "#67E8F9",
        },
        bio: {
          DEFAULT: "#10B981", // health emerald — aşı/başarı
          soft: "#34D399",
          deep: "#059669",
        },
        // epitop tipi renk triadı (rapor stiliyle uyumlu)
        epi: {
          b: "#22D3EE", // B-hücre
          i: "#A78BFA", // MHC-I (violet)
          ii: "#34D399", // MHC-II (emerald)
        },
        line: "rgba(255,255,255,0.08)",
        fg: {
          DEFAULT: "#E6EDF6",
          muted: "#94A3B8",
          faint: "#64748B",
        },
        warn: "#F59E0B",
        danger: "#F43F5E",
      },
      fontFamily: {
        display: ["var(--font-exo)", "ui-sans-serif", "system-ui"],
        sans: ["var(--font-inter)", "ui-sans-serif", "system-ui"],
        mono: ["var(--font-mono)", "ui-monospace", "monospace"],
      },
      borderRadius: {
        xl: "0.9rem",
        "2xl": "1.25rem",
      },
      boxShadow: {
        glow: "0 0 0 1px rgba(34,211,238,0.15), 0 8px 40px -8px rgba(34,211,238,0.25)",
        card: "0 1px 0 0 rgba(255,255,255,0.04) inset, 0 12px 40px -20px rgba(0,0,0,0.8)",
      },
      keyframes: {
        "pulse-dot": {
          "0%,100%": { opacity: "1", transform: "scale(1)" },
          "50%": { opacity: "0.4", transform: "scale(0.85)" },
        },
        drift: {
          "0%,100%": { transform: "translate(0,0)" },
          "50%": { transform: "translate(3%, -4%)" },
        },
      },
      animation: {
        "pulse-dot": "pulse-dot 1.6s ease-in-out infinite",
        drift: "drift 18s ease-in-out infinite",
      },
    },
  },
  plugins: [],
};

export default config;
