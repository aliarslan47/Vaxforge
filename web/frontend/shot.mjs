import { chromium } from "playwright";

const BASE = process.env.BASE || "http://127.0.0.1:3000";
const OUT = "/tmp/vfshots";
import { mkdirSync } from "fs";
mkdirSync(OUT, { recursive: true });

const browser = await chromium.launch();

async function shot(name, path, { width = 1440, height = 900, full = true, wait = 1400 } = {}) {
  const ctx = await browser.newContext({
    viewport: { width, height },
    deviceScaleFactor: 2,
    colorScheme: "dark",
  });
  const page = await ctx.newPage();
  await page.goto(BASE + path, { waitUntil: "networkidle", timeout: 30000 });
  await page.waitForTimeout(wait);
  await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: full });
  console.log("✓", name, `${OUT}/${name}.png`);
  await ctx.close();
}

await shot("landing-desktop", "/", { full: true });
await shot("landing-hero", "/", { full: false });
await shot("landing-mobile", "/", { width: 390, height: 844, full: true });
await shot("run-desktop", "/run", { full: true });

await browser.close();
console.log("bitti");
