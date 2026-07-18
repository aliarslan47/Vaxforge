// VaxForge backend istemcisi.
// GET uçları Next rewrites ile aynı origin'den proxy'lenir.
// AMA streaming /api/run — Next rewrite proxy'si SSE'yi BUFFER'lar (event'ler
// topluca gelir, canlı akmaz). Bu yüzden run çağrısı backend'e DOĞRUDAN gider
// (backend CORS'u açık). Backend origin: aynı host + :8011 (env ile ezilebilir).
function sseBase(): string {
  const env = process.env.NEXT_PUBLIC_API_BASE;
  if (env) return env;
  if (typeof window !== "undefined") {
    return `${window.location.protocol}//${window.location.hostname}:8011`;
  }
  return "";
}

export interface ToolStatus {
  module: string;
  label: string;
  available: boolean;
}
export interface HostInfo {
  name: string;
  label: string;
  source: string;
  n_mhc_i: number;
  n_mhc_ii: number;
  default: boolean;
}
export interface PipelineStep {
  id: string;
  module: string;
  step_tr: string;
  step_en: string;
  method: string;
  real: string;
  fallback: string | null;
  installed: boolean;
}
export interface AdjuvantInfo {
  key: string;
  label_tr: string;
  label_en: string;
  tlr: string;
  citation: string;
  desc_tr: string;
  desc_en: string;
  length: number;
}
export interface AppConfig {
  profiles: string[];
  default_profile: string;
  hosts: HostInfo[];
  default_hosts: string[];
  tools: ToolStatus[];
  steps?: PipelineStep[];
  adjuvants?: AdjuvantInfo[];
  default_adjuvant?: string;
  config_file: string;
  hosts_file: string;
}
export interface RunSummary {
  id: string;
  input: string;
  profile: string;
  timestamp: string;
  molecule: string;
  lang: string;
  hosts: string[];
  n_input: number | null;
  n_candidates: number;
  has_pdf: boolean;
  has_html: boolean;
  has_xlsx: boolean;
}
export interface Candidate {
  [key: string]: string | number;
}
export interface MevComponent {
  seq: string;
  role: "adjuvant" | "linker" | "epitope" | string;
  source?: string;
  kind?: string;
}
export interface MevData {
  sequence: string;
  cassette?: string;
  adjuvant?: string;
  adjuvant_label?: string;
  n_by_kind?: Record<string, number>;
  components?: MevComponent[];
  properties?: Record<string, any>;
  citation_keys?: string[];
}
export interface RunDetail {
  input?: string;
  profile?: string;
  timestamp?: string;
  candidates?: Candidate[];
  mev?: MevData;
  population_coverage?: any;
  iedb_match?: any;
  [key: string]: unknown;
}
export interface SSEEvent {
  phase: string;
  status: string;
  msg: string;
  data: any;
}
export interface ActiveRun {
  job_id: string;
  filename: string;
  profile: string;
  status: "running" | "done" | "error" | "cancelled";
  run_id: string | null;
  phase: string | null;
  msg: string | null;
  n_events: number;
  elapsed: number;
}

export async function getConfig(): Promise<AppConfig> {
  const r = await fetch("/api/config", { cache: "no-store" });
  if (!r.ok) throw new Error("config alınamadı");
  return r.json();
}

export async function getRuns(): Promise<RunSummary[]> {
  const r = await fetch("/api/runs", { cache: "no-store" });
  if (!r.ok) throw new Error("koşular alınamadı");
  return (await r.json()).runs;
}

export async function getRun(id: string): Promise<RunDetail> {
  const r = await fetch(`/api/runs/${id}`, { cache: "no-store" });
  if (!r.ok) throw new Error("koşu bulunamadı");
  return r.json();
}

export function fileUrl(id: string, name: string): string {
  return `/api/runs/${id}/file/${name}`;
}

export async function deleteRun(id: string): Promise<void> {
  const r = await fetch(`/api/runs/${id}`, { method: "DELETE" });
  if (!r.ok) throw new Error("koşu silinemedi");
}

/** Bir SSE yanıt gövdesini ayrıştırıp her event'i onEvent ile verir (ortak). */
async function consumeSSE(resp: Response, onEvent: (ev: SSEEvent) => void): Promise<void> {
  if (!resp.ok || !resp.body) throw new Error("akış başlatılamadı");
  const reader = resp.body.getReader();
  const decoder = new TextDecoder();
  let buf = "";
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    buf += decoder.decode(value, { stream: true });
    const parts = buf.split("\n\n");
    buf = parts.pop() ?? "";
    for (const part of parts) {
      const line = part.trim();
      if (!line.startsWith("data:")) continue;
      const json = line.slice(5).trim();
      if (!json) continue;
      try {
        onEvent(JSON.parse(json) as SSEEvent);
      } catch {
        /* yarım parça — yok say */
      }
    }
  }
}

/**
 * Pipeline'ı POST /api/run ile başlatır ve SSE event'lerini onEvent ile akıtır.
 * İlk event `__job__` job_id taşır (yenileme sonrası reconnect için saklanmalı).
 */
export async function runPipeline(
  form: FormData,
  onEvent: (ev: SSEEvent) => void,
  signal?: AbortSignal,
): Promise<void> {
  const resp = await fetch(`${sseBase()}/api/run`, { method: "POST", body: form, signal });
  await consumeSSE(resp, onEvent);
}

/** Devam eden (ya da yeni bitmiş) bir job'ın akışına imleçten yeniden bağlan. */
export async function streamJob(
  jobId: string,
  onEvent: (ev: SSEEvent) => void,
  cursor = 0,
  signal?: AbortSignal,
): Promise<void> {
  const resp = await fetch(`${sseBase()}/api/jobs/${jobId}/stream?cursor=${cursor}`, { signal });
  await consumeSSE(resp, onEvent);
}

/** Devam eden koşuların listesi (UI 'çalışıyor' rozeti için). */
export async function getActiveRuns(): Promise<ActiveRun[]> {
  const r = await fetch("/api/active-runs", { cache: "no-store" });
  if (!r.ok) return [];
  return (await r.json()).runs ?? [];
}

/** Tek bir job'ın güncel durumu (yenileme sonrası kontrol; 404 → null). */
export async function getJob(jobId: string): Promise<ActiveRun | null> {
  const r = await fetch(`/api/jobs/${jobId}`, { cache: "no-store" });
  if (!r.ok) return null;
  return r.json();
}

/** Devam eden koşuyu iptal et. */
export async function cancelRun(jobId: string): Promise<void> {
  const r = await fetch(`/api/jobs/${jobId}/cancel`, { method: "POST" });
  if (!r.ok) throw new Error("koşu iptal edilemedi");
}
