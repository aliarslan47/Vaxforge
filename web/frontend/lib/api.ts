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
export interface AppConfig {
  profiles: string[];
  default_profile: string;
  hosts: HostInfo[];
  default_hosts: string[];
  tools: ToolStatus[];
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
export interface RunDetail {
  input?: string;
  profile?: string;
  timestamp?: string;
  candidates?: Candidate[];
  [key: string]: unknown;
}
export interface SSEEvent {
  phase: string;
  status: string;
  msg: string;
  data: any;
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

/**
 * Pipeline'ı POST /api/run ile çalıştırır ve SSE event'lerini onEvent ile akıtır.
 * fetch + ReadableStream ile (EventSource POST desteklemediği için).
 */
export async function runPipeline(
  form: FormData,
  onEvent: (ev: SSEEvent) => void,
  signal?: AbortSignal,
): Promise<void> {
  const resp = await fetch(`${sseBase()}/api/run`, { method: "POST", body: form, signal });
  if (!resp.ok || !resp.body) throw new Error("pipeline başlatılamadı");
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
