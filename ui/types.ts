export type Value = number | null | undefined;
export interface Run {
  id: string;
  name: string;
  path: string;
  case_path: string | null;
  status: string;
  step: Value;
  time: Value;
  dt: Value;
  mpi: Value;
  mesh: (number | null)[] | null;
  scheme: string | null;
  coupling: string | null;
  models: string[];
  seconds_per_step: Value;
  target_steps: Value;
  alive: boolean | null;
  updated_at: Value;
  has_restart: boolean;
  has_fields?: boolean;
  case_binding?: string;
  process_status?: string;
  capabilities: {
    output: boolean;
    pause: boolean;
    resume: boolean;
    start: boolean;
  };
  pending_controls?: string[];
  error?: string;
}
export interface Point {
  step: number;
  time: Value;
  seconds: Value;
  dt?: Value;
  continuity: Value;
  energy: Value;
  eos: Value;
  cfl: Value;
  pressure_iterations: Value;
  outer: Value;
  cd?: Value;
  cl?: Value;
}
export interface Detail {
  run: Run;
  history: Point[];
  history_meta: {
    sampled: boolean;
    available_rows: number;
    returned_rows: number;
    timing_scope: string;
    tail_limited?: boolean;
    first_step?: number;
    last_step?: number;
  };
  budgets: Record<string, unknown>;
  resources: Record<string, unknown> | null;
  logs: string[];
  controls: {
    request: string;
    step: number;
    time: number;
    published: string;
  }[];
  files: { name: string; size: number }[];
  case: Record<string, unknown> | null;
}
export interface CaseItem {
  id: string;
  name: string;
  path: string;
  config: Record<string, unknown>;
  capabilities: { start: boolean };
}
export interface Health {
  ok: boolean;
  program: string;
  program_available: boolean;
  case_roots: string[];
  state: string;
  api_version: number;
}
export type Page = "home" | "cases" | "mesh" | "monitor" | "post" | "reports";
