import { useCallback, useEffect, useRef, useState } from "react";
import {
  Activity,
  ArrowDownToLine,
  ArrowRight,
  Box,
  Check,
  ChevronDown,
  Clock3,
  Cpu,
  FileText,
  FolderOpen,
  Gauge,
  Home,
  Layers3,
  LoaderCircle,
  Menu,
  Pause,
  Play,
  Plus,
  RefreshCw,
  Search,
  Server,
  ShieldCheck,
  Timer,
  Upload,
  X,
} from "lucide-react";
import type { CaseItem, Detail, Health, Page, Run } from "./types";
import { api, bytes, download, label, number } from "./api";
import Chart from "./Chart";
import Field from "./Field";

const nav = [
  ["home", "工作台", Home],
  ["cases", "算例", FolderOpen],
  ["mesh", "网格与边界", Box],
  ["monitor", "计算监看", Activity],
  ["post", "后处理", Layers3],
  ["reports", "报告中心", FileText],
] as const;
const residuals = [
  { key: "continuity" as const, label: "连续性", color: "#1769ec" },
  { key: "energy" as const, label: "能量", color: "#f19a41" },
  { key: "eos" as const, label: "状态方程", color: "#18a899" },
];
const timing = [
  { key: "seconds" as const, label: "计算步时 / s", color: "#1599b5" },
];
function title(run: Run) {
  const parts = run.path.split("/");
  return `${parts.at(-2)} / ${run.name}`;
}
function Badge({ run }: { run: Run }) {
  return (
    <span
      className={`badge ${run.alive ? "live" : run.status === "failed" ? "warn" : "neutral"}`}
    >
      <i />
      {run.alive
        ? run.process_status === "external"
          ? "运行中 · 只读"
          : "运行中"
        : label(run.status)}
    </span>
  );
}
function Empty({ title: text, detail }: { title: string; detail?: string }) {
  return (
    <div className="empty">
      <Box size={38} />
      <h3>{text}</h3>
      <p>{detail}</p>
    </div>
  );
}
function Metric({
  name,
  value,
  unit,
  detail,
  icon: Icon,
}: {
  name: string;
  value: string;
  unit?: string;
  detail?: string;
  icon: typeof Clock3;
}) {
  return (
    <div className="metric">
      <span className="metric-icon">
        <Icon size={22} />
      </span>
      <div className="metric-label">{name}</div>
      <div className="metric-value">
        {value}
        <small>{unit}</small>
      </div>
      <div className="metric-detail">{detail || "来自当前运行记录"}</div>
    </div>
  );
}
const budgetNames: Record<string, string> = {
  mass_kg: "总质量 / kg",
  internal_energy_J: "内能 / J",
  kinetic_energy_J: "动能 / J",
  mass_balance_defect_kg_s: "质量收支缺陷 / kg·s⁻¹",
  total_energy_balance_defect_W: "能量收支缺陷 / W",
  cumulative_mass_defect_kg: "累计质量缺陷 / kg",
  cumulative_energy_defect_J: "累计能量缺陷 / J",
  continuity_equation_sum_kg_s: "连续性方程和 / kg·s⁻¹",
  energy_equation_sum_W: "能量方程和 / W",
  mass_budget_identity_kg_s: "质量核算差额 / kg·s⁻¹",
  energy_budget_identity_W: "能量核算差额 / W",
};
function Budgets({ data }: { data: Detail }) {
  const entries = Object.entries(data.budgets).filter(
    ([k, v]) => k in budgetNames && typeof v === "number",
  );
  return (
    <section className="panel budgets">
      <div className="panel-head">
        <h2>物理收支</h2>
        <ShieldCheck size={19} />
      </div>
      {entries.length ? (
        <>
          <dl className="readings">
            {entries.slice(0, 8).map(([k, v]) => (
              <div key={k}>
                <dt>{budgetNames[k]}</dt>
                <dd>{number(v, 5)}</dd>
              </div>
            ))}
          </dl>
          <p className="footnote">原始收支与方程核算分别展示</p>
        </>
      ) : (
        <Empty
          title="等待收支记录"
          detail="开启诊断输出后可查看质量与能量收支"
        />
      )}
    </section>
  );
}
function Resources({ data }: { data: Detail }) {
  const r = data.resources;
  return (
    <section className="panel">
      <div className="panel-head">
        <h2>计算资源</h2>
        <Cpu size={19} />
      </div>
      <dl className="readings">
        <div>
          <dt>MPI 进程</dt>
          <dd>{number(data.run.mpi, 0)}</dd>
        </div>
        <div>
          <dt>进程峰值 RSS</dt>
          <dd>{bytes(r?.max_rank_rss_bytes)}</dd>
        </div>
        <div>
          <dt>节点峰值 RSS</dt>
          <dd>{bytes(r?.max_node_rss_bytes)}</dd>
        </div>
        <div>
          <dt>结构化存储</dt>
          <dd>{bytes(r?.structured_bytes)}</dd>
        </div>
      </dl>
      <p className="footnote">内存高水位来自运行记录</p>
    </section>
  );
}
function Logs({ data }: { data: Detail }) {
  return (
    <section className="panel log-panel">
      <div className="panel-head">
        <h2>运行日志</h2>
        <span className="subtle">最近 {data.logs.length} 行</span>
      </div>
      <pre className="console" tabIndex={0}>
        {data.logs.length
          ? data.logs.join("\n")
          : "当前运行目录尚未提供标准输出日志\n\n时间步与状态可在上方监测图查看。"}
      </pre>
    </section>
  );
}
function Checkpoint({
  data,
  onResume,
}: {
  data: Detail;
  onResume: () => void;
}) {
  const receipt = data.controls.at(-1);
  return (
    <section className="panel">
      <div className="panel-head">
        <h2>检查点与控制</h2>
        <FolderOpen size={19} />
      </div>
      <div className="checkpoint-icon">
        <FolderOpen size={27} />
      </div>
      <h3>{data.run.has_restart ? "已发现原生检查点" : "等待检查点输出"}</h3>
      <p className="subtle">
        {data.run.pending_controls?.length
          ? "请求已提交，等待已接受时间步回执"
          : receipt
            ? `最近回执：第 ${number(receipt.step, 0)} 步 · ${receipt.published}`
            : "续算前执行原生检查与文件校验"}
      </p>
      <button
        className="button small"
        disabled={!data.run.capabilities.resume}
        onClick={onResume}
      >
        <Play size={15} />
        从检查点继续
      </button>
    </section>
  );
}
function Meta({ data }: { data: Detail }) {
  return (
    <div className="meta-strip">
      <span>
        <Box size={15} />
        {data.run.mesh?.every((x) => x != null)
          ? data.run.mesh.join(" × ")
          : "网格信息待关联"}
      </span>
      <span>
        <Layers3 size={15} />
        {data.run.scheme || "时间格式待关联"}
      </span>
      <span>
        <Activity size={15} />
        {data.run.coupling || "耦合方式待关联"}
      </span>
      {data.run.models.map((model) => (
        <span key={model}>{model}</span>
      ))}
    </div>
  );
}

export default function App() {
  const [page, setPage] = useState<Page>("home"),
    [mobile, setMobile] = useState(false),
    [search, setSearch] = useState("");
  const [runs, setRuns] = useState<Run[]>([]),
    [cases, setCases] = useState<CaseItem[]>([]),
    [health, setHealth] = useState<Health | null>(null);
  const [selected, setSelected] = useState(
      () => localStorage.getItem("hundun.run.v1") || "",
    ),
    [detail, setDetail] = useState<Detail | null>(null);
  const [error, setError] = useState(""),
    [notice, setNotice] = useState(""),
    [loading, setLoading] = useState(true),
    [connected, setConnected] = useState(false);
  const [dialog, setDialog] = useState<{
      kind: "import" | "start" | "resume" | "pause";
      caseId?: string;
    } | null>(null),
    [busy, setBusy] = useState(false);
  const [artifact, setArtifact] = useState<Record<string, unknown> | null>(
    null,
  );
  const selectedRef = useRef(selected);
  selectedRef.current = selected;
  const refresh = useCallback(async (signal?: AbortSignal) => {
    const [list, h, c] = await Promise.all([
      api<{ runs: Run[] }>("/runs", undefined, signal),
      api<Health>("/health", undefined, signal),
      api<{ cases: CaseItem[] }>("/cases", undefined, signal),
    ]);
    const sorted = [...list.runs].sort(
      (a, b) =>
        Number(b.alive) - Number(a.alive) ||
        (b.updated_at || 0) - (a.updated_at || 0),
    );
    setRuns(sorted);
    setHealth(h);
    setCases(c.cases);
    setConnected(true);
    if (!sorted.some((r) => r.id === selectedRef.current))
      setSelected(
        (
          sorted.find((r) => r.alive) ||
          sorted.find((r) => r.seconds_per_step != null) ||
          sorted.find((r) => r.has_fields) ||
          sorted[0]
        )?.id || "",
      );
    setLoading(false);
  }, []);
  useEffect(() => {
    const controller = new AbortController();
    let timer: ReturnType<typeof setTimeout>;
    let alive = true;
    const poll = async () => {
      try {
        await refresh(controller.signal);
      } catch (e) {
        if (!controller.signal.aborted) {
          setConnected(false);
          setError((e as Error).message);
          setLoading(false);
        }
      } finally {
        if (alive) timer = setTimeout(poll, 5000);
      }
    };
    void poll();
    return () => {
      alive = false;
      controller.abort();
      clearTimeout(timer);
    };
  }, [refresh]);
  useEffect(() => {
    setDetail(null);
    setArtifact(null);
    if (!selected) return;
    localStorage.setItem("hundun.run.v1", selected);
    const controller = new AbortController();
    let timer: ReturnType<typeof setTimeout>;
    let alive = true;
    const poll = async () => {
      try {
        const d = await api<Detail>(
          `/runs/${selected}`,
          undefined,
          controller.signal,
        );
        if (alive) setDetail(d);
      } catch (e) {
        if (!controller.signal.aborted) setError((e as Error).message);
      } finally {
        if (alive) timer = setTimeout(poll, 5000);
      }
    };
    void poll();
    return () => {
      alive = false;
      controller.abort();
      clearTimeout(timer);
    };
  }, [selected]);
  const run = detail?.run || runs.find((r) => r.id === selected),
    visible = runs.filter((r) =>
      (r.path + " " + r.name).toLowerCase().includes(search.toLowerCase()),
    );
  const active = runs.filter((r) => r.alive).length;
  const select = (id: string, next?: Page) => {
    setSelected(id);
    setMobile(false);
    if (next) setPage(next);
  };
  const navigate = (next: Page) => {
    setPage(next);
    setMobile(false);
  };
  const control = async (action: "output" | "pause") => {
    if (!run) return;
    setBusy(true);
    try {
      await api(`/runs/${run.id}/control`, {
        action,
        request_id: crypto.randomUUID(),
      });
      setNotice(
        action === "pause"
          ? "保存并停止请求已提交，正在等待程序回执"
          : "场输出请求已提交，正在等待程序回执",
      );
      await refresh();
      setDialog(null);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  };
  const requestAction = async (values: {
    path?: string;
    caseId?: string;
    steps?: number;
    ranks?: number;
  }) => {
    setBusy(true);
    try {
      if (dialog?.kind === "import") {
        const result = await api<{ run: Run }>("/runs/import", {
          path: values.path,
        });
        await refresh();
        select(result.run.id, "monitor");
        setNotice("已登记真实运行目录");
      } else if (dialog?.kind === "pause") {
        await control("pause");
        return;
      } else {
        const isResume = dialog?.kind === "resume";
        const result = await api<{ run: { id: string } }>(
          isResume
            ? `/runs/${selected}/resume`
            : `/cases/${values.caseId}/start`,
          {
            steps: values.steps,
            ranks: values.ranks,
            request_id: crypto.randomUUID(),
          },
        );
        await refresh();
        select(result.run.id, "monitor");
        setNotice("作业已提交，正在读取原生启动与加载结果");
      }
      setDialog(null);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  };
  const exportReport = () => {
    if (!detail) return;
    const d = detail;
    download(
      `hundun-${d.run.name}-report.md`,
      [
        `# Hundun 运行报告`,
        ``,
        `生成时间：${new Date().toISOString()}`,
        `运行目录：${d.run.path}`,
        `算例目录：${d.run.case_path || "待关联"}`,
        `状态：${label(d.run.status)}`,
        `已接受步：${number(d.run.step, 0)}`,
        `物理时间：${number(d.run.time, 9)} s`,
        `平均计算步时：${number(d.run.seconds_per_step)} s`,
        ``,
        `## 监测范围`,
        `记录数：${d.history_meta.returned_rows}；抽样：${d.history_meta.sampled ? "是" : "否"}。`,
        `计时口径：最大进程推进时间，输出成本另计。`,
        ``,
        `## 最新原始收支`,
        `\`\`\`json`,
        JSON.stringify(d.budgets, null, 2),
        "```",
        ``,
        `## 数据说明`,
        `本报告整理所选运行的原生记录；物理精度与统计充分性依据对应算例的评价目标判断。`,
        artifact ? `最新截面：${JSON.stringify(artifact, null, 2)}` : "",
      ].join("\n"),
    );
    setNotice("运行报告已导出");
  };
  const current = nav.find((x) => x[0] === page)!;
  return (
    <div className="app-shell">
      {mobile && <div className="nav-scrim" onClick={() => setMobile(false)} />}
      <aside className={`sidebar ${mobile ? "open" : ""}`}>
        <a
          className="brand"
          href="#"
          onClick={(e) => {
            e.preventDefault();
            navigate("home");
          }}
        >
          <span className="brand-mark">H</span>
          <span>
            HUNDUN<small>FLOW WORKBENCH</small>
          </span>
        </a>
        <div className="nav-label">模拟工作空间</div>
        <nav>
          {nav.map(([key, text, Icon]) => (
            <button
              key={key}
              aria-current={page === key ? "page" : undefined}
              className={`nav-item ${page === key ? "selected" : ""}`}
              onClick={() => navigate(key)}
            >
              <Icon size={20} />
              {text}
              {page === key && <span className="nav-dot" />}
            </button>
          ))}
        </nav>
        <div className="sidebar-bottom">
          <div className="current-project">
            <span>当前运行</span>
            <strong title={run?.path}>
              {run ? title(run) : "选择运行记录"}
            </strong>
            <p>{run ? label(run.status) : "登记已有算例或开始新算"}</p>
            <button onClick={() => navigate("cases")}>
              <RefreshCw size={14} />
              切换算例
            </button>
          </div>
          <div className="local-state">
            <span className={connected ? "dot green" : "dot"} />
            {connected ? "本机服务已连接" : "正在连接本机服务"}
            <Server size={15} />
          </div>
        </div>
      </aside>
      <main className="main">
        <header className="topbar">
          <div className="topbar-title">
            <button
              className="icon-button mobile-menu"
              aria-label="打开导航"
              onClick={() => setMobile(true)}
            >
              <Menu />
            </button>
            <span>仿真工作台</span>
            <span className="topbar-divider">/</span>
            <span className="topbar-current">{current[1]}</span>
          </div>
          <div className="search">
            <Search size={17} />
            <input
              aria-label="搜索运行"
              placeholder="搜索算例或运行目录"
              value={search}
              onChange={(e) => {
                setSearch(e.target.value);
                if (e.target.value) setPage("cases");
              }}
            />
          </div>
          <div className="host-state">
            <span className={connected ? "dot green" : "dot"} />
            本机计算节点
            <span className="host-divider" />
            {active} 个运行中
          </div>
        </header>
        <div className="page-content">
          {(error || notice) && (
            <div
              role={error ? "alert" : "status"}
              className={`notification ${error ? "error" : "success"}`}
            >
              <span>{error || notice}</span>
              <button
                className="icon-button"
                aria-label="关闭提示"
                onClick={() => {
                  setError("");
                  setNotice("");
                }}
              >
                <X size={17} />
              </button>
            </div>
          )}
          <div className="page-heading">
            <div>
              <h1>{page === "home" ? "工作台" : current[1]}</h1>
              <p>
                {page === "home"
                  ? "从真实运行记录出发，准备、监测与分析每一次计算。"
                  : page === "cases"
                    ? "在已登记算例上开始新算，或从已保存状态继续。"
                    : run
                      ? title(run)
                      : "选择运行，查看对应的计算数据。"}
              </p>
            </div>
            <div className="heading-actions">
              {page === "cases" || page === "home" ? (
                <>
                  <button
                    className="button"
                    onClick={() => setDialog({ kind: "import" })}
                  >
                    <FolderOpen size={17} />
                    登记运行
                  </button>
                  <button
                    className="button primary"
                    disabled={!cases.length}
                    onClick={() => setDialog({ kind: "start" })}
                  >
                    <Plus size={18} />
                    开始新算
                  </button>
                </>
              ) : page === "reports" ? (
                <button
                  className="button primary"
                  disabled={!detail}
                  onClick={exportReport}
                >
                  <ArrowDownToLine size={17} />
                  导出报告
                </button>
              ) : (
                <>
                  <button
                    className="button"
                    disabled={busy || !run?.capabilities.output}
                    title={
                      run?.capabilities.output
                        ? "在下一个已接受时间步输出场"
                        : "控制能力取决于本机程序与进程身份匹配"
                    }
                    onClick={() => void control("output")}
                  >
                    <Upload size={16} />
                    输出当前场
                  </button>
                  <button
                    className="button"
                    disabled={busy || !run?.capabilities.pause}
                    onClick={() => setDialog({ kind: "pause" })}
                  >
                    <Pause size={16} />
                    保存并暂停
                  </button>
                  <button
                    className="button primary"
                    disabled={busy || !run?.capabilities.resume}
                    onClick={() => setDialog({ kind: "resume" })}
                  >
                    <Play size={16} />
                    从检查点继续
                  </button>
                </>
              )}
            </div>
          </div>
          {loading ? (
            <div className="page-loading">
              <LoaderCircle className="spin" />
              <p>正在读取运行目录</p>
            </div>
          ) : (
            <>
              {page !== "cases" && (
                <div className="run-toolbar">
                  <div className="run-select">
                    <FolderOpen size={17} />
                    <select
                      aria-label="当前运行"
                      value={selected}
                      onChange={(e) => select(e.target.value)}
                    >
                      {runs.length === 0 && (
                        <option value="">暂无运行记录</option>
                      )}
                      {runs.map((r) => (
                        <option key={r.id} value={r.id}>
                          {title(r)}
                        </option>
                      ))}
                    </select>
                    <ChevronDown size={15} />
                  </div>
                  {run && <Badge run={run} />}
                  <span className="toolbar-note">每 5 秒更新记录</span>
                  <button
                    className="icon-button"
                    title="刷新数据"
                    aria-label="刷新数据"
                    onClick={() =>
                      void refresh().catch((e) => setError(e.message))
                    }
                  >
                    <RefreshCw size={16} />
                  </button>
                </div>
              )}
              {page === "cases" ? (
                <Cases
                  runs={visible}
                  cases={cases}
                  select={select}
                  start={(id) => setDialog({ kind: "start", caseId: id })}
                />
              ) : !detail ? (
                <section className="panel">
                  <Empty
                    title={
                      runs.length ? "正在读取所选运行" : "连接你的第一个算例"
                    }
                    detail="从已配置的目录登记运行，或选择原生算例开始计算"
                  />
                </section>
              ) : (
                <>
                  {(page === "home" || page === "monitor") && (
                    <>
                      <div className="metrics">
                        <Metric
                          name="当前步"
                          value={number(detail.run.step, 0)}
                          detail={
                            detail.run.target_steps
                              ? `目标 ${number(detail.run.target_steps, 0)} 步`
                              : "已接受时间步"
                          }
                          icon={Activity}
                        />
                        <Metric
                          name="物理时间"
                          value={number(detail.run.time, 6)}
                          unit="s"
                          detail={`当前 Δt ${number(detail.run.dt, 8)} s`}
                          icon={Clock3}
                        />
                        <Metric
                          name="平均计算步时"
                          value={number(detail.run.seconds_per_step, 2)}
                          unit="s"
                          detail="近期记录 · 输出成本另计"
                          icon={Timer}
                        />
                        <Metric
                          name="局部 CFL"
                          value={number(detail.history.at(-1)?.cfl, 3)}
                          detail={`MPI 进程 ${number(detail.run.mpi, 0)}`}
                          icon={Gauge}
                        />
                      </div>
                      <Meta data={detail} />
                      <div className="monitor-grid">
                        <Field runId={selected} onRendered={setArtifact} />
                        <div className="monitor-right">
                          <section className="panel">
                            <div className="panel-head">
                              <h2>求解监测</h2>
                              <span className="subtle">方程指标</span>
                            </div>
                            <Chart
                              data={detail.history}
                              series={residuals}
                              log
                              height={218}
                              caption="归一化残差 · 已接受时间步"
                            />
                          </section>
                          <section className="panel">
                            <div className="panel-head">
                              <h2>计算步时</h2>
                              <Clock3 size={17} />
                            </div>
                            <Chart
                              data={detail.history}
                              series={timing}
                              height={190}
                              caption="最大进程推进耗时，包含重试"
                            />
                          </section>
                        </div>
                      </div>
                      <div className="bottom-grid">
                        <Budgets data={detail} />
                        <Resources data={detail} />
                        <Checkpoint
                          data={detail}
                          onResume={() => setDialog({ kind: "resume" })}
                        />
                      </div>
                      <Logs data={detail} />
                    </>
                  )}
                  {page === "post" && (
                    <>
                      <div className="post-layout">
                        <Field
                          runId={selected}
                          expanded
                          onRendered={setArtifact}
                        />
                        <section className="panel">
                          <div className="panel-head">
                            <h2>分析范围</h2>
                            <Layers3 size={18} />
                          </div>
                          <dl className="readings">
                            <div>
                              <dt>场来源</dt>
                              <dd>已发布 VTK 输出</dd>
                            </div>
                            <div>
                              <dt>坐标</dt>
                              <dd>物理坐标 · m</dd>
                            </div>
                            <div>
                              <dt>可视化</dt>
                              <dd>指定平面截面</dd>
                            </div>
                            <div>
                              <dt>固体与 Ghost</dt>
                              <dd>按原始标记遮罩</dd>
                            </div>
                          </dl>
                          <p className="section-note">
                            选择物理量与输出时刻，生成真实截面。色标范围随当前场数据确定，图像包含步号、物理时间和单位。
                          </p>
                          <button
                            className="button full"
                            onClick={() =>
                              download(
                                `${detail.run.name}-monitor.csv`,
                                [
                                  "step,time_s,dt_s,advance_s,cfl,continuity,energy,eos",
                                  ...detail.history.map((r) =>
                                    [
                                      r.step,
                                      r.time,
                                      r.dt,
                                      r.seconds,
                                      r.cfl,
                                      r.continuity,
                                      r.energy,
                                      r.eos,
                                    ]
                                      .map((v) => v ?? "")
                                      .join(","),
                                  ),
                                ].join("\n"),
                                "text/csv;charset=utf-8",
                              )
                            }
                          >
                            <ArrowDownToLine size={16} />
                            导出监测数据
                          </button>
                        </section>
                      </div>
                      <section className="panel">
                        <div className="panel-head">
                          <h2>CFL 历史</h2>
                          <span className="subtle">
                            {detail.history_meta.sampled
                              ? "显示抽样记录"
                              : "显示已读取记录"}
                          </span>
                        </div>
                        <Chart
                          data={detail.history}
                          series={[
                            { key: "cfl", label: "局部 CFL", color: "#1769ec" },
                          ]}
                          height={250}
                        />
                      </section>
                    </>
                  )}
                  {page === "mesh" && <Configuration data={detail} />}
                  {page === "reports" && (
                    <>
                      <div className="report-heading panel">
                        <div className="report-icon">
                          <FileText size={32} />
                        </div>
                        <div>
                          <h2>运行记录报告</h2>
                          <p>汇总实际输入、时间步、原始收支和产物来源。</p>
                          <span className="subtle">{detail.run.path}</span>
                        </div>
                        <button
                          className="button"
                          onClick={() =>
                            download(
                              `${detail.run.name}-record.json`,
                              JSON.stringify(
                                { ...detail, latest_plot: artifact },
                                null,
                                2,
                              ),
                              "application/json",
                            )
                          }
                        >
                          <ArrowDownToLine size={17} />
                          完整 JSON
                        </button>
                      </div>
                      <div className="bottom-grid">
                        <Budgets data={detail} />
                        <Resources data={detail} />
                        <Checkpoint
                          data={detail}
                          onResume={() => setDialog({ kind: "resume" })}
                        />
                      </div>
                      <section className="panel">
                        <div className="panel-head">
                          <h2>运行产物</h2>
                          <span className="subtle">
                            {detail.files.length} 个文件
                          </span>
                        </div>
                        <div className="table-wrap">
                          <table>
                            <thead>
                              <tr>
                                <th>文件</th>
                                <th>大小</th>
                                <th>位置</th>
                              </tr>
                            </thead>
                            <tbody>
                              {detail.files.map((file) => (
                                <tr key={file.name}>
                                  <td>
                                    <FileText size={15} />
                                    {file.name}
                                  </td>
                                  <td>{bytes(file.size)}</td>
                                  <td className="path-cell">
                                    {detail.run.path}/{file.name}
                                  </td>
                                </tr>
                              ))}
                            </tbody>
                          </table>
                        </div>
                      </section>
                    </>
                  )}
                </>
              )}
            </>
          )}
          <footer className="page-footer">
            <span>HUNDUN-FLOW · 本机仿真工作台</span>
            <span>
              <ShieldCheck size={13} />
              数据来自原生运行记录
            </span>
          </footer>
        </div>
      </main>
      {dialog && (
        <ActionDialog
          kind={dialog.kind}
          cases={cases}
          initialCase={dialog.caseId}
          run={run}
          health={health}
          busy={busy}
          error={error}
          close={() => setDialog(null)}
          submit={requestAction}
        />
      )}
    </div>
  );
}

function Cases({
  runs,
  cases,
  select,
  start,
}: {
  runs: Run[];
  cases: CaseItem[];
  select: (id: string, page: Page) => void;
  start: (id: string) => void;
}) {
  return (
    <>
      <section className="panel">
        <div className="panel-head">
          <h2>
            运行记录 <span className="count">{runs.length}</span>
          </h2>
          <span className="subtle">历史记录与正在运行的作业</span>
        </div>
        <div className="table-wrap">
          <table className="run-table">
            <thead>
              <tr>
                <th>运行 / 算例</th>
                <th>状态</th>
                <th>当前步</th>
                <th>物理时间</th>
                <th>平均步时</th>
                <th />
              </tr>
            </thead>
            <tbody>
              {runs.map((r) => (
                <tr key={r.id}>
                  <td>
                    <button
                      className="text-button run-name"
                      onClick={() => select(r.id, "monitor")}
                    >
                      <span className="file-tile">
                        <FolderOpen size={18} />
                      </span>
                      <span>
                        {title(r)}
                        <small title={r.path}>{r.path}</small>
                      </span>
                    </button>
                  </td>
                  <td>
                    <Badge run={r} />
                  </td>
                  <td>{number(r.step, 0)}</td>
                  <td>{number(r.time, 6)} s</td>
                  <td>{number(r.seconds_per_step, 2)} s</td>
                  <td>
                    <button
                      className="icon-button"
                      aria-label={`查看 ${r.name}`}
                      onClick={() => select(r.id, "monitor")}
                    >
                      <ArrowRight size={17} />
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          {runs.length === 0 && (
            <Empty
              title="当前搜索没有匹配记录"
              detail="可以调整搜索条件，或登记其他运行目录"
            />
          )}
        </div>
      </section>
      <section className="panel">
        <div className="panel-head">
          <h2>
            原生算例 <span className="count">{cases.length}</span>
          </h2>
          <span className="subtle">新算沿用已有 case.json 与物理资产</span>
        </div>
        <div className="case-grid">
          {cases.map((c) => (
            <article className="case-item" key={c.id}>
              <div className="file-tile">
                <Box size={21} />
              </div>
              <div>
                <h3>{c.path.split("/").slice(-2).join(" / ")}</h3>
                <p title={c.path}>{c.path}</p>
              </div>
              <button
                className="button small"
                disabled={!c.capabilities.start}
                onClick={() => start(c.id)}
              >
                <Play size={14} />
                新算
              </button>
            </article>
          ))}
        </div>
      </section>
    </>
  );
}
function Configuration({ data }: { data: Detail }) {
  const config = data.case;
  return (
    <>
      <div className="config-intro panel">
        <Box size={32} />
        <div>
          <h2>网格与物理配置</h2>
          <p>来自所关联的原生算例。配置身份在启动与恢复时核对。</p>
          <span className="subtle">
            {data.run.case_path || "运行尚待关联唯一算例目录"}
          </span>
        </div>
        {config && (
          <button
            className="button"
            onClick={() =>
              download(
                "case.json",
                JSON.stringify(config, null, 2),
                "application/json",
              )
            }
          >
            <ArrowDownToLine size={16} />
            导出配置
          </button>
        )}
      </div>
      {config ? (
        <div className="config-grid">
          {Object.entries(config).map(([key, value]) => (
            <section className="panel" key={key}>
              <div className="panel-head">
                <h2>
                  {(
                    {
                      mesh: "网格与几何",
                      boundary: "边界条件",
                      boundaries: "边界条件",
                      time: "时间推进",
                      solver: "求解器",
                      reaction: "化学与燃烧",
                      spray: "液相与喷雾",
                      turbulence: "湍流模型",
                      resources: "计算资源",
                      physics: "物性",
                      ibm: "浸没边界",
                    } as Record<string, string>
                  )[key] || key}
                </h2>
                <span className="subtle">{key}</span>
              </div>
              <pre className="config-code">
                {JSON.stringify(value, null, 2)}
              </pre>
            </section>
          ))}
        </div>
      ) : (
        <section className="panel">
          <Empty
            title="等待算例关联"
            detail="运行记录保留原始数据；唯一算例目录确定后展示网格和边界配置"
          />
        </section>
      )}
    </>
  );
}
function ActionDialog({
  kind,
  cases,
  initialCase,
  run,
  health,
  busy,
  error,
  close,
  submit,
}: {
  kind: "import" | "start" | "resume" | "pause";
  cases: CaseItem[];
  initialCase?: string;
  run?: Run;
  health: Health | null;
  busy: boolean;
  error: string;
  close: () => void;
  submit: (values: {
    path?: string;
    caseId?: string;
    steps?: number;
    ranks?: number;
  }) => void;
}) {
  const ref = useRef<HTMLDialogElement>(null);
  const [path, setPath] = useState(""),
    [caseId, setCaseId] = useState(initialCase || cases[0]?.id || ""),
    [steps, setSteps] = useState(10),
    [ranks, setRanks] = useState(run?.mpi || 1);
  useEffect(() => {
    ref.current?.showModal();
  }, []);
  const titles = {
    import: "登记已有运行",
    start: "开始新算",
    resume: "从检查点继续",
    pause: "保存并暂停",
  };
  return (
    <dialog
      className="action-dialog"
      ref={ref}
      onCancel={(e) => {
        if (busy) e.preventDefault();
        else close();
      }}
    >
      <form
        onSubmit={(e) => {
          e.preventDefault();
          submit({ path, caseId, steps, ranks });
        }}
      >
        <div className="dialog-head">
          <h2>{titles[kind]}</h2>
          <button
            type="button"
            className="icon-button"
            aria-label="关闭对话框"
            disabled={busy}
            onClick={close}
          >
            <X size={20} />
          </button>
        </div>
        {kind === "import" ? (
          <>
            <p className="dialog-copy">选择已配置根目录内的真实运行目录。</p>
            <label className="form-label">
              运行目录
              <input
                autoFocus
                required
                placeholder="/path/to/run"
                value={path}
                onChange={(e) => setPath(e.target.value)}
              />
            </label>
            <div className="form-hint">
              允许的根目录：{health?.case_roots.join("、")}
            </div>
          </>
        ) : kind === "pause" ? (
          <>
            <p className="dialog-copy">
              在下一个已接受时间步保存完整检查点，然后结束当前计算。
            </p>
            <div className="operation-summary">
              <FolderOpen size={19} />
              <span>
                {run?.path}
                <small>界面将等待程序回执确认保存结果</small>
              </span>
            </div>
          </>
        ) : (
          <>
            <p className="dialog-copy">
              {kind === "resume"
                ? "沿用源算例的物理输入和完整模型历史，输出写入独立的新运行目录。"
                : "使用已登记的原生算例与物理资产，输出写入独立的新运行目录。"}
            </p>
            {kind === "start" ? (
              <label className="form-label">
                算例
                <select
                  required
                  value={caseId}
                  onChange={(e) => setCaseId(e.target.value)}
                >
                  {cases.map((c) => (
                    <option key={c.id} value={c.id}>
                      {c.path}
                    </option>
                  ))}
                </select>
              </label>
            ) : (
              <div className="operation-summary">
                <FolderOpen size={19} />
                <span>
                  {run?.path}
                  <small>Restart/current</small>
                </span>
              </div>
            )}
            <div className="form-row">
              <label className="form-label">
                推进步数
                <input
                  type="number"
                  min="1"
                  max="10000000"
                  required
                  value={steps}
                  onChange={(e) => setSteps(e.target.valueAsNumber)}
                />
              </label>
              <label className="form-label">
                MPI 进程数
                <input
                  type="number"
                  min="1"
                  max="1024"
                  required
                  value={ranks}
                  onChange={(e) => setRanks(e.target.valueAsNumber)}
                />
              </label>
            </div>
            <p className="form-hint">
              提交后先执行原生配置检查{kind === "resume" ? "与检查点校验" : ""}
              ，启动过程记录在新运行日志中。
            </p>
          </>
        )}
        <div role="alert" className={error ? "notification error" : ""}>
          {error}
        </div>
        <div className="dialog-actions">
          <button
            type="button"
            className="button"
            disabled={busy}
            onClick={close}
          >
            返回
          </button>
          <button type="submit" className="button primary" disabled={busy}>
            {busy ? (
              <LoaderCircle size={16} className="spin" />
            ) : (
              <Check size={16} />
            )}{" "}
            {busy
              ? "正在处理…"
              : kind === "pause"
                ? "发送保存请求"
                : kind === "import"
                  ? "登记运行"
                  : "提交计算"}
          </button>
        </div>
      </form>
    </dialog>
  );
}
