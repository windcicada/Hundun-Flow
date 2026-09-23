import { useEffect, useState } from "react";
import {
  Bot,
  Settings,
  X,
  Send,
  Plus,
  Server,
  FilePlus,
  Archive,
} from "lucide-react";
import App from "./App";
import { api } from "./api";
import "./client.css";
type Obj = Record<string, any>;
const status: Record<string, string> = {
  queued: "准备中",
  running: "执行中",
  waiting: "监测作业",
  needs_model: "配置模型后继续",
  needs_input: "等待补充",
  completed: "已完成",
  cancelled: "已取消",
  failed: "需要检查",
};
export default function Client() {
  const [host, setHost] = useState(
    localStorage.getItem("hundun.host") || "local",
  );
  const [hosts, setHosts] = useState<Obj[]>([]),
    [panel, setPanel] = useState(""),
    [profile, setProfile] = useState<Obj>({}),
    [model, setModel] = useState<Obj>({});
  const [tasks, setTasks] = useState<Obj[]>([]),
    [task, setTask] = useState(""),
    [events, setEvents] = useState<Obj[]>([]),
    [message, setMessage] = useState("");
  const [error, setError] = useState(""),
    [notice, setNotice] = useState(""),
    [busy, setBusy] = useState(false),
    [ranks, setRanks] = useState(1),
    [steps, setSteps] = useState(100);
  const [source, setSource] = useState(""),
    [config, setConfig] = useState(""),
    [revisions, setRevisions] = useState<Obj[]>([]),
    [revision, setRevision] = useState<Obj | null>(null);
  const [advanced, setAdvanced] = useState("{}");
  const [exportRun, setExportRun] = useState(""),
    [exportRuns, setExportRuns] = useState<Obj[]>([]),
    [includeCheckpoint, setIncludeCheckpoint] = useState(true),
    [runtime, setRuntime] = useState(""),
    [artifact, setArtifact] = useState<Obj | null>(null);
  useEffect(() => {
    if (panel === "exports")
      api<{ runs: Obj[] }>("/runs")
        .then((d) => {
          setExportRuns(d.runs);
          setExportRun(
            localStorage.getItem("hundun.run.v1." + host) ||
              d.runs[0]?.id ||
              "",
          );
        })
        .catch((e) => setError(e.message));
  }, [panel, host]);
  const current = tasks.find((t) => t.id === task);
  async function loadHosts() {
    const d = await api<{ hosts: Obj[] }>("/hosts");
    setHosts(d.hosts);
    setProfile(d.hosts.find((h) => h.id === host) || d.hosts[0]);
  }
  useEffect(() => {
    loadHosts().catch((e) => setError(e.message));
    api<Obj>("/model")
      .then(setModel)
      .catch((e) => setError(e.message));
  }, [host]);
  useEffect(() => {
    let live = true;
    const update = () =>
      api<{ tasks: Obj[] }>("/tasks")
        .then((d) => {
          if (live) setTasks(d.tasks);
        })
        .catch(() => {});
    update();
    const timer = setInterval(update, 3000);
    return () => {
      live = false;
      clearInterval(timer);
    };
  }, []);
  useEffect(() => {
    if (!task) return;
    setEvents([]);
    const stream = new EventSource(
      "/api/tasks/" + task + "/events?stream=true",
    );
    stream.onmessage = (e) => {
      const item = JSON.parse(e.data);
      setEvents((old) =>
        old.some((x) => x.seq === item.seq) ? old : [...old, item].slice(-500),
      );
    };
    return () => stream.close();
  }, [task]);
  useEffect(() => {
    if (panel === "inputs")
      api<{ revisions: Obj[] }>("/revisions")
        .then((d) => setRevisions(d.revisions))
        .catch((e) => setError(e.message));
  }, [panel, host]);
  useEffect(() => {
    if (!panel) return;
    const old = document.activeElement as HTMLElement;
    const dialog = document.querySelector<HTMLElement>(".client-panel");
    dialog?.focus();
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") setPanel("");
      if (e.key === "Tab" && dialog) {
        const nodes = Array.from(
          dialog.querySelectorAll<HTMLElement>(
            "button,input,select,textarea,a[href]",
          ),
        ).filter((x) => !x.hasAttribute("disabled"));
        const a = nodes[0],
          b = nodes[nodes.length - 1];
        if (e.shiftKey && document.activeElement === a) {
          e.preventDefault();
          b?.focus();
        } else if (!e.shiftKey && document.activeElement === b) {
          e.preventDefault();
          a?.focus();
        }
      }
    };
    document.addEventListener("keydown", onKey);
    return () => {
      document.removeEventListener("keydown", onKey);
      old?.focus();
    };
  }, [panel]);
  async function act(fn: () => Promise<void>) {
    setBusy(true);
    setError("");
    setNotice("");
    try {
      await fn();
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }
  function choose(id: string) {
    localStorage.setItem("hundun.host", id);
    localStorage.setItem(
      "hundun.host.name",
      hosts.find((h) => h.id === id)?.name || id,
    );
    setHost(id);
    setRevision(null);
  }
  async function operation(action: string, params: Obj = {}) {
    return api<Obj>("/operations", {
      host_id: host,
      action,
      params: { ...params, request_id: crypto.randomUUID() },
    });
  }
  async function submit() {
    await act(async () => {
      if (task) {
        await api("/tasks/" + task + "/followup", {
          message,
          request_id: crypto.randomUUID(),
        });
      } else {
        const t = await api<Obj>("/tasks", {
          host_id: host,
          message,
          max_ranks: ranks,
          max_steps: steps,
          request_id: crypto.randomUUID(),
        });
        setTask(t.id);
      }
      setMessage("");
    });
  }
  return (
    <>
      <App key={host} />
      <div className="client-tools">
        <label>
          <Server size={16} />
          <select
            aria-label="计算主机"
            value={host}
            onChange={(e) => choose(e.target.value)}
          >
            {hosts.map((h) => (
              <option key={h.id} value={h.id}>
                {h.name || h.id}
              </option>
            ))}
          </select>
        </label>
        <button
          onClick={() => {
            setPanel("inputs");
            setError("");
          }}
        >
          <FilePlus size={17} />
          输入版本
        </button>
        <button
          onClick={() => {
            setPanel("settings");
            setError("");
            setAdvanced(
              JSON.stringify(
                {
                  tools: profile.tools || {},
                  known_hosts: profile.known_hosts || "",
                  key_file: profile.key_file || "",
                  password_ref: profile.password_ref || "",
                  python: profile.python || "python3",
                },
                null,
                2,
              ),
            );
          }}
          aria-label="连接与模型"
        >
          <Settings size={18} />
        </button>
        <button
          aria-label="报告与归档"
          onClick={() => {
            setPanel("exports");
            setArtifact(null);
            setError("");
          }}
        >
          <Archive size={18} />
        </button>
        <button
          className="ai-open"
          onClick={() => {
            setPanel("ai");
            setError("");
          }}
        >
          <Bot size={19} />
          AI 助手
        </button>
      </div>
      {panel && (
        <div
          className="client-shade"
          onMouseDown={(e) => {
            if (e.target === e.currentTarget) setPanel("");
          }}
        >
          <section
            className="client-panel"
            role="dialog"
            aria-modal="true"
            aria-label={
              panel === "ai"
                ? "AI 模拟助手"
                : panel === "settings"
                  ? "连接与模型"
                  : panel === "exports"
                    ? "报告与归档"
                    : "输入版本"
            }
            tabIndex={-1}
          >
            <header>
              <div>
                <small>HUNDUN CLIENT · 0.1.0</small>
                <h2>
                  {panel === "ai"
                    ? "AI 模拟助手"
                    : panel === "settings"
                      ? "连接与模型"
                      : panel === "exports"
                        ? "报告与归档"
                        : "输入版本"}
                </h2>
              </div>
              <button aria-label="关闭面板" onClick={() => setPanel("")}>
                <X />
              </button>
            </header>
            {error && (
              <p role="alert" className="client-error">
                {error}
              </p>
            )}
            {notice && (
              <p role="status" className="client-notice">
                {notice}
              </p>
            )}
            {panel === "settings" && (
              <div className="client-scroll">
                <h3>计算连接</h3>
                <p>选择工作目录与外部程序。已有结果和输入可独立管理。</p>
                <label>
                  主机配置
                  <select
                    value={profile.id || ""}
                    onChange={(e) => {
                      const h = hosts.find((h) => h.id === e.target.value)!;
                      setProfile(h);
                      setAdvanced(
                        JSON.stringify(
                          {
                            tools: h.tools || {},
                            known_hosts: h.known_hosts || "",
                            key_file: h.key_file || "",
                            password_ref: h.password_ref || "",
                            python: h.python || "python3",
                          },
                          null,
                          2,
                        ),
                      );
                    }}
                  >
                    {hosts.map((h) => (
                      <option key={h.id} value={h.id}>
                        {h.name || h.id}
                      </option>
                    ))}
                    {!profile.id && <option value="">新连接</option>}
                  </select>
                </label>
                <button
                  onClick={() => {
                    setProfile({
                      name: "计算节点",
                      kind: "ssh",
                      roots: [],
                      max_ranks: 4,
                      max_steps: 100000,
                    });
                    setAdvanced("{}");
                  }}
                >
                  <Plus size={16} />
                  添加连接
                </button>
                <div className="client-grid">
                  <label>
                    名称
                    <input
                      value={profile.name || ""}
                      onChange={(e) =>
                        setProfile({ ...profile, name: e.target.value })
                      }
                    />
                  </label>
                  <label>
                    类型
                    <select
                      value={profile.kind || "local"}
                      onChange={(e) =>
                        setProfile({ ...profile, kind: e.target.value })
                      }
                    >
                      <option value="local">本机 Linux</option>
                      <option value="ssh">SSH Linux</option>
                    </select>
                  </label>
                </div>
                {profile.kind === "ssh" && (
                  <>
                    <div className="client-grid">
                      {["hostname", "username", "port"].map((k) => (
                        <label key={k}>
                          {
                            (
                              {
                                hostname: "主机地址",
                                username: "SSH 用户",
                                port: "端口",
                              } as Obj
                            )[k]
                          }
                          <input
                            value={profile[k] || ""}
                            onChange={(e) =>
                              setProfile({
                                ...profile,
                                [k]:
                                  k === "port"
                                    ? Number(e.target.value)
                                    : e.target.value,
                              })
                            }
                          />
                        </label>
                      ))}
                    </div>
                    <label>
                      远端客户端工作目录
                      <input
                        placeholder="/home/user/client"
                        value={profile.work_dir || ""}
                        onChange={(e) =>
                          setProfile({ ...profile, work_dir: e.target.value })
                        }
                      />
                    </label>
                  </>
                )}
                <label>
                  算例工作目录（每行一个）
                  <textarea
                    rows={2}
                    value={(profile.roots || []).join("\n")}
                    onChange={(e) =>
                      setProfile({
                        ...profile,
                        roots: e.target.value.split("\n").filter(Boolean),
                      })
                    }
                  />
                </label>
                <label>
                  外部求解器入口
                  <input
                    placeholder="/opt/hundun/run"
                    value={profile.program || ""}
                    onChange={(e) =>
                      setProfile({
                        ...profile,
                        program: e.target.value || null,
                      })
                    }
                  />
                </label>
                <label>
                  MPI 入口
                  <input
                    placeholder="/opt/hundun/mpirun"
                    value={profile.mpi || ""}
                    onChange={(e) =>
                      setProfile({ ...profile, mpi: e.target.value || null })
                    }
                  />
                </label>
                <div className="client-grid">
                  <label>
                    MPI 上限
                    <input
                      type="number"
                      min={1}
                      value={profile.max_ranks || 1}
                      onChange={(e) =>
                        setProfile({
                          ...profile,
                          max_ranks: Number(e.target.value),
                        })
                      }
                    />
                  </label>
                  <label>
                    单次步数上限
                    <input
                      type="number"
                      min={1}
                      value={profile.max_steps || 100}
                      onChange={(e) =>
                        setProfile({
                          ...profile,
                          max_steps: Number(e.target.value),
                        })
                      }
                    />
                  </label>
                </div>
                <details>
                  <summary>SSH 身份与外部工具</summary>
                  <p>
                    known_hosts、key_file、password_ref
                    使用本机路径或凭据引用；tools 使用命令参数数组，{"{case}"}{" "}
                    表示输入目录。
                  </p>
                  <textarea
                    rows={7}
                    aria-label="高级主机配置"
                    value={advanced}
                    onChange={(e) => setAdvanced(e.target.value)}
                  />
                </details>
                <div className="client-actions">
                  <button
                    disabled={busy}
                    onClick={() =>
                      act(async () => {
                        const h = await api<Obj>("/hosts", {
                          ...profile,
                          ...JSON.parse(advanced),
                        });
                        await loadHosts();
                        setProfile(h);
                        setNotice("主机配置已保存");
                      })
                    }
                  >
                    保存连接
                  </button>
                  <button
                    disabled={busy || !profile.id}
                    onClick={() =>
                      act(async () => {
                        const d = await api<Obj>(
                          "/hosts/" + profile.id + "/probe",
                          {},
                        );
                        setNotice(
                          `连接已就绪 · ${d.cpu_count || "—"} 个逻辑核\n读取与绘图已启用 · ${d.program_available ? "计算已启用" : "配置程序后启用计算"}\n${d.version || ""}`,
                        );
                      })
                    }
                  >
                    探测能力
                  </button>
                </div>
                <h3>模型服务</h3>
                <label>
                  兼容端点
                  <input
                    placeholder="http://localhost:8000/v1"
                    value={model.base_url || ""}
                    onChange={(e) =>
                      setModel({ ...model, base_url: e.target.value })
                    }
                  />
                </label>
                <label>
                  模型名称
                  <input
                    value={model.model || ""}
                    onChange={(e) =>
                      setModel({ ...model, model: e.target.value })
                    }
                  />
                </label>
                <label>
                  凭据引用
                  <input
                    placeholder="env:HUNDUN_MODEL_KEY 或 keyring:my-model"
                    value={model.credential_ref || ""}
                    onChange={(e) =>
                      setModel({ ...model, credential_ref: e.target.value })
                    }
                  />
                </label>
                <p>凭据由启动客户端的环境变量或系统凭据库提供。</p>
                <button
                  disabled={busy}
                  onClick={() =>
                    act(async () => {
                      await api("/model", model);
                      setNotice("模型配置已保存，可回到任务面板继续");
                    })
                  }
                >
                  保存模型
                </button>
              </div>
            )}
            {panel === "inputs" && (
              <div className="client-scroll">
                <p>
                  当前主机：{hosts.find((h) => h.id === host)?.name}
                  。每次编辑保存为独立输入版本，保留来源、资产校验与修改差异。
                </p>
                <label>
                  导入目录
                  <input
                    placeholder="已有 case.json 所在目录；留空使用原生模板"
                    value={source}
                    onChange={(e) => setSource(e.target.value)}
                  />
                </label>
                <label>
                  完整 case.json（可选）
                  <textarea
                    rows={8}
                    value={config}
                    onChange={(e) => setConfig(e.target.value)}
                    placeholder="导入后可在这里编辑配置"
                  />
                </label>
                <button
                  disabled={busy}
                  onClick={() =>
                    act(async () => {
                      const d = await operation("prepare", {
                        source,
                        ...(config ? { config: JSON.parse(config) } : {}),
                        reason: "工作台输入编辑",
                      });
                      setRevision(d.revision);
                      const list = await operation("revisions");
                      setRevisions(list.revisions);
                      setNotice("已创建独立输入版本");
                    })
                  }
                >
                  创建输入版本
                </button>
                <h3>版本记录</h3>
                {revisions.map((r) => (
                  <button
                    className="client-row"
                    key={r.id}
                    onClick={() => {
                      setRevision(r);
                      setSource(r.path);
                      setConfig(JSON.stringify(r.config, null, 2));
                    }}
                  >
                    <span>
                      {r.id}
                      <small>{r.reason}</small>
                    </span>
                    <span>{new Date(r.created * 1000).toLocaleString()}</span>
                  </button>
                ))}
                {revision && (
                  <>
                    <h3>修改差异</h3>
                    <pre>{revision.diff || "与来源配置相同"}</pre>
                    <div className="client-actions">
                      <button
                        disabled={busy}
                        onClick={() =>
                          act(async () => {
                            const d = await operation("case_check", {
                              case_id: revision.id,
                              ranks,
                            });
                            setNotice(d.output);
                          })
                        }
                      >
                        原生检查
                      </button>
                      <button
                        disabled={busy}
                        onClick={() =>
                          act(async () => {
                            const d = await operation("start", {
                              case_id: revision.id,
                              ranks,
                              steps,
                            });
                            setNotice("新算已派发：" + d.run.id);
                          })
                        }
                      >
                        新算
                      </button>
                    </div>
                    <div className="client-grid">
                      <label>
                        MPI
                        <input
                          type="number"
                          min={1}
                          value={ranks}
                          onChange={(e) => setRanks(Number(e.target.value))}
                        />
                      </label>
                      <label>
                        步数
                        <input
                          type="number"
                          min={1}
                          value={steps}
                          onChange={(e) => setSteps(Number(e.target.value))}
                        />
                      </label>
                    </div>
                  </>
                )}
              </div>
            )}
            {panel === "exports" && (
              <div className="client-scroll">
                <p>导出真实运行记录与选定检查点。完整归档在计算端生成。</p>
                <label>
                  运行
                  <select
                    value={exportRun}
                    onChange={(e) => {
                      setExportRun(e.target.value);
                      setArtifact(null);
                    }}
                  >
                    {exportRuns.map((r) => (
                      <option key={r.id} value={r.id}>
                        {r.name} · {r.id.slice(0, 8)}
                      </option>
                    ))}
                  </select>
                </label>
                <label>
                  检查点
                  <select
                    value={includeCheckpoint ? "yes" : "no"}
                    onChange={(e) =>
                      setIncludeCheckpoint(e.target.value === "yes")
                    }
                  >
                    <option value="yes">包含已接受检查点（保存停止后）</option>
                    <option value="no">输入与运行记录</option>
                  </select>
                </label>
                <label>
                  外部程序发行包目录（可选）
                  <input
                    value={runtime}
                    onChange={(e) => setRuntime(e.target.value)}
                  />
                </label>
                <div className="client-actions">
                  <button
                    disabled={busy || !exportRun}
                    onClick={() =>
                      act(async () =>
                        setArtifact(
                          await operation("report", { run_id: exportRun }),
                        ),
                      )
                    }
                  >
                    导出运行记录
                  </button>
                  <button
                    disabled={busy || !exportRun}
                    onClick={() =>
                      act(async () =>
                        setArtifact(
                          await operation("pack", {
                            run_id: exportRun,
                            checkpoint: includeCheckpoint,
                            runtime,
                          }),
                        ),
                      )
                    }
                  >
                    生成算例归档
                  </button>
                </div>
                {artifact && (
                  <>
                    <p>产物：{artifact.path}</p>
                    {artifact.sha256 && <pre>SHA256 {artifact.sha256}</pre>}
                    {artifact.url ? (
                      <a href={artifact.url} download>
                        下载产物
                      </a>
                    ) : (
                      <p>通过计算端 SFTP 下载大型归档。</p>
                    )}
                  </>
                )}
              </div>
            )}
            {panel === "ai" && (
              <>
                <div className="client-taskbar">
                  <select
                    aria-label="AI 任务"
                    value={task}
                    onChange={(e) => setTask(e.target.value)}
                  >
                    <option value="">新任务</option>
                    {tasks.map((t) => (
                      <option key={t.id} value={t.id}>
                        {status[t.status] || t.status} ·{" "}
                        {t.message.slice(0, 32)}
                      </option>
                    ))}
                  </select>
                  <button
                    aria-label="新建 AI 任务"
                    onClick={() => {
                      setTask("");
                      setEvents([]);
                    }}
                  >
                    <Plus />
                  </button>
                </div>
                {current && (
                  <div className="client-taskinfo">
                    <strong>{status[current.status] || current.status}</strong>
                    <span>
                      主机 {current.host_id} · {current.jobs.length} 个关联作业
                    </span>
                    <button
                      disabled={busy || current.status === "cancelled"}
                      onClick={() =>
                        act(async () => {
                          await api("/tasks/" + task + "/cancel", {
                            request_id: crypto.randomUUID(),
                          });
                          setNotice("AI 已取消；计算作业继续按自身状态运行");
                        })
                      }
                    >
                      取消 AI 任务
                    </button>
                    {current.jobs.map((id: string) => (
                      <button
                        key={id}
                        onClick={() => {
                          choose(current.host_id);
                          localStorage.setItem(
                            "hundun.run.v1." + current.host_id,
                            id,
                          );
                          window.location.reload();
                        }}
                      >
                        查看运行 {id.slice(0, 8)}
                      </button>
                    ))}
                  </div>
                )}
                <div className="client-events" aria-live="polite">
                  {!task && (
                    <div className="client-welcome">
                      <Bot size={36} />
                      <h3>从模拟目标开始</h3>
                      <p>
                        助手通过与工作台相同的任务接口准备输入、运行计算、绘图和归档。计算端定时监测进度，在状态变化后交给助手分析。
                      </p>
                    </div>
                  )}
                  {events.map((e) => (
                    <article key={e.seq} className={"event-" + e.type}>
                      <small>
                        {(
                          {
                            created: "任务目标",
                            needs_model: "模型连接",
                            cancelled: "任务取消",
                            operation: "执行结果",
                            tool: "调用工具",
                            tool_done: "工具完成",
                            answer: "助手回复",
                            text: "助手",
                            completed: "任务完成",
                            question: "待补充",
                            failed: "执行提示",
                            user: "补充目标",
                          } as Obj
                        )[e.type] || e.type}{" "}
                        · {new Date(e.time * 1000).toLocaleTimeString("zh-CN")}
                      </small>
                      <div>{e.text || e.message || e.name || e.action}</div>
                      {e.arguments && (
                        <pre>{JSON.stringify(e.arguments, null, 2)}</pre>
                      )}
                      {e.result && (
                        <details>
                          <summary>操作结果</summary>
                          <pre>{JSON.stringify(e.result, null, 2)}</pre>
                          {e.result.url && (
                            <a
                              href={e.result.url}
                              target="_blank"
                              rel="noreferrer"
                            >
                              打开产物
                            </a>
                          )}
                        </details>
                      )}
                    </article>
                  ))}
                </div>
                <form
                  className="client-compose"
                  onSubmit={(e) => {
                    e.preventDefault();
                    submit();
                  }}
                >
                  {!task && (
                    <div className="client-grid">
                      <label>
                        MPI 上限
                        <input
                          type="number"
                          min={1}
                          value={ranks}
                          onChange={(e) => setRanks(Number(e.target.value))}
                        />
                      </label>
                      <label>
                        单次步数上限
                        <input
                          type="number"
                          min={1}
                          value={steps}
                          onChange={(e) => setSteps(Number(e.target.value))}
                        />
                      </label>
                    </div>
                  )}
                  <textarea
                    aria-label="模拟目标或追问"
                    placeholder="描述目标、算例目录和物理条件…"
                    rows={3}
                    value={message}
                    onChange={(e) => setMessage(e.target.value)}
                  />
                  <button disabled={busy || !message.trim()} type="submit">
                    <Send size={16} />
                    {task ? "补充并继续" : "开始任务"}
                  </button>
                </form>
              </>
            )}
          </section>
        </div>
      )}
    </>
  );
}
