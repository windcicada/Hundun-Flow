import { useEffect, useRef, useState } from "react";
import {
  Download,
  Image as ImageIcon,
  Layers3,
  LoaderCircle,
  RefreshCw,
  SlidersHorizontal,
} from "lucide-react";
import { api, number } from "./api";
import "./plots.css";

type Frame = {
  id: string;
  step: number | null;
  time: number | null;
  status: string;
  partitions?: number;
};
type Variable = { name: string; label: string; units: string; field: string };
type Fields = {
  status: string;
  frames: Frame[];
  variables: Variable[];
  images: { path: string; source: string }[];
};
type Options = {
  frame: string;
  variable: string;
  normal: string;
  coordinate: "center" | number;
};
type Rendered = {
  url: string;
  step: number | null;
  frame: string;
  time: number | null;
  variable: string;
  field: string;
  units: string;
  normal: string;
  coordinate: number;
  minimum: number;
  maximum: number;
  cache_hit: boolean;
  provenance?: {
    mask_fields?: string[];
    sources?: string[];
    time_source?: string;
  };
};
type Saved = { runId: string; key: string; options: Options; image: Rendered };
const recent = new Map<string, Saved>();
const images = new Map<string, Rendered>();
const names: Record<string, string> = {
  speed: "速度幅值",
  pressure: "压力",
  temperature: "温度",
};
const unitLabel = (units: string) =>
  units === "unspecified" ? "原始单位" : units;
const signature = (id: string, options: Options) =>
  JSON.stringify([id, options]);
const isAbort = (error: unknown) =>
  error instanceof DOMException && error.name === "AbortError";

export default function Field({
  runId,
  expanded = false,
  onRendered,
}: {
  runId: string;
  expanded?: boolean;
  onRendered?: (value: Record<string, unknown>) => void;
}) {
  const [data, setData] = useState<Fields | null>(null);
  const [loading, setLoading] = useState(true);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  const [frame, setFrame] = useState("");
  const [variable, setVariable] = useState("speed");
  const [normal, setNormal] = useState("z");
  const [center, setCenter] = useState(true);
  const [coordinate, setCoordinate] = useState("0");
  const [rendered, setRendered] = useState<Saved | null>(null);
  const [reload, setReload] = useState(0);
  const fieldRequest = useRef<AbortController | null>(null);
  const plotRequest = useRef<AbortController | null>(null);
  const callback = useRef(onRendered);
  callback.current = onRendered;

  useEffect(() => {
    const controller = new AbortController();
    fieldRequest.current?.abort();
    plotRequest.current?.abort();
    fieldRequest.current = controller;
    plotRequest.current = null;
    setData(null);
    setRendered(null);
    setError("");
    setLoading(true);
    setBusy(false);
    api<Fields>(
      `/runs/${encodeURIComponent(runId)}/fields`,
      undefined,
      controller.signal,
    )
      .then((value) => {
        if (controller.signal.aborted) return;
        setData(value);
        const ready = value.frames.filter((item) => item.status === "ready");
        const saved = recent.get(runId);
        if (
          saved &&
          ready.some((item) => item.id === saved.options.frame) &&
          value.variables.some((item) => item.name === saved.options.variable)
        ) {
          setFrame(saved.options.frame);
          setVariable(saved.options.variable);
          setNormal(saved.options.normal);
          setCenter(saved.options.coordinate === "center");
          setCoordinate(
            saved.options.coordinate === "center"
              ? "0"
              : String(saved.options.coordinate),
          );
          setRendered(saved);
        } else {
          setFrame(ready.at(-1)?.id ?? "");
          setVariable(value.variables[0]?.name ?? "speed");
          setNormal("z");
          setCenter(true);
          setCoordinate("0");
        }
      })
      .catch((reason) => {
        if (!isAbort(reason) && !controller.signal.aborted)
          setError(reason instanceof Error ? reason.message : "场数据读取失败");
      })
      .finally(() => {
        if (!controller.signal.aborted) setLoading(false);
      });
    return () => {
      controller.abort();
      plotRequest.current?.abort();
    };
  }, [runId, reload]);

  const options: Options = {
    frame,
    variable,
    normal,
    coordinate: center ? "center" : Number(coordinate),
  };
  const key = signature(runId, options);
  const image =
    rendered?.runId === runId && rendered.key === key
      ? rendered.image
      : images.get(key);
  const ready = data?.frames.filter((item) => item.status === "ready") ?? [];
  const selectedFrame = ready.find((item) => item.id === frame);
  const selectedVariable = data?.variables.find(
    (item) => item.name === variable,
  );
  const valid = Boolean(
    frame &&
      selectedVariable &&
      (center ||
        (coordinate.trim() !== "" && Number.isFinite(Number(coordinate)))),
  );

  async function generate() {
    if (!valid || busy) return;
    const requested = { ...options };
    const requestedKey = signature(runId, requested);
    const controller = new AbortController();
    plotRequest.current?.abort();
    plotRequest.current = controller;
    setBusy(true);
    setError("");
    try {
      const result = await api<Rendered>(
        `/runs/${encodeURIComponent(runId)}/plot`,
        requested,
        controller.signal,
      );
      if (controller.signal.aborted || plotRequest.current !== controller)
        return;
      const saved = {
        runId,
        key: requestedKey,
        options: requested,
        image: result,
      };
      images.set(requestedKey, result);
      recent.set(runId, saved);
      if (images.size > 24) images.delete(images.keys().next().value!);
      if (recent.size > 12) recent.delete(recent.keys().next().value!);
      setRendered(saved);
      callback.current?.(result as unknown as Record<string, unknown>);
    } catch (reason) {
      if (!isAbort(reason) && !controller.signal.aborted)
        setError(reason instanceof Error ? reason.message : "截面生成失败");
    } finally {
      if (plotRequest.current === controller && !controller.signal.aborted)
        setBusy(false);
    }
  }

  return (
    <section
      className={`panel field-panel ${expanded ? "field-expanded" : ""}`}
      aria-label="场数据后处理"
    >
      <header className="field-heading">
        <div>
          <span className="field-heading-icon">
            <Layers3 size={18} />
          </span>
          <h2>场数据截面</h2>
          <span className="field-live-label">真实输出</span>
        </div>
        <button
          className="field-icon-button"
          title="刷新可用输出"
          aria-label="刷新可用输出"
          onClick={() => setReload((value) => value + 1)}
          disabled={loading || busy}
        >
          <RefreshCw size={16} />
        </button>
      </header>
      {loading ? (
        <div className="field-empty">
          <LoaderCircle className="plot-spin" size={25} />
          <strong>正在读取场数据目录</strong>
          <p>核对可用时间步和物理变量</p>
        </div>
      ) : !data || !ready.length || !data.variables.length ? (
        <div className="field-empty">
          <ImageIcon size={30} />
          <strong>{error ? "场数据暂时不可用" : "等待场数据输出"}</strong>
          <p>
            {error ||
              "运行目录写出 VTK / Visit 场文件后，可在这里选择时间步和截面。"}
          </p>
          <button
            className="field-text-button"
            onClick={() => setReload((value) => value + 1)}
          >
            重新读取
          </button>
        </div>
      ) : (
        <>
          <div className="field-controls">
            <label>
              物理变量
              <select
                value={variable}
                disabled={busy}
                onChange={(event) => setVariable(event.target.value)}
              >
                {data.variables.map((item) => (
                  <option key={item.name} value={item.name}>
                    {names[item.name] ?? item.label} · {unitLabel(item.units)}
                  </option>
                ))}
              </select>
            </label>
            <label>
              输出时间步
              <select
                value={frame}
                disabled={busy}
                onChange={(event) => setFrame(event.target.value)}
              >
                {ready.map((item) => (
                  <option key={item.id} value={item.id}>
                    Step {item.step ?? item.id}
                    {item.time == null ? "" : ` · ${number(item.time, 9)} s`}
                  </option>
                ))}
              </select>
            </label>
            <label>
              截面法向
              <select
                value={normal}
                disabled={busy}
                onChange={(event) => setNormal(event.target.value)}
              >
                <option value="x">X 轴 · YZ 平面</option>
                <option value="y">Y 轴 · XZ 平面</option>
                <option value="z">Z 轴 · XY 平面</option>
              </select>
            </label>
            <label>
              截面位置
              <select
                value={center ? "center" : "custom"}
                disabled={busy}
                onChange={(event) => setCenter(event.target.value === "center")}
              >
                <option value="center">区域中心</option>
                <option value="custom">指定坐标</option>
              </select>
            </label>
            {!center && (
              <label className="field-position">
                {normal.toUpperCase()} 坐标 · m
                <input
                  type="number"
                  step="any"
                  value={coordinate}
                  disabled={busy}
                  onChange={(event) => setCoordinate(event.target.value)}
                  aria-label="截面坐标，单位米"
                />
              </label>
            )}
            <button
              className="field-generate"
              disabled={!valid || busy}
              onClick={generate}
            >
              {busy ? (
                <LoaderCircle className="plot-spin" size={16} />
              ) : (
                <SlidersHorizontal size={16} />
              )}
              {busy ? "生成中" : "生成截面"}
            </button>
          </div>
          <div className="field-selection">
            <span>
              {selectedFrame?.partitions
                ? `${selectedFrame.partitions} 个分区`
                : "场文件已登记"}
            </span>
            <span>字段：{selectedVariable?.field}</span>
            {data.frames.some((item) => item.status !== "ready") && (
              <span>其余输出正在等待完整发布</span>
            )}
          </div>
          {error && (
            <div role="alert" className="plot-alert">
              {error}
            </div>
          )}
          <div
            className={`field-canvas ${image ? "has-image" : ""}`}
            aria-busy={busy}
          >
            {image ? (
              <img
                src={image.url}
                alt={`${image.field}，时间步 ${image.step ?? image.frame}，${image.normal}=${number(image.coordinate, 8)}米截面`}
              />
            ) : (
              <div className="field-empty">
                <div className="field-crosshair">
                  <Layers3 size={28} />
                </div>
                <strong>按真实物理坐标查看流场</strong>
                <p>选择输出步和截面，点击“生成截面”。</p>
                <span>常用物理量与原始标量 · 真实字段读取</span>
              </div>
            )}
            {busy && (
              <div className="field-progress" role="status">
                <LoaderCircle className="plot-spin" size={25} />
                <strong>读取分区并计算截面</strong>
                <span>大规模网格需要一些时间，结果生成后自动显示</span>
                <div className="plot-progress-line" />
              </div>
            )}
          </div>
          {image && (
            <footer className="field-caption">
              <div>
                <strong>
                  {image.field} <span>[{unitLabel(image.units)}]</span>
                </strong>
                <p>
                  Step {image.step ?? image.frame} ·{" "}
                  {image.time == null
                    ? "输出未提供时间记录"
                    : `t = ${number(image.time, 9)} s`}{" "}
                  · {image.normal} = {number(image.coordinate, 8)} m
                </p>
                <p>
                  范围 {number(image.minimum, 6)} — {number(image.maximum, 6)}{" "}
                  {unitLabel(image.units)}
                  {image.provenance?.mask_fields?.length
                    ? " · 使用输出中的固体与边界掩膜"
                    : " · 输出未附固体掩膜"}
                </p>
              </div>
              <a className="field-download" href={image.url} download>
                <Download size={15} />
                下载原图
              </a>
            </footer>
          )}
        </>
      )}
    </section>
  );
}
