import { useEffect, useMemo, useRef, useState } from "react";
import { number } from "./api";
import type { Point } from "./types";
import "./plots.css";

type Series = { key: keyof Point; label: string; color: string };
const left = 64,
  right = 20,
  top = 14,
  bottom = 34;
export default function Chart({
  data,
  series,
  log = false,
  height = 230,
  caption,
}: {
  data: Point[];
  series: Series[];
  log?: boolean;
  height?: number;
  caption?: string;
}) {
  const host = useRef<HTMLDivElement>(null);
  const [W, setWidth] = useState(740);
  useEffect(() => {
    const node = host.current;
    if (!node) return;
    const observer = new ResizeObserver((entries) =>
      setWidth(Math.max(280, Math.round(entries[0].contentRect.width))),
    );
    observer.observe(node);
    return () => observer.disconnect();
  }, []);
  const [hidden, setHidden] = useState<Set<keyof Point>>(() => new Set());
  const [hover, setHover] = useState<number | null>(null);
  const active = series.filter((item) => !hidden.has(item.key));
  const plot = useMemo(() => {
    const sorted = data
      .filter((row) => Number.isFinite(row.step))
      .sort((a, b) => a.step - b.step);
    const numeric = (row: Point, key: keyof Point) => {
      const value = row[key];
      return typeof value === "number" &&
        Number.isFinite(value) &&
        (!log || value > 0)
        ? value
        : null;
    };
    const valid = sorted.flatMap((row) =>
      active.flatMap((item) => {
        const value = numeric(row, item.key);
        return value == null ? [] : [log ? Math.log10(value) : value];
      }),
    );
    if (!valid.length) return null;
    let xmin = sorted[0].step,
      xmax = sorted.at(-1)!.step;
    if (xmin === xmax) {
      xmin -= 0.5;
      xmax += 0.5;
    }
    let ymin = Math.min(...valid),
      ymax = Math.max(...valid);
    if (ymin === ymax) {
      const padding = log ? 0.5 : Math.max(Math.abs(ymin) * 0.1, 0.5);
      ymin -= padding;
      ymax += padding;
    } else {
      const padding = (ymax - ymin) * 0.09;
      ymin -= padding;
      ymax += padding;
    }
    const x = (value: number) =>
      left + ((value - xmin) / (xmax - xmin)) * (W - left - right);
    const y = (value: number) =>
      top +
      ((ymax - (log ? Math.log10(value) : value)) / (ymax - ymin)) *
        (height - top - bottom);
    const lines = active.map((item) => {
      const segments: string[] = [];
      let segment = "";
      const dots: { x: number; y: number; index: number; value: number }[] = [];
      sorted.forEach((row, index) => {
        const value = numeric(row, item.key);
        if (value == null) {
          if (segment) segments.push(segment);
          segment = "";
          return;
        }
        const px = x(row.step),
          py = y(value);
        segment += `${segment ? "L" : "M"}${px.toFixed(2)},${py.toFixed(2)} `;
        dots.push({ x: px, y: py, index, value });
      });
      if (segment) segments.push(segment);
      return { ...item, segments, dots };
    });
    const xticks =
      sorted[0].step === sorted.at(-1)!.step
        ? [sorted[0].step]
        : [
            ...new Set(
              [0, 1, 2, 3, 4].map((tick) =>
                Math.round(xmin + (tick / 4) * (xmax - xmin)),
              ),
            ),
          ];
    return { sorted, lines, x, y, xmin, xmax, ymin, ymax, numeric, xticks };
  }, [data, series, hidden, log, height, W]);
  const hoverRow = plot && hover != null ? plot.sorted[hover] : null;
  function toggle(key: keyof Point) {
    setHidden((current) => {
      const next = new Set(current);
      if (next.has(key)) next.delete(key);
      else next.add(key);
      return next;
    });
    setHover(null);
  }
  return (
    <div className="history-chart" ref={host}>
      <div className="chart-legend">
        {series.map((item) => (
          <button
            key={item.key}
            type="button"
            className={hidden.has(item.key) ? "muted" : ""}
            aria-pressed={!hidden.has(item.key)}
            onClick={() => toggle(item.key)}
          >
            <span style={{ background: item.color }} />
            {item.label}
          </button>
        ))}
        {log && <span className="chart-scale">对数刻度 · 正值</span>}
      </div>
      {!plot ? (
        <div className="chart-empty" style={{ height }}>
          {!active.length
            ? "点击图例选择显示曲线"
            : log
              ? "当前记录尚无可显示的正值数据"
              : "等待数值记录"}
          <span>曲线来自所选运行的原始记录</span>
        </div>
      ) : (
        <div className="chart-svg-wrap">
          <svg
            className="history-svg"
            viewBox={`0 0 ${W} ${height}`}
            role="img"
            aria-label={`${series.map((item) => item.label).join("、")}随时间步变化${log ? "，纵轴为对数刻度" : ""}`}
            onMouseLeave={() => setHover(null)}
            onMouseMove={(event) => {
              const rect = event.currentTarget.getBoundingClientRect();
              const px = ((event.clientX - rect.left) / rect.width) * W;
              let nearest = 0,
                distance = Infinity;
              plot.sorted.forEach((row, index) => {
                const d = Math.abs(plot.x(row.step) - px);
                if (d < distance) {
                  distance = d;
                  nearest = index;
                }
              });
              setHover(nearest);
            }}
          >
            {[0, 1, 2, 3, 4].map((tick) => {
              const fraction = tick / 4;
              const value = plot.ymax - fraction * (plot.ymax - plot.ymin);
              const yy = top + fraction * (height - top - bottom);
              return (
                <g key={tick}>
                  <line
                    x1={left}
                    x2={W - right}
                    y1={yy}
                    y2={yy}
                    className="chart-grid"
                  />
                  <text
                    x={left - 10}
                    y={yy + 4}
                    textAnchor="end"
                    className="chart-axis"
                  >
                    {log ? `10^${number(value, 1)}` : number(value, 3)}
                  </text>
                </g>
              );
            })}
            {plot.xticks.map((value) => {
              const xx = plot.x(value);
              return (
                <text
                  key={value}
                  x={xx}
                  y={height - bottom + 18}
                  textAnchor="middle"
                  className="chart-axis"
                >
                  {number(value, 0)}
                </text>
              );
            })}
            <text
              x={W - right}
              y={height - 2}
              textAnchor="end"
              className="chart-axis"
            >
              时间步
            </text>
            {plot.lines.map((line) => (
              <g key={line.key}>
                {line.segments.map((segment, index) => (
                  <path
                    key={index}
                    d={segment}
                    fill="none"
                    stroke={line.color}
                    strokeWidth={1.8}
                    strokeLinejoin="round"
                  />
                ))}
                {line.dots.length <= 32 &&
                  line.dots.map((dot) => (
                    <circle
                      key={dot.index}
                      cx={dot.x}
                      cy={dot.y}
                      r={line.dots.length === 1 ? 3.5 : 2}
                      fill={line.color}
                    />
                  ))}
              </g>
            ))}
            {hoverRow && (
              <g>
                <line
                  x1={plot.x(hoverRow.step)}
                  x2={plot.x(hoverRow.step)}
                  y1={top}
                  y2={height - bottom}
                  stroke="#8291a7"
                  strokeDasharray="4 4"
                />
                {plot.lines.map((line) => {
                  const value = plot.numeric(hoverRow, line.key);
                  return value == null ? null : (
                    <circle
                      key={line.key}
                      cx={plot.x(hoverRow.step)}
                      cy={plot.y(value)}
                      r={3.5}
                      stroke="white"
                      strokeWidth={1.5}
                      fill={line.color}
                    />
                  );
                })}
              </g>
            )}
          </svg>
          {hoverRow && (
            <div
              className="chart-tooltip"
              style={{
                left: `${Math.min(70, Math.max(30, (plot.x(hoverRow.step) / W) * 100))}%`,
              }}
            >
              <strong>Step {number(hoverRow.step, 0)}</strong>
              {plot.lines.map((line) => (
                <span key={line.key}>
                  <i style={{ background: line.color }} />
                  {line.label}
                  <b>{number(plot.numeric(hoverRow, line.key), 7)}</b>
                </span>
              ))}
            </div>
          )}
        </div>
      )}
      {caption && <p className="chart-caption">{caption}</p>}
    </div>
  );
}
