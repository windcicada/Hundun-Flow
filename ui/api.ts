export async function api<T>(
  path: string,
  body?: unknown,
  signal?: AbortSignal,
): Promise<T> {
  const response = await fetch("/api" + path + (path.includes("?") ? "&" : "?") + "host_id=" + encodeURIComponent(localStorage.getItem("hundun.host") || "local"), {
    signal,
    ...(body === undefined
      ? {}
      : {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(body),
        }),
  });
  const result = await response.json();
  if (!response.ok)
    throw new Error(result.error?.message ?? `请求失败 (${response.status})`);
  return result as T;
}
export function number(value: unknown, digits = 3): string {
  if (typeof value !== "number" || !Number.isFinite(value)) return "—";
  if (value !== 0 && (Math.abs(value) < 0.001 || Math.abs(value) >= 1e7))
    return value.toExponential(2);
  return value.toLocaleString("zh-CN", { maximumFractionDigits: digits });
}
export function bytes(value: unknown): string {
  if (typeof value !== "number" || !Number.isFinite(value)) return "—";
  const units = ["B", "KiB", "MiB", "GiB", "TiB"];
  const index = Math.min(
    units.length - 1,
    Math.max(0, Math.floor(Math.log2(value || 1) / 10)),
  );
  return `${number(value / 1024 ** index, index ? 1 : 0)} ${units[index]}`;
}
export const phase: Record<string, string> = {
  unknown: "进程状态待核实",
  unverified: "进程状态待核实",
  unavailable: "记录读取异常",
  inactive: "进程已结束",
  starting: "正在启动",
  solving: "计算中",
  writing: "正在写出",
  ready: "运行中",
  retrying: "调整时间步",
  stopped: "已保存停止",
  completed: "已完成",
  failed: "运行结束 · 待检查",
  unreadable: "记录读取异常",
};
export function label(value: string): string {
  return phase[value] ?? value;
}
export function download(
  name: string,
  value: string,
  type = "text/plain;charset=utf-8",
) {
  const url = URL.createObjectURL(new Blob([value], { type }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
