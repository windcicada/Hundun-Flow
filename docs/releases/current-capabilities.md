# 当前版本能力

公开产品版本为 `1.0.0`，内部实现/API 源码线仍命名为 `v0.4`。当前分支在 V1.0 Re=3900 预注册门完成前是 release candidate；版本字符串不等于发布接受。

候选产品路径覆盖 tensor-stretched Cartesian 网格、MPI 分解、单相低马赫理想气体、固定/自适应时间推进、入口/出口/周期/对称等外边界、静止封闭 STL IBM、WALE/Vreman、exact-history Restart、Visit 和 Evidence V6。

压力--焓主闭环在同一目标时间层联合更新 `p/h/rho/T/U` 与最终质量通量；每个接受步分别检查 EOS、continuity、energy、closed mass 和 gauge。动量修正采用 common-face owner AFC，并发布 provisional/committed CFL、retained correction 和唯一面聚合证书。IBM 默认可选 `strict_quadratic` 或 `adaptive_order`，线性降阶必须局部、满秩、可审计，最近点复制不是生产 fallback。

V1.0 的物理发布范围限于预注册的 Re=3900 `20D x 10D x 3D` 周期薄域中短程门，并明确保留 10% blockage、薄展向域和网格分辨率限制。未完成或未声明的能力包括反应流、喷雾、颗粒、移动/相交/非流形 IBM、AMR、GPU、可压缩激波、长程工程统计和任意复杂几何的统一精度保证。
