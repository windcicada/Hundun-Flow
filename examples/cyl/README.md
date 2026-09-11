# Re3900 圆柱绕流

全新初始场，直径 D=0.02 m，入口速度 2.89668 m/s，温度 300 K，出口静压 107138.196 Pa。O2/N2 空气采用 COAST EXEC/Fuels/air/mechanism 的 NASA7 与 Perry 数据，O2 质量分数为 0.23291751145757963。

网格为 448×280×40（5,017,600 单元），圆柱附近约 53–56 单元/D，展向长度 πD/2，展向间距约 0.7854 mm。x 方向采用速度入口与静压出口，y 方向对称，z 方向周期；IBM 使用自适应阶数重构，湍流模型为 Vreman 壁面函数。

网格保留原始网格在 x≈−0.01238～0.09005 m、y≈−0.03016～0.03013 m 内的全部平面节点，远场采用几何渐变。相邻远场单元的增长比上限约 1.038。`grid.py` 以 `174db3d` 的轴文件为输入生成当前网格；展向取原始节点的偶数序号。

时间步由当前流场的局部 CFL 自适应确定，目标为 0.30，`convective_cfl_margin=0.05` 给出 0.25～0.35 的调节区间。全场最大 CFL 位于区间内时保持步长，越界时按 `Δt新=Δt旧×0.30/CFL` 调整；步长增长系数为 1.10，重试缩减系数为 0.85。CFL 使用有限体积单元出流质量通量除以 ρV 的定义，MPI 汇总各子域限制；验收上限为 0.35。动量采用 CN，质量、焓与组分采用 BE。

在仓库根目录构建并运行：

```sh
cmake --build build --target hundun -j8
mpiexec --bind-to core -n 128 build/versions/v0.4/hundun run examples/cyl \
  --output /tmp/cyl --steps 1000 --output-interval 100 --restart-interval 100
```

从该新算例检查点继续：

```sh
mpiexec --bind-to core -n 128 build/versions/v0.4/hundun run examples/cyl \
  --restart /tmp/cyl/Restart --output /tmp/cyl2 --steps 1000 \
  --output-interval 100 --restart-interval 100
```

`meta.json` 记录物性、网格和输入校验值。发展段形成涡街后，使用时间加权统计与网格对比评价 Re3900 结果。

动态步长下的 SIMPLE、PISO 与默认 CN/BE 耗时见 [perf.md](perf.md)。
