# Re3900 圆柱绕流

全新初始场，直径 D=0.02 m，入口速度 2.89668 m/s，温度 300 K，出口静压 107138.196 Pa。O2/N2 空气采用 COAST EXEC/Fuels/air/mechanism 的 NASA7 与 Perry 数据，O2 质量分数为 0.23291751145757963。

网格为 576×320×80（14,745,600 单元），圆柱附近约 53–56 单元/D，展向长度 πD/2。x 方向采用速度入口与静压出口，y 方向对称，z 方向周期；IBM 使用自适应阶数重构，湍流模型为 Vreman 壁面函数。

时间步由当前流场 CFL=0.5 自适应确定，初值随流场约束调整；自动重试用于启动和瞬态变化。动量采用 CN，质量、焓与组分采用 BE。

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
