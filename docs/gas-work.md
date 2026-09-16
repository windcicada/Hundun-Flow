# 气相专项执行台账

入口为 [gas.md](gas.md)。当前基线为 main 8235e59 加本专项本地改动。

| 项目 | 本次实现与证据 | 接续工作 |
|---|---|---|
| G0 | gas-ref.json：131 文件哈希、关键调用表；gas-cfl.json：完整 courant 例程的 259 输入对照 | 实际参考构建、算例后端与端到端计时窗口 |
| G1 | 三种 CFL 局部指标共同计算，默认运行参数分离为 run.json | 方向指标准入／重试／Restart、逐方程时间和固定外迭代模式 |
| G2 | Vreman 无散度反例定向回归 | ICCG、矩阵／通量、壁函数、精度及同场性能 |
| G3 | 继承 S28–S29 动态 TCR 与恢复证据 | 原始 dyn711 八步统计与分支／重新混合对齐，8/16 场动态模型组合（TCR off 复用 S24） |
| G4 | 继承 JL4 原时间步组合检查 | 化学筛选／任务均衡、热源／边界及实场轨迹 |
| I1 | 独立监看、CN 十阶段计时、接受点 stop/output、请求回执、状态查询；2→4 进程恢复通过 | 化学细分和通信包含关系，重试阶段状态通知 |
| I2 | 当前目录启动、run.json、显式覆盖、步数／时间结束，严格模式检查通过 | 简洁物理输入模板及最终方法说明扩展 |
| I3 | legacy FP64 VTK、中心坐标、分区正向重叠、原子时间索引；3/8 进程 VTK 读取与流线通过 | VisIt 应用打开、双状态附加字段目录和迁移序列关联 |
| I4 | Restart 查看、流式校验和、损坏检测通过 | 原版 input.d／历史检查点导入及来源说明 |

本台账逐项保留专项出口。各条证据按实际覆盖范围登记，完整专项验收
由 G0–G4 和 I1–I4 的后续条目共同闭合。

CFL 复现命令：

```sh
mkdir -p check/gc/f4 check/gc/f8
cp /home/wyf/code_dev/src.TCR.dyn711/courant.F90 check/gc/courant.F90
gfortran -O2 -Jcheck/gc/f4 tools/gas_cfl.f90 check/gc/courant.F90 -o check/gc/f4/probe
gfortran -O2 -fdefault-real-8 -Jcheck/gc/f8 tools/gas_cfl.f90 check/gc/courant.F90 -o check/gc/f8/probe
python3 tools/gas_cfl.py --float check/gc/f4/probe --double check/gc/f8/probe \
  --hundun b3/versions/v0.4/tests/v04_gas_cfl_probe --runner 'bash check/jam.sh' \
  --source check/gc/courant.F90 --output docs/gas-cfl.json
```

外围串行驱动把全局最大操作限定为单进程恒等操作，保持原例程的全部
局部循环和公式。默认 REAL 与 FP64 的比较分别记录，程序文件哈希
随结果保存。

接口定向检查：`check/gas-io.log`（5/5）、`check/gas-timing.log`
（CN/BE 监看）、`check/gas-fault.log`（输出路径与 MPI 故障收敛 2/2）。
`check/gas-vtk.log` 与独立 VTK 读取覆盖新版索引顺序保护；日志与
测试目录为本地生成资产，复现入口随源码保存。
