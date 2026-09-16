# 跨主机环境与运行

本页区分源码构建、反应计算依赖和独立运行包。主仓库包含求解器、
小型示例及依赖身份；**截至 2026-09-16，GitHub Release 尚未提供锁定的
Cantera SDK 或独立运行包下载资产**。仅 clone 主仓库不足以运行真实化学。
请先取得维护者提供的配套 SDK/运行包，再执行下面的命令。

## 已验证环境

| 项目 | 反应计算参考组合 |
| --- | --- |
| 系统与架构 | Linux x86-64，Ubuntu 22.04 / glibc 2.35 |
| 编译器 | Clang 15，C++17，GCC 11 的 libstdc++ 开发文件，ABI1 |
| 构建工具 | CMake ≥ 3.21、Ninja、lld 15 |
| 并行 | OpenMPI 4，满足 MPI ≥ 3；编译和运行使用同一实现 |
| 化学 | 锁定的 Cantera 3.2.0 C++ SDK，见下节 |
| 分析工具 | Python 3.8+；迁移/数值分析脚本还需要 NumPy |
| 数值选项 | FP64、Release、`-ffp-contract=off`、ThinLTO；不使用 fast-math |

普通核心构建可使用其他符合 C++17/MPI3 的工具链；这不代表该组合已经
通过本文的反应计算及跨主机验收。不要混用 Clang/libc++、libstdc++ ABI0，
也不要用宿主 MPICH 启动链接到另一套 OpenMPI 的程序。

在标准 Ubuntu 22.04 环境中，可由管理员安装构建工具：

```sh
sudo apt-get update
sudo apt-get install git cmake ninja-build clang-15 lld-15 g++-11 \
  libopenmpi-dev openmpi-bin python3 python3-numpy
```

管理员安装命令不适用于共享主机上的无权限用户；此时使用已准备的隔离
构建环境或完整独立运行包。不要替换宿主 `/lib` 中的 glibc。

## Cantera SDK：版本号相同不等于二进制身份相同

`pip install cantera` 不能替代本项目要求的 C++ SDK。普通 configure 不下载
依赖、不编译第三方库，并对库及许可证校验 SHA-256。锁定信息：

- [SDK 清单](../../third_party/cantera/PREBUILT-LINUX-X86_64.json)
- [上游与依赖身份](../../third_party/cantera/UPSTREAM.json)
- [实际准入检查](../../cmake/HundunPortableCantera.cmake)

SDK 根目录应包含 `include/cantera/`、`lib/`、`share/cantera/data/` 和
`licenses/`。`libcantera_shared.so`、`libcantera_shared.so.3` 必须为指向
`libcantera_shared.so.3.2.0` 的相对符号链接；各依赖许可必须完整保留。

```sh
export HUNDUN_CANTERA_ROOT=/absolute/path/to/cantera-sdk
sha256sum "$HUNDUN_CANTERA_ROOT/lib/libcantera_shared.so.3.2.0"
# 093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760
```

配套独立运行包中的 `ct/` 就是该 SDK。包内 glibc、OpenMPI 组件与启动器
共同构成运行环境；只复制 `hundun` 或 `.so` 文件不足以在旧系统运行。

## 从源码构建反应程序

在仓库根目录运行（SDK 已取得并校验）：

```sh
cmake --preset release \
  -DCMAKE_C_COMPILER=clang-15 -DCMAKE_CXX_COMPILER=clang++-15 \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld-15 \
  -DCMAKE_SHARED_LINKER_FLAGS=-fuse-ld=lld-15
cmake --build --preset release -j 8
b3/versions/v0.4/hundun --version
b3/versions/v0.4/hundun check examples/g
```

若环境提供无版本后缀的 `clang`、`clang++`、`lld`，并且实际为上述版本，
可直接使用 README 的 `cmake --preset release`。更换编译器/MPI/SDK 后
使用新构建目录。生产占满核数时不同时执行编译或回归。

## 新算与 Restart 冒烟检查

以下为 8³ 小型 JL4/四场/Vreman/动态 TCR 示例，不是完整 GTMC 几何。
输出目录应为未使用的新目录。普通用户运行，不要为 MPI 无故使用 root。

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mpiexec -n 1 b3/versions/v0.4/hundun run examples/g \
  --output g1 --steps 2 --output-interval 0 --restart-interval 1 \
  --diagnostics-interval 1 \
  --initial-state 100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001
mpiexec -n 2 b3/versions/v0.4/hundun run examples/g \
  --output g2 --steps 1 --restart g1/Restart \
  --output-interval 0 --restart-interval 1 --diagnostics-interval 1
python3 tools/v04_evidence_validate.py runtime g1/evidence.jsonl
hf_generation=$(tr -d '\r\n' < g1/Restart/current)
python3 tools/v04_evidence_validate.py runtime g2/evidence.jsonl \
  --run-start-manifest "g1/Restart/$hf_generation/manifest.bin"
```

证据检查成功时可无输出，以退出码 0 为准。真实长测还需检查时步是否
塌缩、连续性与质量/能量/元素账本、温度及释热率；启动成功不能代替这些验证。

## 旧系统与独立运行包

Ubuntu 18.04/20.04 的宿主 glibc 通常不能直接加载参考 SDK。采用
[完整独立运行包](../pkg.md)，不更改系统库、不把私有库全局加入
`LD_LIBRARY_PATH`：

```sh
tar -xzf hf.tgz
cd hf
python3 pack.py --verify .
./run check case
./mpirun -n 2 --bind-to core ./run run case --restart restart \
  --output next --steps 1 --restart-interval 1
```

`./run` 和 `./mpirun` 显式选择同一私有加载器和 MPI 组件。包的 manifest
对应其原程序与输入；更新程序必须重新登记身份，不能继续声称原包校验覆盖
新程序。用户要求最新提交时，应在匹配环境独立重建，不直接复用旧包程序。

2026-09-16 的跨机复验使用 `8235e59` 归档、Clang 15.0.6、锁定 SDK，
在 glibc 2.31 主机通过私有运行库完成上面的 1→2 进程恢复和证据审核。
这一结果覆盖环境、原生 I/O 和小型反应组合，不是大规模实场长期验收。

## 旧 COAST 检查点不是原生 Restart

完整 GTMC 需要原网格轴、IBM 标记、21 个入口标签、热化学资产及全部
随机场。`examples/g` 不能代替实场。`tools/v04_pdf_import.cpp` 的转换
入口为 `v04_pdf_import CASE TRANSFER RESTART`；TRANSFER 是经过审核的
二进制转移资产，不是把 COAST Restart 目录改名。

先核对分区覆盖、单位质量摩尔量到质量分数的换算、组分顺序、微小负值
处理及库存变化。缺失的方法历史经显式 V1 恢复重建，再做写盘读回和原生
续算验证。参考场、物理模式和收敛标准保持明确，不能为启动方便重新点火、
丢弃随机场或关闭守恒审核。详见[接入记录](../rc.md)。
