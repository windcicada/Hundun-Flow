# 安装要求

HUNDUN-FLOW 直接从源码构建，不需要安装脚本或在线下载依赖。

必需组件：

- CMake 3.21 或更新版本；
- 支持 C++17 的 C/C++ 编译器；
- MPI 3 实现及其 C++ 开发文件；
- POSIX 线程库。

在 Ubuntu 或 Debian 上，安装 `cmake`、`g++` 和一种 MPI 实现即可。MPI 的运行器与编译时库必须来自同一套安装。若系统同时装有多套 MPI，请先确认 `mpicxx`、`mpiexec` 和 CMake 检测到的库一致。

源码已经包含 yyjson。构建时不访问网络，也不需要 Python。
