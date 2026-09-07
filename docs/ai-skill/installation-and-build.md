# 自动化安装与构建

探测 CMake、C/C++ 编译器和同一套 MPI 后，按[构建指南](../user-guide/build.md)使用独立目录。可执行文件为 `build/release/versions/v0.4/hundun`，不是旧路径 `build/release/src/hundun`。

记录完整源码提交、工作区状态、构建命令、退出状态、二进制 SHA-256、build manifest、编译器及 MPI 版本。构建失败不得继续冒用旧程序；不同候选的性能证据不能混合。

测试与长测共用硬件时遵守已有作业的核数限制。不要自动停止、重编译或替换已冻结的长测程序。
