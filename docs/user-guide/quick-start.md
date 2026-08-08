# 快速开始

以下命令使用随源码提供的最小模板。它用于确认安装、输入解析和基本运行路径，不代表精度或收敛性验证。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 2
build/release/src/hundun examples/minimal/case.json --validate
build/release/src/hundun examples/minimal/case.json --print-resolved
mpiexec -n 1 build/release/src/hundun examples/minimal/case.json
```

相对路径以 `case.json` 所在目录为基准。运行前先用 `--validate` 检查结构，再用 `--print-resolved` 查看规范化后的实际配置。

模板的输出位于 `examples/minimal/output` 和 `examples/minimal/Restart`。重复运行前，请移动或清理上一次运行生成的目录；程序不会把不完整或不匹配的 Restart 当成有效状态。
