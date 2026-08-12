# 当前版本能力

当前源码是 `0.2.0 candidate`。命令行运行路径覆盖结构化网格、三种密度闭合、五类外边界、固定与自适应时间推进、静止 STL IBM、WALE、并行 Checkpoint 和结构化诊断。

profile-1、profile-2、profile-3 是 body-fitted density profiles；profile-4、profile-5、profile-6 加 static IBM；profile-7、profile-8、profile-9 再分别组合 WALE。浸入路径采用两次 PISO pressure corrector，并统一压力、算子、最终通量和力诊断的数据权威。Checkpoint v3 presence 1--9 与这张表一致。

九个 profiles 的实现和小规模测试为 `implemented-and-accepted`。正式 0.2.0 发布仍取决于 S3-V1；生产级 LES、反应流、喷雾、颗粒、移动壁面、GPU、大网格性能和任意复杂几何精度不在当前能力声明内。
