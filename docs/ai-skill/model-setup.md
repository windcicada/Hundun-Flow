# 模型与配置准备

从[当前模板](../../examples/minimal/README.md)和[input schema](../../versions/v0.4/docs/input-schema.md)准备输入。当前单相低马赫路径使用 JSON schema 1，不按历史 density/IBM/WALE 九个 profile 选择程序。

明确网格、物性数据、六面边界、离散与时间控制，以及可选标量、LES、IBM。均匀初场用 `--initial-state`；恢复使用 `--restart`，不能同时给出。

有组分时同时核对组成、物性和 EOS；有被动标量时保留实际允许的范围，不擅自限制为质量分数。冻结用户给定的物性、边界和误差标准，发现矛盾应报告，不得以经验值悄悄替换。

校验通过后仍须做对应数值验证，且不得与冻结长测争抢资源。
