# 控制方程

当前流动路径面向低马赫数、牛顿流体。质量和动量方程写为

\[
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf{u})=0,
\]

\[
\frac{\partial (\rho\mathbf{u})}{\partial t}
+\nabla\cdot(\rho\mathbf{u}\otimes\mathbf{u})
=-\nabla p+\nabla\cdot\boldsymbol{\tau},
\]

其中

\[
\boldsymbol{\tau}=\mu\left(\nabla\mathbf{u}+\nabla\mathbf{u}^{T}
-\frac{2}{3}(\nabla\cdot\mathbf{u})\mathbf{I}\right).
\]

密度可以是常数、随物质标量变化，或由理想气体闭合得到。低马赫数模型中的压力用于满足质量守恒，不求解声波传播。

被动或物质标量使用守恒输运形式，并带分子扩散。具体启用字段由配置决定。
