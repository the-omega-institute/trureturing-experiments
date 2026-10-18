# 理论配套实验

本目录保存 [trureturing](https://github.com/the-omega-institute/trureturing/tree/4128c47e530bc03508c9b11d1edb59a9b5e33f33) 理论正文配套的有限核验程序、结果数据、文献来源与节点注册表。理论正文和 Lean 证明位于主仓库；这里的有限枚举与数值检查不构成全称证明或 Lean 证书。

| 目录 | 内容与适用范围 | 主仓库正文 |
| --- | --- | --- |
| `cfmp/` | 度十星形的组合与有理端点检查；九周期盒的端点检查；Coxeter 有限局部检查、有限群阶枚举和非正规指数 72 作用检查 | [度十星形](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/CFMP_DEGREE10_STAR.md)、[九周期](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/CFMP_CYCLE_DEGREE9.md)、[有限商](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/CFMP_CYCLE_DEGREE9_FINITE_QUOTIENT.md)、[非正规指数 72](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/CFMP_CYCLE_DEGREE9_NONNORMAL_INDEX72.md) |
| `grid-parity-transport/` | S4、关联图、格点步词和有限矩形端口路由的结果数据；核验程序仍位于主仓库的 [verify_parity_transport.py](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/tools/scripts/agent/grid-parity-transport/verify_parity_transport.py) | [矩形网格奇偶输运](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/RECTANGULAR_GRID_PARITY_TRANSPORT.md) |
| `rro-sft-profile-closure/` | 第 29–32 节的精确整数与有理秩检查、有限路径恢复和结果数据 | [SFT 理论](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_SFT.md) |
| `rro-sft-29-35/` | 第 29–35 节的有限分拆图、群代数、标号路径和固定 C2 端点检查；详见其 [报告](rro-sft-29-35/README.md) | [SFT 理论](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_SFT.md) |
| `periodic-tree/` | 141 行节点注册表，包含地址、度量、组合标签与流向字段 | [周期树](https://github.com/the-omega-institute/trureturing/blob/4128c47e530bc03508c9b11d1edb59a9b5e33f33/docs/develop/theory/PERIODIC_TREE.md) |

## 运行入口

CFMP Python 程序仅依赖标准库。`audit_degree10_bounds.py` 按当前工作目录读取 `verify_degree10_star.py`，须从 `experiments/theory/cfmp/` 运行。其余 Python 程序可在同目录按文件名运行。不要使用 `python -O`，检查依赖 `assert`。

`verify_cycle9_coxeter_order.cpp` 使用 GNU C++ 的 `<bits/stdc++.h>`，需 GCC 和 C++17。该程序枚举 11,985,408 个矩阵，编译与执行入口为：

```sh
cd experiments/theory/cfmp
g++ -O3 -std=c++17 verify_cycle9_coxeter_order.cpp -o verify_cycle9_coxeter_order
./verify_cycle9_coxeter_order
```

SFT profile 检查需要 Python 3、NumPy 和 SymPy；运行 `python3 experiments/theory/rro-sft-profile-closure/verify_closure.py` 会在程序目录生成 `verification-summary.json`。第 29–35 节检查要求 Python 3.10 或更新版本，仅依赖标准库，其运行命令见配套报告。

## 来源与许可证

程序与数据的来源路径为主仓库同一修订下的 `docs/develop/theory/`，周期树注册表的来源为 `docs/develop/archive/theory/PERIODIC_TREE_registry.jsonl`。本目录程序和数据适用仓库根目录的 [Apache-2.0 许可证](../../LICENSE)；源仓库的许可证与本仓库许可证字节相同。文献的作者、版本、引用地址和查阅范围保存在 `rro-sft-29-35/literature_ledger.json`。NumPy 与 SymPy 是运行时依赖，本目录不包含其源码。
