# E2a 结果

按 [预登记协议](PROTOCOL.md) 的原判据，**两槽完整流程可行的主张被推翻**：边界来源不能在只读、不构建条件下取得完整的声明级证明值依赖读数，原槽记不可行且不替补。残余槽有一个前瞻闭环，按期结算为 `refuted`。其答案计分链路已被观察到；x 与 mode 的活证明路径判形链路未被观察到。选题政策不裁决。

本结果沿用输入提交 `b64db4da10a725240f4b93adb8e8df6abd6806b5`，主仓基线仍为 `442e8d0a8855794547208a39a9ae17cd9798e572`。协议、候选框、事件、预测、阈值和终局均未改写。原 README 的前瞻框结果保留其当时范围；本文件结算其后的执行证据。

## 三类主张

| 主张 | 按协议结算 | 证据与边界 |
| --- | --- | --- |
| 流程可行 | **两槽完整可行被推翻**；仅残余槽闭环被观察到 | [compiled-availability.json](compiled-availability.json)、[decisions.json](decisions.json) 和完整边界框保留缺失读数及不可行结果。协议明确缺失来源消耗原槽；另一个单元不能补足。残余闭环证据见下表。 |
| 读数可信 | 答案计分与显式缺失处理在该单元上得到支持；不作全链路已验证声明 | 在线与离线 L2 均退出 0 且字节相同；独立 evaluator 核对事件、结算、时序、论证判形与算术。一个按期反驳提供 a=0，答案计分链路已观察。没有已编译活证明路径，x 与 mode 未测，该段链路未观察。 |
| 选题政策 | **不裁决，仅描述** | 两来源中仅一个形成事件；两预测均给真命题低概率，廉价探针内反驳，论证被评为 bind-only。单个已计分单元不提供政策效果、校准、长期价值或因果比较证据。 |

## 残余单元的前瞻闭环

载体为 [issue #18101](https://github.com/the-omega-institute/trureturing/issues/18101)。其 issue 对象与全部四条评论原始 REST 字段保存在 [issue-18101.rest.json](results/issue-18101.rest.json)，可离线读取下列全部公开证据。

| 协议义务 | 结果证据 |
| --- | --- |
| 事件先冻结且保持原文 | `decisions.json` 在 `07:14:23Z` 冻结残余第 3 项；`93d43316` 保存该事件。预登记、两席预测、终局及 L2 的 statement/quantifiers 与冻结对象逐字相等。 |
| 预测早于首个数学探针 | [预登记评论 6231250391](https://github.com/the-omega-institute/trureturing/issues/18101#issuecomment-6231250391) 的 REST `created_at=updated_at=07:32:48Z`；两席 p=0.001、0.02 均已固定。solver 报告 `execution_started_at=07:33:24Z`，晚 36 秒；原 proposer 的无探针声明保存在 `decisions.json`。实际执行时刻是工作流自报，非独立进程遥测。 |
| 固定截止与资源上限 | 预选发布秒 `07:32:47Z` 对应原截止 `08:32:47Z`。实际创建晚 1 秒，可用窗口 3,599 秒，满足协议的正值且不超过 3,600 秒规则；未编辑或重发重置。solver 自报 91 秒完成，处于最初 300 秒探针及 1,800 秒求解上限内。evaluator 的固定判词与完成时刻为 `07:38:52Z`、`07:39:49Z`；由 solver 证据发表 `07:37:10Z` 至 evaluator 证据发表 `07:40:57Z` 的外包络也只有 227 秒，小于 600 秒。 |
| solver 加载 L1 求解 | [solver 评论 6231287063](https://github.com/the-omega-institute/trureturing/issues/18101#issuecomment-6231287063) 标明独立 codex-cli flight、加载 `skills/autonomous-research/` 且未读预测；其返回结论确认读取两个技能并按已解决探针停止。来源与必要事实见 [sources.json](results/sources.json)。 |
| 独立 evaluator | [评论 6231318594](https://github.com/the-omega-institute/trureturing/issues/18101#issuecomment-6231318594) 独立核对反例、原载体/读出/目标、时序及 §3.2 论证判形。其第 1–4 项判词在读预测前固定，再单独计算 Brier；未先读 L2 计分输出。 |
| 终局早于截止 | `settled_at=07:34:55Z`，solver 证据发表 `07:37:10Z`，evaluator 证据发表 `07:40:57Z`；[终局评论 6231318725](https://github.com/the-omega-institute/trureturing/issues/18101#issuecomment-6231318725) 发表 `07:40:58Z`，比原截止早 3,109 秒。终局保持 `refuted`、`shape_evidence: null`。 |
| 角色分离及无替补 | 两预测 flight、solver `6ad324be35d500a2c9d6a7a5`、evaluator `6ad325bc5d51934501461e0e` 是四个不同席位；proposer 只生成框与事件。输入隐藏及角色独立是已归属的工作流声明，不等于先验无污染。边界仍无事件，残余仍为原第 3 项；未新增事件、替代失败槽或扩大样本。 |

全部时刻均为 UTC。既有反驳是 p=(0,2)、q=(2,0)、z=(0,1)：独立 evaluator 确认两位置 Grundy 读出同为 2，到同一 z 的合法性不同。完整定义推导已保存在 solver 评论中，本次结果整理未重新求解或做新的数学探针。该结果不反驳既有 `row_exists_bound`，也不是当前 Lean 编译、kernel 认证或新冻结交付。

## L2 读数及复现

实际在线读取在主仓只读树运行，命令为：

```sh
PATH='/Users/auricstudio/.local/share/uv/python/cpython-3.12.13-macos-aarch64-none/bin':"$PATH" make -s -C tools research-readouts ISSUES="18101" AS_OF="2026-10-17T07:41:07Z"
```

Python 3.12.13 的 PATH 仅作用于该命令；未修改宿主配置或构建主仓。实际 stdout 分别为 [readouts-online.json](results/readouts-online.json) 与 [readouts-offline.json](results/readouts-offline.json)。[data.json](results/data.json) 保存采集时点、命令参数、退出码、分页完整性与哈希，不保存命令流水。REST issue 报告 4 条评论，完整分页返回 1 页共 4 条；所有终局引用已包含在快照中，外部 evidence 数为 0。

两份输出的 SHA-256 都是：

```text
3bf407e42620a0552ea0bc10aa3a27d4385c69765a8a9e31518128088df814e5
```

| 读数 | 结果 |
| --- | --- |
| 单元数 / 按期计分数 | 1 / 1；未形成事件的边界槽不伪造成 L2 单元 |
| 答案 | a=0（refuted） |
| 两席 Brier 损失 | 1/1000000、1/2500 |
| 等席平均概率 / 等席平均损失 surprise | 21/2000 / 401/2000000 |
| 原阈值 / y / z | 9/16 / 0 / 0 |
| x / mode / evidence_status | null / null / missing-shape |
| score_coverage / mode_coverage / mode_denominator | 1 / 0 / 0 |
| five_counts | p0=0、p2=0、p5=0、p25=0、p3=0 |
| X / Y / Z | 全部 null（完整模式分母为 0） |
| missing_counts | prediction_units=0、answer_units=0、shape_units=1、mode_units=1 |

这些分母描述唯一形成的单元，不把一槽计分覆盖率 1 写成原两槽完成率。平均损失是 `(1/1000000 + 1/2500)/2 = 401/2000000 < 9/16`，与 evaluator 的独立算术完全一致；它不是平均概率的 Brier 损失，也不是数学价值读数。

从指定主仓只读树重放（省略 ISSUES，并向 stdin 提供快照）：

```sh
PATH='/Users/auricstudio/.local/share/uv/python/cpython-3.12.13-macos-aarch64-none/bin':"$PATH" make -s -C tools research-readouts AS_OF="2026-10-17T07:41:07Z" < /Users/auricstudio/trureturing-experiments-0addbbde-e2a/docs/reports/autonomous-research-e2a/results/issue-18101.rest.json
```

必要核验程序 [check_results.py](results/check_results.py) 以同一快照重放钉版 L2，比较字节、算术、事件文本、时间约束及受保护文件差异；不计算事件的数学答案。运行方式（实验仓根目录）：

```sh
python3.12 -B docs/reports/autonomous-research-e2a/results/check_results.py --main-repo '/Users/auricstudio/Library/Application Support/heca/worktrees/568cd92c/lively-pebble'
```

核验结果在 [verification.json](results/verification.json)。离线输入保留 REST 原字段，但离线模式本身不能认证 GitHub 元数据或发现删去的评论；归属与完整性来自本次在线分页获取及报告条数核对。

## 发现与边界

| 现象 | 证据位置 | 边界 |
| --- | --- | --- |
| 边界来源需要编译后的声明级证明值依赖读数；现有 harness 在本次只读、不构建条件下不提供 | `compiled-availability.json`、`boundary-frame.json.gz`；钉版 `CompiledAxioms.lean:9` 区分 type 与 value，但没有完整已编译 D5 环境/依赖图 | 未测用户数保持 null；未用 imports、源码名称命中或头部声称替代。仅说明本次输入与约束下不可取得，不声称所有环境都不能测。 |
| 现行契约下，本内部候选的全 bind-only 结算按 CLAUDE.md §3.2 不交付，因此其 x 结构上不可测 | solver 与 evaluator 评论；钉版 `CLAUDE.md:198–206`；L2 契约的 Terminal block and evidence；终局 `shape_evidence=null` | evaluator 判的是结算论证，没有已编译声明活路径；L2 所需来源不存在。不得据此填 x=0 或 mode=p0。该候选无真实辅助声明消费者，也非外部具名开放问题例外；不推广到所有 bind-only 情形。 |
| 两预测席自报同属 GPT-6 族，未取得多族 | 预登记两条 `model_family: GPT-6`，L2 `model_families=["GPT-6"]`；`sources.json` 的 role_evidence | 载体分别为 codex-cli 与 nyxid-oracle，不能据载体差异声称模型族多样性；自报标签不是供应商独立认证，也不能证明先验独立。 |
| 第四轮面板的 oracle 席两次未完成，改派同族 codex | `sources.json` 的 panel_four，明确归属于 caller 实施 brief 验收项 4(d) 所提供的观察 | 未独立核验该面板的载体遥测，不推断失败原因或质量影响；此处指面板 oracle，不是已有预登记预测的残余单元 oracle。不保留日志或回执副本。 |

没有新增机制、配额建议、价值函数或政策裁决。残余充分性、六坐标效果、长期研究价值及模型校准均未由本试点测得。两槽可行性失败与答案计分读数可信是不同的结算，不相互抵消。

## 交付与保留范围

本次仅新增结果正文、必要数据与核验程序，并在 README 追加入口。`git diff --exit-code b64db4da10a725240f4b93adb8e8df6abd6806b5 --` 加上 `verification.json` 中逐项列出的 9 个受保护路径返回 0；README 原字节保留为前缀。结果提交不包含日志、实施日记、命令流水或 runner 回执副本。issue 快照和两份 L2 stdout 是本任务明确要求保留的复现数据。

原 README 已披露的存储边界仍未解决：指定仓库 REST metadata 返回 `private=false`，因此不能声称成果已存入私有仓库。仓库身份、当前分支提交与推送由 caller 明确授权；未改变可见性设置。推送只更新既有 PR #87，不合并。
