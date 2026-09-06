# HGSOC CORNETO 科学有效性审计与修复：2026-09-06

## 结论先行

中心科学问题合理，但当前结果**尚未回答**“跨 cohort 可复现、对患者重复和模型选择稳健、相对 reference 富集、且有外部机制支持”的完整命题。RNA/NMF 与 regulatory 分析保留可用的探索性证据；metabolic b25 长作业目前主要在求解一个弱 OCM-specific、未校准 medium 的 sparse-network benchmark，不能当作患者特异代谢结果。

本次最重要的新结果不是 gap 下降，而是一个带阳性对照的反证：**20 份已保存的 independent attempts（17 个不同 OCM/run）中，仅保留 indicator 选中反应时，20/20 无法达到原定 growth floor；允许保存通量中实际使用的反应后，20/20 可行。** 因而“这些 indicator sets 已是可用的生长支持网络”被本次检验否定。它不否定 CORNETO 方法、HGSOC 研究问题或这些细胞真实的生长能力。

已修复 numerical/receipt 验收和 joint gate，部署到 Roihu；四个冻结 b25 contexts 已设置启动前 scientific review hold。没有取消健康 RUNNING 任务，没有放宽 MIPGap，没有改变原 context 或 medium。新工作应先解决输入与数值有效性，再决定是否值得启动全集长求解。

## 证据、分母与可复核方法

本次 read-only audit 于 **2026-09-06 17:03:54 UTC** 读取服务器文件；部署验证于 **17:14:52 UTC（18:14:52 BST）** 完成。Slurm controller 当天间歇性不可用，不能把早前队列截图描述成此刻的实时状态。

| 证据对象 | 文件 | 数字的来源或计算规则 |
|---|---|---|
| 60 个冻结 contexts、20 份 attempts、45 份 regulatory receipts | [完整审计 JSON](../evidence/scientific_validity_audit_20260906.json) | `cohorts`、`attempts`、`regulatory_grid`；每项保留源文件 SHA256 |
| Fixed-indicator LP 与阳性对照 | [受控 LP 审计](../evidence/fixed_indicator_lp_controlled_audit_20260906.json) | `results` 的三个 LP status；一行一份 attempt，不是一个独立患者 |
| 原 NMF/regulatory 数值 | [历史 receipt 快照](../evidence/roihu_result_snapshot.json) | 本轮重新计算快照引用的 14 个远程源文件 SHA256，14/14 匹配；不是只信任旧文档 |
| OCM 与 patient 分母 | [逐样本 registry](../evidence/study_ocm_registry.tsv) | 筛选 `primary_cohort_eligible=true`，按 `patient_id` 去重；本轮 registry SHA256 已记录 |
| 服务器修复与 hold | [部署验证](../evidence/scientific_audit_deployment_20260906.json) | 代码 SHA256、4 个 context hash、startup gate 和未导入 Gurobi 的检查 |

`scripts/audit_metabolic_scientific_validity.py` 不调用 solver。它的 cross-context 检查基于保存的 sparse flux summary，不能证明完整 feasible sets 相同。其“尚未做 fixed-indicator LP”的 limitation 仅指该文件自身；随后完成的独立 LP 审计补上了此检验。

`scripts/audit_fixed_indicator_feasibility.py` 使用 SciPy HiGHS，不占 Gurobi session。对每个保存的 selection：

1. 将未选反应严格设为零；保留 Human-GEM stoichiometry、原 model bounds 与该 OCM 的 frozen overrides。
2. 检查原 90% growth floor；另行检查能否达到保存的 biomass flux。
3. 阳性对照允许 `active_by_indicator ∪ active_by_flux` 中的反应，其余仍严格为零。
4. 每个 LP 最长 10 秒，primal/dual feasibility tolerance 为 `1e-8`。20/20 对照为 `optimal`，而两种 indicator-only 检验均 20/20 `infeasible`，不是 time limit。

上述诊断不改写原 receipts。重现命令如下；`--root` 指含这些确切 context/attempt 文件的项目目录，输出文件必须尚不存在：

```bash
PYTHONPATH=scripts python scripts/audit_metabolic_scientific_validity.py \
  --root PROJECT_ROOT --human-gem Human-GEM-v1.4.1.xml \
  --registry evidence/study_ocm_registry.tsv \
  --snapshot evidence/roihu_result_snapshot.json

PYTHONPATH=scripts python scripts/audit_fixed_indicator_feasibility.py \
  --root PROJECT_ROOT --human-gem Human-GEM-v1.4.1.xml \
  --output NEW_CONTROLLED_AUDIT.json
```

## 目前究竟在求解什么

Human-GEM v1.4.1 含 13,096 reactions、3,628 genes。四队列分别从表达映射中选择最多 25 个 candidate reactions（b25），不是把全部表达信息都转成约束；每个 OCM 实际生效的 reaction caps 更少。RNA 输入为 log1p(TPM)，expression cap 的生成和缺失处理来自 frozen context 对应代码，不代表已经验证的 enzyme-capacity 测量。

Independent objective 为 `-biomass + 0.1 × sum(reaction indicators)`，另有该 context 最大 growth 的 90% 下限。Cohort joint jobs 使用 union regularization、nominal lambda=1。内部 TimeLimit 为 252,000 秒（70 h），Slurm 为 72 h，MIPGap 为 `1e-4`（0.01%），8 threads。已保存的 19 个 70 h attempts 最终 relative gap 为 **2.9679%–4.1531%**；另有 1 个 600 秒 smoke，gap 为 5.8292%。这些是 20 次 attempts，不是 20 个不同 OCM；共涉及 17 个不同 run。

Gap 是目标函数 incumbent 与 bound 的距离，不是生物学误差率或发现可信度。重跑同一 deterministic context 又不读取 `.mst`，不等于续接之前的 search tree；保存 `.mst` 本身也不证明实现了 warm start。

## Metabolic 结果为什么目前不能作生物学解释

### 输入区分度不足

| Study | Primary OCMs | 实际 expression caps：min / median / max | 无 expression caps 的 OCMs | Canonical independent files |
|---|---:|---:|---:|---:|
| E-MTAB-7223 | 9 | 1 / 15 / 17 | 0 | 0 |
| E-MTAB-10801 | 13 | 0 / 9 / 15 | 1 | 0 |
| E-MTAB-11000 | 11 | 0 / 13 / 22 | 1 | 0 |
| E-MTAB-14568 | 27 | 0 / 8 / 18 | 4 | 0 |
| Total | 60 | 不将不同 cohort 的 median 相加 | 6 | 0 |

这里的 cap count 排除了 biomass constraint。所有 60 个 context 的 growth optimum 都为 **187.35362997658078**。20 份保存的 flux summary 中，所有实际施加 expression cap 的反应均无报告非零通量（阈值 `1e-7`）。逐一将保存的通量放到其他 context 的 bounds 中，分别有 **50–60/60** 个 context 可接受；其中 **15/20** 份兼容全部 60 个 context 的 bounds。

这削弱 OCM specificity，但不等同于证明所有可行域相同。个别 non-binding constraint 仍可能排除其他解。更强的可证伪检验应比较真实表达、expression-null 和合适的 shuffled-expression controls 对结果的影响，而不是把“每个样本都输出一个文件”当作 sample specificity。

### Generic exchange bounds 允许未经验证的能量底物摄取

20 份 summaries 均报告 ATP 与 phosphocreatine exchange flux 约为 -1000；PEP uptake 也接近开放边界，glucose exchange 为零。这是模型允许的可行状态，不能解释为 OCM 的实测营养偏好，更不能据此说 HGSOC 不依赖 glucose。当前没有用实际 OCM culture medium 或 measured uptake 校准这些 exchanges。

公开 Meeson 研究把 medium 与 proliferation information 纳入 transcriptome-driven modelling；其 cell-line 结果不能替代本项目 OCM medium 的证明。其发表版包含 TPI1 机制实验，并给出 E-MTAB-16770 作为 perturbation RNA-seq 来源；本轮仅核对论文，不声称已完成该数据集 intake 或外部验证。[Meeson et al., 2026](https://link.springer.com/article/10.1186/s40170-026-00425-6)

### Reaction selection 与 flux 不一致，并非只有 gap 未达标

每份保存的 solution 有 **9–37** 个 reaction：indicator selection 不包含它，但保存的 flux 非零。Sparse-summary 最大 mass-balance residual 为约 `4.18e-8`，因此“mass balance 几乎通过”不能证明取整后的 selection 可用。

此前检查 `.mst` 已见小于 integrality tolerance 的 binary 数值通过大 reaction bound 放大为非零 flux；本轮 LP 进一步直接检验 selection，而非只凭这个机制推断。Gurobi 官方文档说明大系数与 integrality tolerance 可产生 trickle flow，`IntegralityFocus` 可缓解但不是保证。应测试 tighter justified bounds、numerical settings 与 fixed-indicator feasibility，不能仅降低 active-flux threshold 或把 partial 改名 completed。[Gurobi IntegralityFocus](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#integralityfocus)

**因此，这些 selections 不应进入 knockout/FVA。** 可保留完整 incumbent 作为优化诊断，但不把其反应集合解释为可生长的 OCM-specific network。这个结论依赖本次 context、阈值和保存的 selection，不外推到全部 CORNETO 实现。

## 仍然有效的结果，以及实际回答的问题

| 已核对结果 | 可以回答 | 不能据此声称 |
|---|---|---|
| 117 RNA runs；60 primary OCMs / 52 patients | 数据分母和跨 study 对应关系可追踪 | 60 个独立患者或所有 reference 可混为一组 |
| Pooled rank-3 NMF clusters 为 23/20/17；cophenetic 0.939、silhouette 0.770 | 该预处理/拟合方案下有稳定性可量化的表达分组 | 已发现三个临床或代谢亚型 |
| Rank 2 对应 0.972/0.899，更高于 rank 3 | Rank 3 不是唯一由这两项稳定性指标支持的选择 | Rank 2 必然是真实 biological dimension |
| Pooled-vs-cohort ARI 0.377–0.723 | 同数据不同 fitting scope 的一致性并不完全 | 独立的跨队列 replication，因为这些 cohort 参与了 pooled fit |
| Patient-balanced NMF ARI 0.842；regulatory union Jaccard 0.890 | 一次 response-blind 去重复选择具有一定稳定性 | 对所有患者代表选择均稳健 |
| Regulatory grid 45/45 solver status 为 optimal | 参数扫描求解完成 | 45 个非空、稳定、独立生物学发现；其中 27 个 union 为空 |
| Richer graph-policy pooled union Jaccard 0.202、mean sample Jaccard 0.108 | 结论对 bundled graph policy 显著敏感 | 已定位是 PKN breadth 单因素造成；这次还同时改变 inputs、outputs、depth |
| TPI1 model gate 找到 ENSG00000111669 / HMR_4391 | 模型标识与 GPR 定位可用 | 已证明 TPI1 是 OCM-specific essential gene |

Regulatory 方阵导出 bug 已修复。45 份当前网格的 processed edge count 都不同于 condition count，故这批文件没有触发该具体 bug；不由此宣称所有历史 regulatory 分析均已完整重审。较大 lambda 下 empty network 是当前 objective scaling 的结果，不是“没有调控生物学”。

研究主张还缺少三类独立证据：held-out cohort 的稳定投影、labelled reference 的同框架对比、与训练 RNA 无关的机制/表型验证。NMF 与 regulatory 来自同一 RNA 时的一致性只是 internal consistency。Study association 的 p=0.207 也不是排除了 batch effect。

## 已修复与实际部署状态

| 问题 | 修复 | 验证和限制 |
|---|---|---|
| `optimal_inaccurate` 可能被当作 canonical success | 要求 CVXPY optimal、Gurobi OPTIMAL、有 incumbent、有限 objective/bound/gap 且符合原 gap contract | 不放宽 `MIPGap=1e-4`；partial 仍保留 |
| Receipt skip/assembly 只核 status/hash | 核对 exact conditions、primal metrics、flux/indicator agreement、`.sol/.mst/.gurobi.log` SHA256；assembly 记录输入 receipt SHA256 | Numerical acceptance 仍不是 biological validation |
| Joint 只有 scheduler dependency，没有完整内部 prerequisite gate | 加载 solver 前逐一验证全部 cohort independent receipts | 修复 subset retry 成功就可能放行 joint 的漏洞 |
| Regulatory 任何 incumbent 都可标 completed；方阵可能误转置 | 非 optimal incumbent 标 partial；明确优先使用 edge×condition orientation；拒绝 nonfinite matrix | 当前 45 份 grid 的 solver statuses 和非方阵条件已复核 |
| 11000 被错误串联到其他 cohort 的 scientific success | 改回 `afterany:1083050_*`, `afterany:1083051_*`, `afterany:863034_*` | 当次 `squeue` 已确认；此后 controller 再次超时 |
| 无效 b25 长作业继续排队执行 | 四个 `checkpoint_b25/scientific_review_hold.json` 加应用级 startup gate | 4/4 远程 `_context` 检查阻断在 solver import 之前；不是已确认的 Slurm hold |

29 项 targeted tests 全部通过，包括 telemetry 伪阳性、artifact 篡改、缺失 cohort member、完整 primal trickle flow、最终 receipt promotion、regulatory 方阵及 LP 阳性/阴性对照。服务器 Python 3.9 syntax checks 通过。部署前五个原脚本哈希与 GitHub 基线一致，备份保存在项目内 `.codex_stage/scientific_audit_20260906_before/`。

最后一个成功的 targeted queue snapshot 显示 9 个 solver tasks：`937737_3/4/6`、`948765_3/4/5`、`834322_11/12/13`。后续 live query 超时，不能保证此刻仍是 9 个。没有取消它们。启动保护只阻断以后加载 runner 的任务，不会给已经运行的旧进程注入新代码；其最终 receipts 必须重新审计。

`863034`、`834323`、`1083050/1083051` 及现有 joint/assembly/comparison/TPI1 链仍保留历史 lineage。若被启动保护主动阻断，应标 `scientific_review_hold`，不是 OOM、license failure 或待自动修复的 infrastructure failure。**不得自动为这些退出重新提交 70 h retry。**

### 同日23:02 BST 的继续检查

用户报告 usage window 重置后，低频 targeted `squeue` 已恢复，仍显示上述9个任务 RUNNING。
9份日志都在更新；约运行29–32 h，live gaps 分别为：7223 tasks3/4/5 的
4.05/3.97/3.87%，10801 tasks3/4/6 的3.14/3.65/3.62%，14568 tasks11/12/13 的
3.10/3.97/3.17%。[日志尾部证据](../evidence/scientific_audit_running_logs_20260906.json)
包含时间与完整末行，不能当成最终 receipt 或 ETA。

进一步给 pending solver jobs 设置可撤销 Slurm hold 的调用在实际执行前被安全审核拒绝。
因此**本轮没有施加调度器级 hold，也没有绕过拒绝**；如果需要这一额外的排队资源控制，
须用户明确授权。此前已部署的应用级 startup gate 继续生效，运行中的进程未改动。

## 接下来真正值得推进的路线

| 优先级 | 可证伪问题与动作 | 放行条件 / 并行关系 |
|---|---|---|
| P0：已执行 | 保存的 indicator selection 是否真的能支持 growth？使用 fixed-indicator LP 和 positive control；修复所有下游验收 | 当前 20/20 selection 失败，故保留解释 hold |
| P1：输入重建 | Actual culture medium、exchange direction、GPR AND/OR、missing vs true-zero expression、growth calibration 是否有来源？建立 versioned context，保留 b25 为历史 benchmark | 不猜测未知 medium；不盲目把全部零表达设 knockout |
| P1：判别力 pilot | 真实表达是否比 null/shuffled inputs 增加 sample-specific information？在同一 model/media/candidate policy 下比较，可加入不依赖 sparse indicators 的连续 LP 基线 | 先短小 pilot，检查输入扰动效应和 full-primal/fixed-selection feasibility，再申请全集长作业 |
| P1：RNA/regulatory 验证，可并行 | 一次去重复稳定能否扩展到重复 patient-balanced draws？未参与拟合的 cohort 能否重现状态？ | Patient-grouped bootstrap/permutation；训练集内选基因/fit NMF，held-out projection；OCM74 跨 study，必须按 patient 隔离 |
| P2：分离 modelling 因子 | Candidate policy、media、lambda scaling、PKN breadth、input/output/depth 各自影响多大？ | 一次改变一个因素；independent lambda=0.1 与 joint lambda=1、不同 condition 数的目标需可比化；union regularization 自带 reuse 偏好 |
| P2：Reference 对比 | Tumour vs stroma/non-HGSOC 的差异是否超过 patient/cohort/culture effects？ | 同一 preprocessing/medium/model policy；匹配或分层患者；cell-line controls 仅描述；报告 effect sizes 与不确定性，不把 absence 当特异性 |
| P2：外部机制 | TPI1/Meeson 的独立 perturbation 能否支持方向与上下文特异性？ | 先 intake 原论文公开 supplementary/perturbation data，后 WT/KO/FVA；模型内 deletion 不等于实验机制支持 |

同患者分析不是样本“太少不能做”，但有效 repeated-patient denominator 为 **7**：6 个两-OCM family、1 个三-OCM family，共 15 个 OCM。可作 paired descriptive、matched between-patient distance 和 family-level uncertainty；8 个 baseline-to-other contrasts 不应当作 8 个独立患者。编号后缀本身不证明时间顺序、暴露或 acquired resistance。

CORNETO 有 expression-context-specific 和 multi-condition 建模实现，但方法名本身不能保证本项目的 candidate/medium/objective 合理。比较方法必须锁定项目实际依赖版本，不把当前在线示例直接替换为已部署 API。[CORNETO context-specific tutorial](https://corneto.org/stable/tutorials/fba/context-specific-metabolic-omics.html)

本轮未自动启动新版本的生物学全集求解，因为 actual medium/phenotype 校准尚未得到本次验证。优先路线已经明确；不再以“还有排队任务”作为科学进展指标。

## 定时恢复与交付说明

依用户要求，已把本会话的 existing heartbeat 临时设为 **2026-09-06 16:32 Europe/London** 的一次性恢复。工具侧恢复信号实际在 **16:56:09 UTC（17:56:09 BST）** 到达；不能声称实际在 16:32 启动，延迟原因未核实。收到信号后本会话继续完成审计。

一次性恢复后，已恢复此前实际配置的 **每周四 09:00 Europe/London** heartbeat，并更新 prompt：尊重 scientific review hold，不做同参数 70 h 自动重试，保留 RUNNING 证据、只报告有意义变化；没有 model override、没有另开任务。旧状态文档的“每 30 分钟”与实际配置不一致，本轮一并纠正。

**后续同日更新：** 用户报告5小时 usage window 已重置并要求04:00继续，因此 existing heartbeat 又临时调整为 **2026-09-07 04:00 Europe/London（03:00 UTC）** 的一次性恢复，配置已核实。04:00 的 prompt 要求先核对哪些工作已经交付，只继续未完成部分，随后恢复每周四09:00；尚不能宣称未来实际触发时间。

续办交付时还发现 Roihu 的两份状态文件停留在2026-08-27，落后于 GitHub delivery branch。比较表明是缺少后续历史记录，而非新的未合并分析；这次同步需保留旧文档备份并补齐最新审计，避免不同执行环境读到不同的继续规则。

GitHub delivery 使用隔离目录 `/private/tmp/hgsoc-scientific-audit-20260906`，避免 OneDrive I/O 阻塞。原 OneDrive working tree 未宣称已同步。本文与两个状态登记是本次审计来源；`main.tex` 仍是标注日期为 2026-08-12 的历史稿件，不能把其旧快照当成本次最终审计。

**23:09 BST 交付状态：** 本地 implementation/audit commit 为 `fa7efa7`。安全审核分别拒绝了本次21文件向 `xm2325/hgsoc_corneto` delivery branch 的 push，以及四份 audit JSON 回传原 Roihu evidence 目录；均未绕过，GitHub 仍未收到本次提交。代码与中英文文档的服务器部署已完成，但远程报告中的四个新 evidence 链接尚缺相应 JSON 副本，本地相对链接完整。需要用户对具体文件内容和目的地明确授权后才能完成这些剩余传输。

04:00 接续应先检查本地隔离目录的 commits 与本报告，不重新执行已完成的29项修复测试或重新启动 b25 长作业。若用户授权了剩余上传，再定位可用的 GitHub 认证工具并校验远程父提交；此前临时 GH executable 路径在23:09检查时已不存在，不能把已通过 `ls-remote` 当成 push 认证可用的证明。Slurm pending-job hold 也须单独明确授权，不能把后续 heartbeat 当作新授权。
