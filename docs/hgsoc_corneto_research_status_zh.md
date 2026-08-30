# HGSOC CORNETO 研究状态与依赖登记（中文对应版）

最后运行更新：2026-08-30 10:03 BST（12:03 EEST）。本文件与
`docs/hgsoc_corneto_research_status.md` 对应，记录研究范围、已完成证据、排队分析、失败尝试、依赖关系与可声明范围。
仅有 Slurm `COMPLETED` 不足以证明科学分析完成；只有输出 `receipt` 通过相应内容验证后，结果才算科学上完成。

## 核心科学问题（Central scientific question）

> HGSOC OCM 中是否存在跨 cohort 可复现、对重复患者和 modelling choices
> 稳健的 metabolic 与 regulatory states；这些状态相对于 stroma/reference
> samples 是否富集，并能否得到外部 mechanistic evidence 支持？

主要证据单位是 60 个 primary HGSOC tumour OCM。它们来自 52 位患者和四个
study，因此所有分析都必须区分 OCM-level、patient-balanced、cohort-stratified
和 pooled 结果。其余公开 RNA runs 是 reference/QC 数据，不得静默合并进 primary cohort。

## 冻结的样本范围（Frozen sample universe）

| Study | 全部 RNA runs | Primary HGSOC tumour OCM | 主要 reference runs |
|---|---:|---:|---|
| E-MTAB-7223 | 36 | 9 | 16 stroma、2 cell-line controls、non-HGSOC/ambiguous 及其他排除项 |
| E-MTAB-10801 | 36 | 13 | 17 stroma、6 non-HGSOC |
| E-MTAB-11000 | 12 | 11 | 1 个 non-primary tumour run |
| E-MTAB-14568 | 33 | 27 | 6 non-HGSOC |
| **合计** | **117** | **60 OCM / 52 patients** | **57 reference/excluded runs** |

reference categories 具有不同科学角色，必须分开分析：stroma 是最强的正式
secondary contrast；confirmed non-HGSOC 是 exploratory histology contrast；
ambiguous samples 不是 controls；两个 cell-line controls 只能用于 descriptive
QC；later-passage 与 non-Tighe samples 用于 stability 或 response-blind replication，
不能用于 population-level inference。

## 分析层级与 claim limits

| Evidence tier | 分析 | 允许的解释 |
|---|---|---|
| Primary | 四个 cohort-specific primary-HGSOC metabolic/regulatory models；cross-cohort recurrence | 可复现的 model-predicted metabolic/regulatory state |
| Robustness | pooled-60、patient-balanced-52、lambda、PKN、NMF rank、alternative optima | 对 modelling choices 和重复患者的敏感性 |
| Reference | tumour vs stroma；confirmed HGSOC vs confirmed non-HGSOC；cell-line/passage QC | enrichment/shared core，不能声称 absolute presence/absence |
| Mechanistic | Meeson/TPI1、WT vs delta-TPI1、FVA、order sensitivity | 与模型/已知机制的一致性，不是 drug-response validation |
| Phenotype-linked | exact AUC、raw GI50、cumulative exposure、grouped CV | **在 exact phenotype tables 通过 intake QC 前保持 blocked** |

所有 expression-derived flux 都是模型预测的 feasible flux states，不是 measured
flux。由同一 RNA profiles 得到的 regulatory/NMF/metabolic agreement 是 internal
consistency，不是 independent validation。

**解释限制（Interpretation hold，2026-08-28 登记）：** 已保存的 b25 partial solutions
暴露出 sample-specific expression constraints 过少、exchange bounds 未校准，以及
flux/indicator 数值不一致。即使以后 solver 收敛，也不能自动解决这些问题。冻结的 b25
runs 仍属于 optimization benchmarks；patient-specific metabolic claims 需要单独的
input 与 numerical-quality audit。本次 monitor 没有修改任何 scientific parameter。

### 同患者覆盖范围（Within-patient coverage，2026-08-27 核实）

对 `evidence/study_ocm_registry.tsv` 筛选 `primary_cohort_eligible=true`，再按
`patient_id` 分组，得到七位 repeated patients、15 个 primary OCM：六组各两个 OCM，
一组有三个 OCM。加上45位 singleton patients，恰好对应60个 OCM、52位患者。
每组取一个 baseline 与其余 OCM 比较，共八个 contrasts，不是八位独立患者。

| Patient | Primary OCM IDs | 已提交的 independent array tasks |
|---|---|---|
| OCM66 | OCM66-1; OCM66-5 | 834320_5; 834320_4 |
| OCM74 | OCM74-1; OCM74-3; OCM74-5 | 834320_7; 834320_6; 834322_26 |
| OCM110 | OCM110-1; OCM110-9 | 834321_2; 834321_1 |
| OCM288 | OCM288-4; OCM288-7 | 834322_4; 834322_5 |
| OCM296 | OCM296-3; OCM296-5 | 834322_8; 834322_9 |
| OCM327 | OCM327-1; OCM327-3 | 834322_15; 834322_16 |
| OCM333 | OCM333-1; OCM333-3 | 834322_18; 834322_19 |

14568 的相应 indices 也由现有863034 repair 覆盖。这些是 independent sample solves，
不是本次新提交的 patient-level joint model。834324-834327 仍是 cohort-level joint。
已重新读取 regulatory job 592019 的 `regulatory_longitudinal_joint_l0p001.json`：
它记录这七个 families 的八个 response-blind comparisons，不是 treatment causality。
其中的 baseline labels 本身不能独立证明 collection chronology。

Patient-level joint sensitivity 已讨论，但尚未启动。检验 within-patient similarity 需要
independent solutions 与匹配的 between-patient comparisons；joint union regularization
本身就偏好 reaction reuse，不能单独证明 patient effect。OCM74 跨7223与14568，
cohort-specific candidate selection 是额外 confound。时间顺序、治疗以及
same-biopsy/spatial relationships 必须依据 source metadata，不能只按编号后缀推断。

## 已完成且可审计的结果（Completed and auditable results）

### RNA 与 pooled-input gates

- 四个 RNA studies 均已下载、定量并聚合，且通过严格 run/receipt 检查。
- Pooled primary expression 包含 60 runs、60 OCM、52 patients 和 60,609 genes；
  study counts 为 9/13/11/27。
- Pooled metabolic input gate job **599942** 已完成，验证 matrix、manifest、
  Human-GEM 和全部 SHA256 provenance。Receipt：
  `data/processed/corneto/pooled_primary_60/metabolic_input_gate.json`。

### NMF

- Primary cohort 与 pooled NMF jobs **591326-591328** 已完成。
- Pooled rank-3 clusters 包含 23/20/17 samples；cophenetic correlation 为
  0.939，silhouette 为 0.770。Rank 2 更稳定（0.972/0.899），因此 rank 3
  是冻结的 comparability benchmark，不是唯一最优 rank。
- Pooled-vs-cohort rank-3 ARI：7223 为 0.723、10801 为 0.529、11000 为
  0.377、14568 为 0.569；mapped assignment agreement 为 0.889/0.846/0.727/0.741。
- pooled state 与 study 的关联未达到清晰统计证据（chi-square p=0.207、
  Cramer's V=0.265），但仍需要 study-aware sensitivity。
- Patient-balanced-52 rank-3 comparison 已完成：mapped agreement 0.942、
  ARI 0.842、NMI 0.825。这支持对 repeated-patient weighting 的稳健性，但仅是一次 deterministic patient-balanced selection。

### Regulatory CORNETO

- True multi-condition normalized-lambda grid 与 retries 已完成；final summary
  包含九个 nominal lambdas 下的 45/45 pooled/cohort receipts。
- Patient-balanced regulatory analysis 在 52 个共同患者上保留了相似网络：
  pooled/balanced union Jaccard 为 0.890，mean per-sample Jaccard 为 0.899。
- Narrow-vs-richer PKN sensitivity 已完成。Pooled union Jaccard 为 0.202，
  mean sample Jaccard 为 0.108；cohort union Jaccard 约为 0.108-0.143。因此
  network conclusions 明显受 PKN 影响，必须报告 stable cores 与 uncertain alternatives。
- Regulatory longitudinal summary 覆盖 60 runs 和八个 within-family transitions。
  这是 response-blind 分析；acquired-resistance 解释仍需 exact exposure 和 phenotype。
- Regulatory x NMF state integration 已覆盖四个 cohorts。它是基于相同 input 的
  descriptive integration，不是 independent validation，也不是 subtype discovery。
- E-MTAB-14568 regulatory alternative-optima ensemble 与 serial repairs 已完成。
  正确的证据对象是 solution ensemble，而不是单个 sparse optimum。

### Human-GEM 与 TPI1 model gate

- Public Meeson TPI1 table audit job **591424** 已完成；它审计 published/model
  values，不是 OCM-specific knockout result。
- Independent model-level TPI1 gate job **599873** 已完成且有效：Human-GEM
  含 13,096 reactions 和 3,628 genes；TPI1 映射到 `ENSG00000111669` 与 reaction
  `HMR_4391`；没有调用 solver。Receipt：
  `data/processed/corneto/tpi1_model_gate.json`。
- 过时的 pending job **591049** 在启动前已取消，因为它会在 original baselines
  之后覆盖同一 receipt。其职责现已拆分为 599873（model only）和 599950（strict
  receipt-dependent preflight）。

## 终态任务审计与 retry lineage

下表记录最终 evidence，而不是把每个诊断性 Slurm attempt 当成独立实验计数。

| Analysis family | 最终 evidence 状态 | 失败/被取代尝试及处理 |
|---|---|---|
| RNA quantification 与 aggregation | 四个 aggregation receipts 均为 `completed`：36/36/12/33 runs；每个均有 60,609 genes、227,462 transcripts | 早期逐-run download、OOM 与 aggregation failures 均已修复；没有剩余 RNA retry |
| Primary 与 pooled NMF | 591326-591328 已完成；patient-balanced NMF/compare 592020、592094 已完成 | 592083 comparison 失败，已由 592094 替代 |
| True multi-condition regulatory lambda grid | Final retry 591593 与 summary 591595 已完成；pooled + 四 cohorts、九个 lambdas 共 45/45 receipts 有效 | Initial 591416 tasks 遇到 Gurobi session limit；失败证据保留，仅重跑受影响 labels |
| Regulatory alternative optima | Final summary 591572 验证 E-MTAB-14568 全部 27 samples：26 个 nonempty completed、1 个 zero-edge blocked | 六个初始 session-cap errors 已由 591569 串行 retry |
| Richer-PKN sensitivity | 592021/592023/grid jobs 与 comparison 592118 均完成 | 无未解决 retry |
| Patient-balanced regulatory sensitivity | 592143 与 592149 在 52 common patients 上完成 | 无未解决 retry |
| Longitudinal regulatory summary | 592019 完成：60 runs、八个 response-blind within-family comparisons | 早期 prototype 588876 失败，已由 588883/592019 lineage 替代 |
| Regulatory x NMF integration | Final v4 job 592053 完成，覆盖 9/13/11/27；576 个 BH-adjusted edge-state tests 均无 q<0.05 finding | 592018、592040、592046 在 interface/path 修正阶段失败，已由 592053 替代 |
| Meeson public evidence | 591424 完成；toy joint/order/global-retention receipts 验证 algorithmic behavior | 7223/10801/11000 的 cohort-specific order/ensemble jobs 因 dependency 被取消，尚未重新绑定有效 metabolic retries；14568 tasks 仍 pending |
| Metabolic growth-fraction sensitivity | 588286-588289 没有有效 receipt | 四个任务均达到原 4 h limit；primary baselines 验证前不排 retry |

对 normalized regulatory grid 而言，高 lambda 下 solver completed 不等于正向 network
finding：pooled network 在 nominal lambda 0.05 及以上变为空，大部分 cohort network
也在 0.05-0.1 collapse。这是当前 scaling 下 over-regularisation 的证据，不是 biological
absence 的证据。lambda 0.001 时 pooled-vs-merged-cohort edge-union Jaccard 为 0.746；
lambda 0.01 时降至 0.286。这些仍是 response-blind technical results。

## Metabolic baseline：运行中、失败与排队任务

冻结的 primary settings：Human-GEM v1.4.1、raw TPM 经 `log1p` 转换、primary
tumour only、candidate budget 25、growth fraction 0.9、independent lambda 0.1、
joint lambda 1.0，以及不允许 fallback 的 explicit Gurobi。

所有 monolithic cohort attempts 均已终态，且都没有有效 final receipt。
588250/588252 与 600004/600005/600007 达到 wall-time limits；588251 是已确认的
64G OOM；588253 与 600006 失败但没有 canonical scientific receipt。后续的 11000
retry 727583 也达到 72 h。这证明只在最后写一次 JSON 的策略不适合本 MILP。

替代设计把每个 cohort 的 expression-derived bounds 与 objectives 冻结在
`checkpoint_b25/context.json`；一个 array task 只求解一个 independent OCM；全部
independent canonical receipts 成功后才求 joint cohort，最后组装
`full_direct_b25.json`。本次更新状态如下：

| Study | 必需 OCM receipts | 此前 r4/legacy 终态 | 已取消的排队中 64G tasks | 当前 recovery array |
|---|---:|---|---|---:|
| E-MTAB-7223 | 9 | 805860_5：70 h partial incumbent；其余 OOM | 没有剩余 pending task | **834320**，128G，tasks 0-2 运行；3-8 排队 |
| E-MTAB-10801 | 13 | 805861_1：70 h partial incumbent；其它已启动 tasks OOM | 805861 tasks 6-12 | **834321**，128G，tasks 0-2 运行；3-12 排队 |
| E-MTAB-11000 | 11 | 805003_0：72 h Slurm TIMEOUT，无 receipt | 805862 tasks 0-10 | **834323**，128G，`0-10%3`，等待其它三个 r5 arrays |
| E-MTAB-14568 | 27 | r4 已启动 tasks 均 OOM | 805863 task 8 与 tasks 12-26 | **834322**，128G，tasks 1/3/5 运行、6-26 排队；r5 tasks 0/2 为 OOM、4 为 SIGBUS；**863034**，256G，等待 repair |

00:24 EEST live check 时，六个24-hour legacy tasks 已运行13 h 13 min，72-hour
11000 task 已运行11 h 44 min。Solver-step CPU efficiency 约94-96%，disk counters
继续增加，peak RSS 为15.2-23.8 GiB，相对于64G request 没有 memory pressure。这支持
“仍在 active solving”，而不是 OOM 或 idle hang；但不能证明 feasibility、optimality 或
scientific correctness。这些 pre-instrumentation processes 没有 live Gurobi progress log，
也没有写 canonical/partial receipt、`.sol` 或 `.mst`，因此当前无法取得它们的 incumbent
与 gap。六个24-hour tasks 尚余约10 h 47 min，仍可能重复此前 timeout pattern。只有
instrumented jobs 提供 live solver logs，并在 solver 正常返回时写出 internal-limit
telemetry 与 partial receipt；OS 强制终止或 filesystem 写入失败仍可能阻止保存。

749576/749580/749584 的六个健康 RUNNING tasks 与健康的 805003_0 均未取消。
749576/749580/749584 的 tasks 0-1 此前已在 24 h timeout，且没有 independent
receipt。instrumented recovery 在运行前验证并跳过 context hash 一致的 canonical
receipt，因此 legacy task 若成功，其结果不会被重复计算。

Instrumented arrays 启动后，7223、10801 与14568 的多批 r4 elements 在原64G request
下通常运行1-9分钟即被 Slurm memory cgroup 终止。Python steps 被强制杀死前没有机会写
canonical 或 partial scientific receipt。12:09 EEST 时仍有7个健康的64G tasks 在运行；
这些任务被明确保留。所有仍在 pending 的64G tasks 已取消，包括尚未启动的整个805862
array。independent-job 默认内存已提高到128G，throttle 为3的 recovery arrays
**834320-834323** 已被 Slurm 接受。它们会验证并跳过 context/provenance 匹配的 canonical
receipts，因此只重算 missing/OOM indices。

r5 dependencies 使用显式 Slurm array wildcard（`jobid_*`），而不是只依赖 array master
ID，因此只要仍有一个保留的 r4 task 在运行，recovery 就不会启动。11000 array 还额外等待
其它三个 recovery arrays 与 legacy job 805003，从而把计划中的 solver load 控制在10个
Gurobi licence sessions 以内。这只是 operational recovery；目前没有新增、已验证的 biological
result。

03:42 EEST 检查时，r4 tasks 805863_10 与805863_11 也在运行约8 h 6 min 后以
`OUT_OF_MEMORY` 结束；现有834322 recovery 已覆盖这两个 indices，因此没有提交额外
retry。五个 instrumented r4 tasks 仍健康且持续求解：805860_5、805861_0、805861_1、
805861_5 与805863_0；其 live relative gap 为3.30-3.91%。legacy 805003_0 运行至
43 h 30 min，但仍没有 live telemetry。四个 cohort 均没有 canonical independent receipt
或 joint receipt，因此这些 incumbents 仍只能作为 optimization diagnostics。全部 r5 arrays
均正确处于 dependency 阻塞状态。

17:11 EEST 检查时，805861_0 与805863_0 随后分别在约27 h 与30 h 后以
`OUT_OF_MEMORY` 结束。现有128G recovery arrays 已覆盖这两个 indices，因此没有增加
重复 retry。最后一个14568 r4 task 终止后，834322_0 正常启动，并在约12秒时得到首个
feasible incumbent，relative gap 为17.7%；其余14568 tasks 因 scheduler priority 排队，
并非 dependency 或 licence failure。仍在运行的 r4 tasks 805860_5、805861_1 与805861_5
的 live gaps 为3.26-3.87%。legacy 805003_0 已运行52 h 32 min。四个 cohort 的
canonical independent 与 joint receipt 计数仍为零，因此没有释放任何 biological
interpretation。

随后834322_0 尽管请求128G，仍在5 min 42 s 后发生第二次 cgroup
`OUT_OF_MEMORY`。被杀前 incumbent 为-134.95360、bound 为-143.05700、relative gap
为6.00%，但没有写 canonical 或 partial receipt；这些数值只能用于诊断。因此提交了覆盖
全部27个 indices 的256G repair array **863034**，throttle 为3，frozen context 与 solver
参数均不变。它等待834322全部 tasks 终态，验证并跳过匹配的 canonical receipts，只重算
missing/noncanonical indices。14568 joint job 834326 现依赖
`afterok:863034_*`；这移除了因 r5 task 失败而永远无法满足的旧 dependency，同时保持
fail-closed progression。

09:50 EEST 检查时，805861_5 也在31 h 5 min 后以 `OUT_OF_MEMORY` 结束，834322_2
则在4 h 后以 `OUT_OF_MEMORY` 结束。现有 recovery arrays 834321 与863034 已分别覆盖
相应 indices，因此不需要新增 successor。active solver tasks 805860_5、805861_1、
834322_1、834322_3 与834322_4 保持健康，live gaps 为3.39-3.83%；legacy 805003_0
已运行69 h 11 min，但没有 live telemetry。其余 r5 14568 tasks 受 array throttle 限制，
r6 正确等待 r5 终态。四个 cohort 的 canonical independent 与 joint receipt 计数仍为零。

### 2026-08-27：完成审计的 full-duration partial outputs

两个 instrumented r4 tasks 达到252000秒 solver limit 后正常返回，并在 Slurm 强制终止前
保存了结果。两份 receipts 均明确记录 `partial_incumbent`、`scientific_success=false`、
Gurobi `TIME_LIMIT`，且 frozen-context SHA256 匹配。Slurm `FAILED 2:0` 是 partial
output 按设计触发的 fail-closed exit，不是新的 OOM。

| Array task / RNA run | Incumbent objective | Best bound | Relative gap | cohort 的 `checkpoint_b25/instrumented_attempts/` 下的 receipt 文件 |
|---|---:|---:|---:|---|
| 805860_5 / ERR2808261 | -136.0536082251 | -141.1906456601 | 3.775745% | `005_ERR2808261_job833065_task5.json` |
| 805861_1 / ERR6389069 | -136.2536081782 | -141.1262291334 | 3.576141% | `001_ERR6389069_job832834_task1.json` |

两份 receipts 引用的 `.sol`、`.mst` 与 `.gurobi.log` 均存在且非空，`summary_error` 均为
null。这验证了 result persistence，而非 optimality 或 biological conclusions。gap 仍高于
未改变的 `MIPGap=0.0001`（0.01%）标准。已部署 runner 与所审代码的 hash 一致：
runner 会写 `.mst`，但 retry 不会重新加载它。因此目前 recovery 复用的是 canonical
completed samples，而非 partial incumbent 或保存的 branch-and-bound tree；不能声称
已实现 warm-start reuse。

现有834320/834321 arrays 已自动启动，并覆盖这两个 partial tasks，本轮未追加 retry。
legacy 805003_0 在72 h timeout 后仍无 receipt，由834323 覆盖。当前9个128G solver
tasks 正在运行：834320、834321、834322 各3个，live gaps 为3.19-4.65%。使用实际
step JobId 查询的 `sstat` RSS 约5.2-29.6 GiB，CPU/I/O counters 非零；这些 container
accounting 采样不能保证后续不会 OOM。四个 cohort 仍没有 canonical independent 或
joint receipts，audit、comparison 与 TPI1/FVA gates 均保持关闭。会话内 monitor 已明确
转向 r5/r6 chain，保留原 schedule，且没有 model override。

### 2026-08-27 partial-solution scientific audit（2026-08-28 登记）

本次 read-only audit 使用上表两份70-hour attempt receipts、对应 `.sol` 和 `.mst`、
冻结的 `context.json`、已部署的 objective/indicator code，以及 Human-GEM v1.4.1。
Model SHA256 仍为
`57d1b137f0c90d83a3e4f9a8225d74d37523594e6ee99f622b160a014d9f7050`。
Context SHA256 分别为
`7ea9d2268ec7647bfa0f47f8215913442cafc2b6f76faefa13a40c402b7fcb1b`
（7223）和
`1475eee9397af6644fcfaa6500594fbe1441b7d4d4acfbf714a66bac91f0792c`
（10801）。以下内容是 model diagnostics，不是已验证的 biological findings。

| Saved-solution quantity | OCM66-1 / ERR2808261 | OCM110-9 / ERR6389069 |
|---|---:|---:|
| 实际施加的 expression-derived reaction bounds，不含 biomass | 15 | 1 |
| Selected indicators，阈值 >=0.5 | 513 | 511 |
| Nonzero fluxes，绝对值 >1e-7 | 544 | 542 |
| Biomass flux | 187.3536299766 | 187.3536299766 |
| 有 nonzero flux 但 indicator <0.5 的反应数 | 31 | 32 |

- Cohort candidate budget 为25，不代表每个 OCM 都有25条有效 expression bounds。
  `scripts/run_corneto_14568_pilot.py` 的 `_candidate_sets` 仅保留
  `proposed_upper > 0` 的 candidates，`_reaction_bounds` 会跳过 missing candidates。
  因此 zero-expression candidates 不会自动关闭反应。修订 biological analysis 前需要
  区分 missing 与 zero expression；两者均不能成为未经审查就设置零边界的理由。
- 所有 expression-capped reactions 在保存的 flux reporting threshold 下均为零。
  两份 flux vectors 交换到对方的 expression 与 biomass bounds 下，在1e-7 tolerance
  下也均无违规。SBML mass-balance 最大 absolute residual 小于3.34e-9，且没有超过
  1e-6 的 model-bound violation。这验证了 reciprocal flux feasibility，不能证明完整
  feasible sets 相同，也不能证明网络只能对应某一位患者。
- 两份解均通过 `EX_atp[e]`、`EX_pep[e]`、`EX_pcreat[e]` 摄取 ATP、
  phosphoenolpyruvate 与 phosphocreatine，三个 flux 均为-1000。未修改的模型允许这些
  exchanges；当前 context 未根据 measured uptake 或实际 culture medium 校准。
  因此 biomass 与 energy-pathway 结果不能解释为已测量的 OCM physiology。
- 31/32条 flux-indicator discrepancies 的 flux magnitude 约0.00386-0.00906，
  `.mst` 中 binary values 约3.86e-6-9.06e-6。已部署约束为
  `lb*y <= v <= ub*y`，边界可达1000；这些现象符合 integer tolerance 被放大造成的
  trickle flow。Flux mass balance 通过，不等于 indicators 取整后的网络通过验证。
  活跃反应分类前需要有依据的 tightened bounds 与 fixed-indicator feasibility audit。
  参见 [Gurobi IntegralityFocus 说明](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#integralityfocus)。
- 已部署的 single-sample objective 为 `-biomass + 0.1*sum(indicators)`。
  用 `.sol` 重算与 incumbent objective 一致；约0.2的差异来自 sparsity，而非 biomass
  不同。两份 indicator sets 的 intersection 为402、union 为622（Jaccard 0.6463），
  但这不是 HGSOC-specific core。解使用了 PPP reactions `G6PDH2c`/`PGLc` 和
  TPI1 reaction `HMR_4391`；使用某反应不等于 essentiality、enrichment 或 drug response。

Input/media revision、expression-null controls、integrality checks、alternative-solution/FVA
analysis 与 patient-balanced contrasts，仍是这些解释尚未完成的 follow-up gates。
只降低 MIP gap 不能替代这些检查。本次 audit 没有修改 media、expression policy、lambda、
integrality tolerance 或任何 running job。

### 2026-08-28 11:15-11:17 BST：live recovery 与已覆盖的 SIGBUS failure

以下 snapshot 来自 targeted `squeue`/`sacct`、当前 solver-log tails，以及使用实际
numeric solver-step IDs 的 `sstat`。Gap 为 live log 中的舍入值，不是 final receipt，
也不是完成时间预测。

| Array task | Actual numeric JobId | 检查时 elapsed | Live gap |
|---|---:|---|---:|
| 834320_0 | 890339 | 27 h 54 min | 3.06% |
| 834320_1 | 890340 | 27 h 54 min | 3.72% |
| 834320_2 | 890341 | 27 h 54 min | 3.37% |
| 834321_0 | 890279 | 28 h 3 min | 3.47% |
| 834321_1 | 890280 | 28 h 3 min | 3.65% |
| 834321_2 | 890281 | 28 h 3 min | 4.19% |
| 834322_1 | 862558 | 68 h 1 min | 3.15% |
| 834322_3 | 864522 | 66 h 35 min | 3.42% |
| 834322_5 | 900491 | 21 h 19 min | 3.27% |

九份 logs 均有近期更新。Actual-step RSS 采样约18.0-45.6 GiB，CPU/I/O counters 非零；
与此前一样，container accounting 不能证明128G allocation 不会 OOM。14568 两个
最久的 solves 正接近 internal 70 h limit，不代表即将保证收敛。

本次新增审计失败 **834322_4**，actual JobId **866509**，对应
**ERR13907041 / OCM288-4**，于2026-08-27 15:57:28 EEST 结束，运行
42 h 38 min 29 s。Accounting 为 `FAILED 7:0`；最终输出
`logs/met-inst-ind-14568-r5-834322_4.out` 明确记录 Singularity wrapper 的
**Bus error** 与 srun exit 135。这确认了 SIGBUS，不是已确认的 OOM 或
Gurobi session-cap failure；底层 runtime/node/filesystem 原因仍未确定。
Solver log 最后记录 incumbent -136.25361、bound -141.20775、gap 3.64%。
该 attempt 没有生成 canonical/attempt JSON、`.sol` 或 `.mst`，仅保留 progress-log
证据。现有 write-on-return instrumentation 仍不能保护进程突然退出时的解文件。

已提交的 **863034** 请求256G、throttle 为3，既覆盖此前 OOM indices，也覆盖本次
missing/noncanonical index。其 dependency 仍为 `afterany:834322_*`，joint **834326**
仍等待 `afterok:863034_*`。增加内存不是已验证的 SIGBUS 修复。本次没有追加或重复 retry，
没有取消任何健康任务。其它 joint、assembly、comparison 与 TPI1/FVA gates 的依赖保持
不变。11000 array **834323** 仍用 explicit wildcards 等待三个 r5 independent arrays。

Canonical independent receipts 仍为 **0/9、0/13、0/11、0/27**；四个 joint 与
`full_direct_b25.json` 均不存在。此前两份70-hour partial receipts 的 context 仍匹配，
artifacts 仍保留；此 snapshot 没有新的 r5 partial receipt。没有释放 biological
interpretation，也没有提交新的 patient-level 或 pooled research job。

### 2026-08-29 12:05-12:08 EEST：node-local failures 与已覆盖的10801 repair

Targeted audit 显示九个 active solvers，canonical independent receipts 仍为
**0/9、0/13、0/11、0/27**。7223 的 active tasks 834320_0-2 已运行约
50 h 43 min，live gaps 为3.00%、3.59%、3.28%；10801 的834321_7/10/11
已运行约4-6 h，gaps 为3.86%、4.00%、4.24%；14568 的834322_5/6/7
已运行约19-44 h，gaps 为3.21%、3.82%、4.38%。九份 logs 均持续更新，
actual solver-step RSS 约5.7-40.1 GiB，CPU 与I/O counters 非零。

10801 的三个长任务834321_0-2（actual JobIds 890279-890281）均在 node
`rc5140` 运行约44-47 h 后发生 Singularity-wrapper SIGBUS；最后 live gaps
分别为3.43%、3.61%、4.14%，均未写 attempt 或 canonical receipt。随后六个
tasks（834321_3-6、834321_8-9）也落在 `rc5140`，并在7-10秒内失败：
`srun` 无法执行 `/scratch/project_2012997/xiaomei/hgsoc_corneto_env/bin/python`，
报告 `No such file or directory`。这一 common-node pattern 是 node/filesystem/runtime
failure 的 operational evidence，不是 OOM、solver infeasibility 或 Gurobi
session-cap error。

因此提交了一个 fail-closed 10801 repair array **937737**（`0-12%3`、128G），
frozen context 与 solver parameters 不变。它等待 `afterany:834321_*`，排除
`rc5140`，验证并跳过以后可能生成的 canonical receipts，只重算 noncanonical
indices。Joint job 834325 已改为依赖 `afterok:937737_*`。文档所列 recovery chain
中仍 pending 的 elements 与 downstream solver jobs 也排除了 `rc5140`；没有修改或
取消已经运行的健康任务。

14568 的834322_1 与834322_3 都达到252000-second solver limit，随后发生 wrapper
SIGBUS。834322_1 没有写出 attempt receipt；834322_3 在 SIGBUS 前完成原子写入，
经审计为 `partial_incumbent`、`scientific_success=false`、context SHA256 匹配、
Gurobi `TIME_LIMIT`，objective -136.6536067185、best bound -141.3243116334、
relative gap 3.417916%；其 `.sol`、`.mst` 与 Gurobi log 均存在且非空。这只能作为
optimization evidence。现有256G repair array 863034 已覆盖两个 noncanonical
indices，因此没有重复提交14568 retry。全部 joint、assembly、comparison 与
TPI1/FVA scientific gates 继续关闭。

### 2026-08-30 12:00-12:03 EEST：7223 partial receipts 与 repair serialization

7223 r5 的前三个 tasks 都达到未改变的252000-second Gurobi limit。834320_0 与
834320_1 在 Singularity wrapper 发生 SIGBUS 前原子写出了经审计的
`partial_incumbent` receipts。两份 receipts 均为 `scientific_success=false`，
7223 context SHA256 匹配，Gurobi `TIME_LIMIT`，requested `MIPGap=0.0001`，
且 `.sol`、`.mst` 和 solver log 均存在且非空。两者的
objective/bound/gap 分别为 -137.0536041377/-141.1212280161/2.967907% 和
-136.2536156425/-141.1159266378/3.568574%。834320_2 也达到 solver limit，
最后 log 为 objective -136.85361、bound -141.29013、gap 3.24%，但 SIGBUS
发生在 attempt receipt 写入前。三者都不是 canonical scientific result。

由于834320 的 noncanonical indices 没有既有 successor，因此提交了一个 fail-closed
7223 repair array **948765**（`0-8%3`、128G），context 与 scientific/solver
parameters 不变。它等待 `afterany:834320_*`，排除已知故障 node `rc5140`，
验证并跳过 canonical receipts，只重算 noncanonical indices。Joint job 834324
现依赖 `afterok:948765_*`。

11000 independent array 834323 已改为串联等待三个 repair arrays
（`948765_*`、`937737_*`、`863034_*`）。因此最多只有三个 throttle-3 repair
arrays 同时请求 Gurobi sessions，不会再叠加11000 array；这保持 operational ceiling，
没有修改 scientific parameters。另外，10801 task 834321_10 在运行16 h 55 min 后
被其128G Slurm memory cgroup OOM kill；现有937737 repair 已覆盖它，因此没有新增 retry。

本次检查有九个 active solver tasks。7223 tasks 3-5 的 live gaps 为
4.25/4.18/4.24%，10801 tasks 7/11/12 为3.59/3.59/3.72%，14568 tasks 5-7
为3.18/3.58/4.28%。所有 logs 均持续更新；actual-step RSS 约6.0-53.9 GiB，
CPU/I/O 非零。Canonical counts 仍为 **0/9、0/13、0/11、0/27**；没有释放任何
joint、assembly、comparison 或 TPI1/FVA scientific gate。

r5 instrumented independent tasks 请求128G、8 CPU、72 h Slurm limit，同时向 Gurobi
显式传入 `TimeLimit=252000` 秒（70 h）、`MIPGap=1e-4`、8 threads、seed 0，留出
原子写 receipt 的时间。存在 incumbent 时，receipt 记录 objective、best bound、
absolute/relative gap、status、solution count、nodes/work/iterations 与模型维度，并保存
`.sol`、`.mst` 和 solver log。达到 time limit 且有 incumbent 时状态为
`partial_incumbent`，进程非零退出阻断 `afterok`；它明确不是 biological result，也不是
canonical cohort receipt。

Smoke job **805824** 已在 E-MTAB-7223 OCM ERR2808250 验证该机制。600.0 秒时有
10 个 feasible solutions，incumbent objective -134.7536165、best bound
-142.6087313、relative gap 5.8292%，探索 14,986 nodes。任务原子写入
`partial_incumbent` receipt、`.sol` 和 `.mst` 后按设计 exit 2。这证明 observability
与 fail-closed control 有效，但不等于 scientific completion。本次更新时四队列的
canonical independent receipts 仍为 0/9、0/13、0/11、0/27。

替代 instrumented joint jobs 为 **834324-834327**，assembly jobs 为
**834328-834331**。7223/10801/11000 的 joint memory 为196G，14568 为384G，均用
70 h internal solver limit。834326 由256G repair array 863034 显式 gate。旧 audit
**805872** 与 strict comparison **805873** 已被替代
并取消。新 audit **834332** 与 strict comparison **834333** 通过 `afterany` 等待四个 r5
assembly，并对 missing/partial receipts fail closed。TPI1 preflight **834334** 与 full
WT-versus-delta-TPI1 FBA/FVA **834335** 必须通过 strict comparison 的 `afterok`。只有在
Slurm 接受完整 r5 replacement chain 后，旧的 r4 downstream jobs 805864-805875 才被取消。

## Pooled metabolic joint analysis

该分析有意只做 joint-only：四个 cohort jobs 提供 independent/cohort evidence；
不在一个 72 小时 pooled job 中重复 60 个 independent MILPs，以免重复计算并在超时后丢失进度。

Pooled input gate **599942** 仍然有效，但三个 final-only pooled smokes 均发生 operational
failure：600008 在12 h timeout、615508 在36 h timeout、727584 在72 h timeout，均没有
scientific receipt；对应 full-60 与 comparison successors 从未运行即被取消。因此停止重复
相同的 opaque pooled job。

下一次 pooled attempt 必须先构建冻结的 pooled checkpoint context，并使用同一个
instrumented joint runner：先用2-4 conditions 的 bounded smoke 暴露
incumbent/bound/gap 并写 partial receipt，之后60-condition job 才有资格启动。本次没有
提交新的 pooled solver，因为 cohort checkpoint contexts 不能静默替代 pooled context。
在四个 canonical cohort receipts 与 instrumented pooled receipt 都验证前，
pooled-vs-cohort biological interpretation 保持 blocked。

## TPI1/FVA 与 Meeson dependency chain

```mermaid
flowchart LR
    M["599873 model-only TPI1 gate: valid"]
    B["834333 four strict metabolic receipts"] --> P["834334 receipt-dependent TPI1/FVA preflight"]
    P --> V["834335 WT vs delta-TPI1 FBA/FVA"]
```

现有 Meeson order-sensitivity 与 ensemble jobs 绑定在各个 original cohort jobs 上。
如果 original job timeout 但 conditional retry 成功，对应的 cancelled downstream
task 必须改为依赖有效 retry receipt 后重新提交。不能仅依据 dependency state 推断结果。

## 仍需实现/运行的 reference comparisons

1. **Tumour vs stroma（正式 secondary analysis）。** 分别分析 E-MTAB-7223 与
   E-MTAB-10801；身份允许时使用 paired contrasts，然后 meta-analyse prevalence/effect
   sizes。报告 tumour-enriched、stroma-enriched、shared stable core、study-specific
   与 unresolved features。
2. **HGSOC vs confirmed non-HGSOC（exploratory）。** 先确认 histology；ambiguous
   samples 仍保持独立类别。
3. **Tumour vs 两个 cell-line controls（descriptive QC）。** 不应作 population-level
   p-values 或 population-level claims。
4. **Passage/non-Tighe stability。** 使用 paired/case-level expression、NMF、regulatory
   与 metabolic retention；不要附加 drug-response labels。
5. 所有 network contrast 都应以 prevalence difference、patient/study-aware uncertainty、
   alternative-optima frequency、lambda/PKN sensitivity 与 multiplicity control 代替简单 set subtraction。

## 缺少的数据与 blocked claims

目前尚未获得通过验证的 intake table，其中包含 exact OCM-linked AUC、raw GI50 units、
cumulative treatment exposure 及其 primary provenance。在通过 frozen phenotype gate 之前，
项目不能声称 Taxol response、intrinsic resistance、acquired resistance prediction 或
response-model accuracy。`chemo_naive_at_biopsy` 不能替代 exact cumulative exposure。

## Operational update rules

1. 保留失败 logs 与 receipts；retry 只有在 canonical result 缺失时写入，已有 receipt 必须先验证后才能跳过。
2. Audit jobs 使用 `afterany`，确保失败可见；scientific downstream jobs 使用 `afterok` 加内部 receipt validation。
3. 每个 receipt 记录 job ID、exact command、model/input hashes、sample contract、solver、
   lambda、candidate budget、objective 与 claim limit。
4. 不得根据最有吸引力的 biological story 选择 lambda、NMF rank、PKN 或 candidate budget；
   primary settings 与 sensitivity grids 必须在 outcome analysis 前冻结。
5. 任何任务进入终态或 validated receipt 改变 evidence state 时，都要更新本登记。
6. 每个 long MILP 必须设置短于 Slurm 的 internal solver time limit，持久化
   incumbent/bound/gap 与 solver artifacts，并区分 `partial_incumbent` 和 canonical
   `completed`。Partial receipt 可指导 recovery，但不能释放 scientific downstream jobs。

## 定时监控（Recurring monitor）

Codex heartbeat **HGSOC CORNETO Roihu pipeline monitor** 已绑定到当前对话，
每 30 分钟运行一次；它没有显式 model 或 reasoning override，因此遵循当前对话/默认设置，
不会创建单独的 standalone monitoring conversation。这个时间间隔的理由是：

- Slurm dependencies 会自动启动有效 successors，因此分钟级轮询不会加速 pipeline。
- 30 分钟足以在 instrumented checkpoint chain 中发现 smoke/startup、license-session、
  timeout 与 receipt failures。
- 对 branch-and-bound 中日志很少的 multi-hour Gurobi jobs，30 分钟不会造成无效的高频查询。

每次运行都检查 scheduler state、resource use、log tails 与 receipt JSON。它只能执行确定且
属于当前范围的修复：不重复提交、不取消健康任务、不改变 scientific parameters，也不在没有
有效 receipt 时作结果声明。OOM/TIMEOUT retry 保留参数并只增加有依据的 resources；Gurobi
session-cap failure 要等 active sessions 少于 8 后再低并发或串行 retry。有实质变化时更新本文件
并 push 到 delivery branch；无变化时只报告简短 checkpoint。所有 cohort、pooled、comparison
与 TPI1/FVA outputs 进入终态并完成 audit 后，monitor 应报告完成并暂停。
