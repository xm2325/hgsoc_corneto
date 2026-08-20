# HGSOC CORNETO 研究状态与依赖登记（中文对应版）

最后运行更新：2026-08-20 14:20 BST（16:20 EEST）。本文件与
`docs/hgsoc_corneto_research_status.md` 对应，记录研究范围、已完成证据、排队分析、失败尝试、依赖关系与可声明范围。
仅有 Slurm `COMPLETED` 不足以证明科学分析完成；只有输出 `receipt` 通过相应内容验证后，结果才算科学上完成。

## 运行纠正：2026-08-20

首轮 checkpoint preparation jobs 727669、727673、727677 在写出 context
receipt 前失败。原因是三个 solver 入口没有显式导出项目 WLS license file，
Gurobi 因而选择 size-limited local license，并拒绝 8,400-row/26,192-column
model。这是 execution-environment defect，不是 OOM、model infeasibility 或
biological result。prepare、independent、joint 三个 sbatch 现均导出并检查与
既有 solver jobs 相同的 `GRB_LICENSE_FILE`。

修正后的 chains 为 749575--749578（E-MTAB-7223）、749579--749582
（E-MTAB-10801）及 749583--749586（E-MTAB-14568）。startup gate 时，
749575、749579 与 749583 均已 exit 0，并分别为 9、13 与 27 个 samples
写出 `prepared` context receipts。日志确认 WLS academic license 2849103。其 independent arrays 同时依赖 context preparation 成功及 active jobs
727583/727584 终态，以保留 solver-session safety margin。当前仍无 scientific
receipt，因此这些 jobs 尚不能支持 biological claim。

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

## Metabolic baseline：终态失败与 checkpointed recovery

冻结的 scientific settings 保持不变：Human-GEM v1.4.1、raw TPM 经 log1p 转换、
primary tumour only、candidate budget 25、growth fraction 0.9、independent lambda
0.1、joint lambda 1.0，以及 explicit Gurobi 且不允许 fallback。

旧 monolithic runner 依次求解全部 independent conditions，再求 joint problem，并且只在
进程结束时写唯一 receipt。终态证据如下：

| Study | Original / retry evidence | 2026-08-19 final receipt |
|---|---|---|
| E-MTAB-7223 | 588250 在 24 h TIMEOUT；600005 在 72 h TIMEOUT（step MaxRSS 52.6 GB） | missing |
| E-MTAB-10801 | 588251 在 64G OOM；600004 在 72 h/196G TIMEOUT（step MaxRSS 60.5 GB） | missing |
| E-MTAB-11000 | 588252 在 24 h TIMEOUT；600006 在 43 h 30 min/128G OOM | missing；targeted 196G retry **727583** 正在运行 |
| E-MTAB-14568 | 588253 在 65 h 32 min/128G OOM；600007 在 72 h/196G TIMEOUT（step MaxRSS 57.4 GB） | missing |

每个 retry log 在终止前都包含精确的 expected independent LP reads（9/13/11/27）。
这支持 execution-stage diagnosis：independent work 已完成，但在进入或尝试 joint stage
时没有被保存；它不构成 biological result。Availability audit **599875** 正确报告
status=incomplete、0/4 valid receipts 与 final_comparison_permitted=false；strict
comparison **599876** 因缺 E-MTAB-7223 receipt 而 fail closed。TPI1 jobs
599950/599951 随后由 dependency 自动取消，没有产生 scientific output。

现已实现不改变 scientific parameters 的 checkpointed recovery。Frozen context 记录
samples、candidate IDs、bounds、model/input SHA256 与 solver；每个 independent OCM
原子写 receipt；joint solve 单独获得完整 72 h；只有全部 condition receipts 与 joint
receipt 均为 optimal 且 context hash 一致，才组装 canonical full_direct_b25.json。
Checkpoint fail-closed assembly tests 与既有 joint-FBA tests 合计 4/4 通过。

| Study | Prepare | Independent array | Joint | Assemble |
|---|---:|---:|---:|---:|
| E-MTAB-7223 | 727669 | 727670（0-8，%2） | 727671（196G/72h） | 727672 |
| E-MTAB-10801 | 727673 | 727674（0-12，%2） | 727675（196G/72h） | 727676 |
| E-MTAB-14568 | 727677 | 727678（0-26，%2） | 727679（384G/72h） | 727680 |

Independent arrays 等待 727583 与 pooled smoke 727584 终态后释放；三组 array 合计
最多 6 个 Gurobi tasks，低于项目规定的 8 个 active solver workflows 安全上限。
Audit **727685** 与 strict comparison **727686** 等待三个 assemble jobs 与 727583，
并保持 fail-closed。

## Pooled metabolic joint analysis

Pooled analysis 保持 joint-only，不重复 60 个 independent MILPs。之前两个
four-condition one-per-study smoke 都没有 receipt：**600008** 在 12 h TIMEOUT，
**615508** 在 36 h TIMEOUT。后者 MaxRSS 仅 11.5 GB，且无 OOM、license-session 或
Python exception，因此只再提交一次 scientific settings 完全相同的 **727584**，
使用最大 72 h wall time 与 128G，并以 canonical smoke receipt 不存在为 duplicate guard。

只有 727584 成功，**727684** 才会启动 pooled-60 joint solve（384G/72h）。
Pooled-vs-cohort comparison **727689** 同时依赖 727684 与 strict four-cohort
comparison 727686。早期 600009/600010 和 615515/615517 因 afterok 失败被取消，
从未运行；它们是 provenance，不是 negative scientific result。

## TPI1/FVA 与 Meeson dependency chain

独立 model-only TPI1 gate 仍有效：Human-GEM 含 13,096 reactions 与 3,628 genes；
TPI1 映射到 ENSG00000111669 和 reaction HMR_4391。它不依赖 cohort receipts，也不是
OCM knockout result。

新的 scientific chain 为 strict comparison **727686** -> receipt-dependent preflight
**727687** -> WT versus delta-TPI1 FBA/FVA **727688**。任何 invalid cohort receipt
都会取消该链。Cohort-specific Meeson order/ensemble analyses 同样必须等待相应
canonical cohort receipt，不能仅凭 scheduler state 释放。

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

## External response-blind validation

External 结果记录在
/scratch/project_2012997/xiaomei/hgsoc_corneto_external/regulatory_external_comparison_v2.json
（status=completed）。Taylor signature 在 external scoring 前按 nominal lambda 0.001
冻结：edge 必须至少出现在 6/52 Taylor patients 与 2/4 cohorts。共有 5 条 signed
edges 通过，signature SHA256 为
e10fb081217ad15ece5119ed003dbc16f76cdf3699a1b78ee7bce41a6cf558e6。

True multi-condition fits 的全部 conditions 均 optimal：GSE277107 ovary/omentum bulk
22/22、GSE189955 patient pseudobulk 59/59、GSE208216 organoid models 14/14。在预先
声明的 50% patient/model prevalence threshold 下，5 条 Taylor edges 中没有任何一条
在任一 external group 成为 consensus feature：ovary n=11、omentum n=11、HGSOC
epithelial-candidate n=12、fibroblast proxy n=12、normal-FT n=6、PDO n=11、FT
organoid n=3。这是对 strict external-consensus claim 的 negative falsification
evidence，并不证明 pathway 不存在。Point prevalences 与 bootstrap intervals 仍只能
descriptive；不能作 drug-response 或 causal claim。

GSE180661 source acquisition 已完成：32,555,276,423-byte HDF5 的 locally observed
SHA256 为 edc5f5d7449478d7dfec1f575c89670bcd1bb041124ef2d53a4df0ecf7e29be6。
因为 GEO 未提供 upstream digest，receipt 正确保留 matrix_identity_frozen=false；
在 review 并冻结 observed hash 前，pseudobulk 仍 blocked。DepMap 同样保持 deliberate
blocked，直到提供明确 quarterly release 及同版 Model、Chronos 与 README files。

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

## 定时监控（Recurring monitor）

Codex heartbeat **HGSOC CORNETO Roihu pipeline monitor** 已绑定到当前对话，
按约两小时的低频 cadence 运行；它没有显式 model 或 reasoning override，因此遵循当前对话/默认设置，
不会创建单独的 standalone monitoring conversation。这个时间间隔的理由是：

- Slurm dependencies 会自动启动有效 successors，因此分钟级轮询不会加速 pipeline。
- 30 分钟足以在 588250/588252 接近 24-hour limits 时发现 smoke/startup、license-session、
  timeout 与 receipt failures。
- 对 branch-and-bound 中日志很少的 multi-hour Gurobi jobs，30 分钟不会造成无效的高频查询。

每次运行都检查 scheduler state、resource use、log tails 与 receipt JSON。它只能执行确定且
属于当前范围的修复：不重复提交、不取消健康任务、不改变 scientific parameters，也不在没有
有效 receipt 时作结果声明。OOM/TIMEOUT retry 保留参数并只增加有依据的 resources；Gurobi
session-cap failure 要等 active sessions 少于 8 后再低并发或串行 retry。有实质变化时更新本文件
并 push 到 delivery branch；无变化时只报告简短 checkpoint。所有 cohort、pooled、comparison
与 TPI1/FVA outputs 进入终态并完成 audit 后，monitor 应报告完成并暂停。
