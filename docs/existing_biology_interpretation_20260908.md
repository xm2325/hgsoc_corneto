# 已有结果：生物学解释、竞争解释与证伪（2026-09-08）

## 结论

中心问题仍合理，但目前只能保留“存在值得验证的表达状态与参考网络候选”的探索性解释，
不能宣布已发现跨队列稳健的 HGSOC 特异代谢—调控状态。下述判断基于
[2026-09-06 科学审计](hgsoc_corneto_scientific_audit_20260906.md)及其链接的
controlled LP audit、RNA registry 和历史 receipt 快照；不是本轮新作业的结果。

## 哪些解释可以保留，怎样推翻它们

| 证据层 | 可保留的解释 | 竞争解释 | 下一步反证检验 |
|---|---|---|---|
| NMF | 表达矩阵中存在可压缩的共同变化模式，值得注释其 pathway composition | cohort、培养条件、患者重复或技术因素产生分组 | 按 patient 隔离留出 cohort，训练集选基因与拟合；比较 held-out reconstruction 与训练均值 baseline；检查 patient balancing |
| 同患者 OCM | 同源样本可能共享患者背景，也可能有空间、时间或培养差异 | 同批实验导致相似；共享编号不等于 longitudinal | 验证 metadata，再比较 within-patient 与 study-pair-matched between-patient distances；按 patient 而非样本对计独立单位 |
| Regulatory | 在指定 PKN、输入输出与惩罚下，存在可解释 RNA 的候选连接 | 被指定为 output 的 TF 必然富集；网络先验驱动；同一 RNA 重复使用 | 同时报告 imposed outputs；逐因素改变 PKN/输出选择；训练折内重新预处理，避免泄漏 |
| Tumour–stroma | 若配对差异稳定，可支持 tumour 相对该 reference 的表达富集 | 谱系与培养条件差异，不是 HGSOC 特异性 | patient-level paired effect、方向一致性、区间及多重检验；随后用独立 tumour/reference 数据验证 |
| Metabolic b25 | 保存的 incumbent 可用于诊断模型与数值问题 | 默认 uptake、弱表达约束与 indicator/flux 不一致产生伪差异 | 首先通过 fixed-indicator feasibility、medium 与 expression-information 对照；之后才解释通路或 TPI1 |

受控 LP 中，20 份 attempts（17 个不同 OCM/run）的 indicator-only 网络均无法达到原 growth floor，
而允许 indicator 与实际 flux-support 并集后均可行。该证据否定这些保存 selection 已是有效
growth-support networks 的说法，不是否定细胞生长或 CORNETO 方法。
数字来自 `evidence/fixed_indicator_lp_controlled_audit_20260906.json` 的逐 attempt 检验，
不是独立患者样本量。Gap 衡量优化界限，不是生物学可信度。

## 本轮新增分析能回答什么

- Patient-grouped NMF 检验的是表达重构的 transferability；通过也不等于疾病亚型已独立复现。
- Hallmark over-representation 为 NMF genes 提供探索性 annotation，不是独立 DE 或 GSEA；
  同一 RNA 上的 NMF–regulatory 一致只能叫 internal consistency。
- Paired tumour–stroma pathway scores 检验相对参考富集，不测量 metabolic flux。
- Real/shuffled expression caps 和 energy-uptake exclusion 检验模型是否感知表达与边界条件；
  未校准的 cap 不等于 enzyme capacity，排除三种 uptake 不等于完整 OCMI medium。

## 14:31 UTC 后新增 receipt 检查

正式 NMF 1121179 的8个 folds 均无 patient overlap，8个输入 hash 全匹配，8个 NPZ 存在。
Rank2/rank3 在四个 held-out cohorts 的 patient-mean reconstruction 改善分别为
12.06%/13.85%、10.62%/14.79%、9.34%/14.26%、8.92%/12.62%
（顺序7223、10801、11000、14568）。支持存在可迁移表达结构；更高 rank 的模型容量
本身也改善拟合，不能仅据改善幅度选定“正确亚型数”。

Biology receipt 的10个顶层输入 hash 已匹配。7个重复患者中6个表现出更低的
within-patient expression distance，OCM327反向；应优先核对该患者样本采集、培养、QC，
不可直接命名为 acquired resistance。OCM74 的同 study 两个样本很近，但跨 study 比较
接近 matched between-patient reference，提醒 cohort/batch 仍可能影响相似性。

17位患者的 tumour–stroma 配对中，EMT 与 angiogenesis 的 percentile-score 差值为
−0.04548、−0.03621；cholesterol homeostasis、MYC targets V2、oxidative phosphorylation
为+0.01191、+0.01061、+0.00345；这些项目的 BH q 均约0.000455。
方向与 tumour/reference 谱系或增殖差异相容，但不是 HGSOC 特异机制验证。
OXPHOS 效应小且基于 RNA 排名，不能称为呼吸通量升高。
Hallmark 名称（如 PANCREAS_BETA_CELLS）只是 gene-set 标签，不证明细胞身份转变。
共50项检验中29项 q<0.05；有限 Monte Carlo resolution 与非独立 gene sets 均应注明。

完整数值来源与 receipt SHA 见 `evidence/post_audit_validation_audit_20260908.json`
及其指向的服务器 receipts。原1121180未产生 LP 结果；修复版1138529/1138550另行验证。
