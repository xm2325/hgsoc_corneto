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

1121178 与 1121182 目前只有完成日志已读；本文件不填入尚未审计的 pathway 名称、效应值或 p 值。
1121180 在 gene-ID 检查阶段失败，没有产生可解释的新 LP 结果。
