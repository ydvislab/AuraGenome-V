# 第三个数据集：iPSC 诱导神经元 MPRA（预注册，2026-11-03）

本文件和分析代码在任何预测打分之前写好并提交。整理数据时只读了位置、序列和轨道元数据，没有把任何效应和预测放在一起看。

## 1. 要回答的问题

前两个数据集的"无匹配轨道"都是把 K562、HepG2 的轨道藏起来模拟的。这次换成一种模型几乎没有覆盖的细胞：
- 细胞：WTC11 iPSC 经 NGN2 诱导的兴奋性神经元；
- AlphaGenome 里与之对应的轨道只有 4 条（neuron CAGE ×2、glutamatergic neuron DNase、motor neuron ATAC）；
- Enformer 里只有 2 条（GM23338 诱导的 bipolar neuron DNase、原代 Neurons CAGE）。

**问题：参数一个不改，TrackPilot 在这种细胞上是否仍满足前两次预注册的结论？**

## 2. 数据

- 来源：Kosicki et al., *Nat Commun* 16:4786 (2025)，ENCODE ENCSR548AQS，CC BY 4.0。
- 文件：ENCFF958ZLG（三个重复合并的元件定量）、ENCFF752NFF（设计序列），md5 已与 ENCODE 核对。
- 设计：270 bp 片段，每隔 4 个碱基引入一个人工颠换突变。
- 整理：`pilot/scripts/kosicki_prepare.py`，规则如下。
  - 效应 = 突变片段的 ratio − 同一实验中参考片段的 ratio。
  - 只用"+"链片段。"−"链片段只有 172 个突变，而且位置与参考序列对不上。
  - 每个突变都用设计序列核对参考碱基。
  - 只保留至少有 30 个定量突变的片段。
- 结果：328 个片段，归并为 159 个区域，22,397 个突变测量，20,486 个不重复的突变。三个抽查片段的序列与 hg38 完全一致。
- 分析单位：**区域**。重叠的片段共享同一组预测，先在区域内取平均，与 Kircher 的做法相同。

## 3. 轨道划分（只用元数据）

- **匹配轨道**（关键词，按"生物样本名 + 轨道名"做子串匹配）：
  - AlphaGenome：`CL:0000540`（neuron）、`glutamatergic neuron`、`motor neuron`；
  - Enformer：`bipolar neuron`、`Neurons,`。
- **相关轨道**（精确名称）：
  - 纳入：名称含 neur、brain、cortex、cerebel、hippocamp、WTC11、glutamat、spinal cord、medulla、putamen、caudate、thalam、pons、substantia 之一；
  - 排除：名称含 renal、kidney、adrenal、astrocyte、pericyte、endothelial、perineur、smooth muscle、crest、glioblastoma、glioma、neuroblastoma、tumor、neuroepithelioma、neurofibroma 之一；
  - 匹配轨道不计入相关轨道。
  - 结果：AlphaGenome 有 17 个相关生物样本，Enformer 有 70 个。
- **三种情景**（与前文相同）：
  - C1：保留全部轨道，锚点为匹配轨道均值；
  - C2：藏起匹配轨道，锚点为相关轨道均值；
  - C3：藏起匹配和相关轨道，无锚点。

## 4. 冻结的方法与运行

```
proxy_select.py pilot --variants data/kosicki/proxy_variants_<model>.tsv --pred-dir predictions/kosicki_<model> \
  --version <ALL_FOLDS|ENFORMER> --common-eval --anchor-prior 0.95 --family-prior --controls --tag kosicki_<model>/<region>
```

- 参数：α = 0.1，β = 0.5，π = 0.95，家族 r = 0.9；可测池 50%；200 次重复；seed 0。
- 汇总：`pilot/scripts/kosicki_report.py`。
- 打分：
  - Enformer：本机运行，`enformer_score.py`，打分方式与 Siraj 相同；
  - AlphaGenome（论文主模型）：用户在自己的终端用自己的 key 运行，打分方式与论文相同（CenterMaskScorer，501 bp，DIFF_LOG2_SUM）。
  - 两个模型分别检验、分别报告。

## 5. 假设（每个模型分别检验，单位为区域，k = 8）

界值与 P1–P4 相同。只有一处调整：每个片段约 68 个突变，可测池约 34 个，所以 N4 用 k = 8，而不是 P4 的 k = 16。

- **N1**（C3）：TrackPilot 减全轨道均值。平均 ≥ 0.03，且单侧 Wilcoxon P < 0.05。
- **N2**（C3）：判别面板减随机面板。单侧 Wilcoxon P < 0.05。
- **N3**（C1）：TrackPilot 减匹配均值。单侧 95% bootstrap 下界 > −0.02，且差值 < −0.05 的区域比例的 Clopper–Pearson 单侧上界 ≤ 15%。
- **N4**（C2）：TrackPilot 减相关轨道均值。单侧 Wilcoxon P < 0.05。

**判定**：N1 和 N3 都成立，判为"在模型几乎未覆盖的细胞类型上成立"。否则如实报告哪一条没有成立。

**描述性报告**（不做检验）：
- 匹配均值、全轨道均值、事后最佳单轨的绝对水平；
- 标签效率（k = 2、4、8、12、16）；
- 增益与默认差距的关系。

## 6. 局限（预先写明）

- 片段只有 270 bp，每个片段约 68 个突变，单个片段内的秩相关噪声较大。
- 人工颠换突变每隔 4 个碱基一个，不是饱和突变。
- 这组数据的匹配轨道不是同一种细胞，只是"神经元"大类。这正是要检验的情形，但匹配均值本身可能就很差。
- 分析单位是区域。同一区域内的重叠片段不独立，先取平均。
