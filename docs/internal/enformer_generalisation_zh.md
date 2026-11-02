# 泛化检验：换成第二个基础模型 Enformer（规格，运行前冻结，2026-11-02）

**状态：已冻结。** 本文件与打分脚本 `pilot/scripts/enformer_score.py` 在同一个 commit 中提交，提交时还没有任何 Enformer 预测。

冻结之后：
- 不再修改方法、参数、假设、界值或判定规则；
- 如果之后发现代码错误必须修正，修正内容和理由另行记录，原判定照样报告。

## 1. 这次检验回答什么

TrackPilot 只用"变异 × 轨道"的预测表和轨道元数据，理论上不依赖 AlphaGenome。这次检验的问题是：**换一个基础模型，参数一个都不改，Siraj 预注册的结论是否依然成立？**

这次检验不回答"Enformer 和 AlphaGenome 哪个更好"，两者的绝对水平只作描述性报告。

## 2. 模型与打分

- **实现**：`enformer-pytorch` 0.8.12，权重 `EleutherAI/enformer-official-rough`（经作者同意于 2026-11-02 下载），human head，在本机 MPS 上运行。
- **轨道**：`targets_human.txt`（calico/basenji，cross2020）中的 CAGE（638）、DNASE（674）和 ATAC（10）轨道，共 1,322 条。不用 ChIP。
- **序列**：GRCh38，以元件中点为中心取 196,608 bp，通过 UCSC REST API 获取。每个变异都检查参考碱基一致。
- **打分**：
  - 单次前向，不做反向互补，也不做平移平均；
  - 变异分数 = log2(1 + Σalt) − log2(1 + Σref)，求和范围是与"以变异为中心的 501 bp"重叠的所有 128 bp 输出 bin。
  - 这对应 AlphaGenome 的 CenterMaskScorer(width=501, DIFF_LOG2_SUM)。
- **批量**：只用批量 1。实测 MPS 在批量大于 1 时输出错误（与批量 1 的最大差 983）；批量 1 与 CPU 一致（相关 0.99999999，最大差 0.005）。
- **耗时**：实测 MPS 约 2.4 s/序列，CPU 约 3.5 s/序列。Siraj 用 MPS 跑，Kircher 用 CPU 并行跑。
- **冒烟测试**：在 HBG1 的前 6 个变异上跑过一次，只看了预测值的分布，没有接触任何标签。

## 3. 数据

**Siraj**
- 标签相对参考寡核苷酸（`siraj_ref`），与论文主分析相同。
- 位点：按字母序排好的 125 个位点，用 `numpy.random.default_rng(0).choice(..., 40, replace=False)` 随机抽取 40 个。
- 抽 40 个的原因是算力：40 个约需 16 小时，全部 125 个约需 50 小时。抽样只用了位点名，没有看任何结果。
- 这 40 个位点：

10_121010256_A_G_wC 10_75928933_T_C_wC 10_81164146_C_T_wC 10_95535240_C_A_wC 11_3855868_C_T_wC_Alt_11_3855917_C
11_76887140_G_A_wC 11_9336093_C_T_wC 12_53440779_C_T_wC 14_65175364_A_C_wC 15_29077460_T_G_wC 15_83479349_T_C_wC
15_96714746_T_G_wC 17_31239645_C_T_wC 19_10347084_T_C_wC 19_1493207_C_T_wC 19_4153566_A_G_wC 19_53360976_C_G_wC
1_11713563_G_A_wC 1_211556131_C_G_wC 1_230849886_T_G_wC 1_86174772_T_C_wC_Alt_1_86174814_T_C 1_93426167_G_A_wC
20_45319505_T_A_wC 20_45944659_G_C_wC 22_42466980_A_T_wC 22_50742721_T_A_wC 2_224874874_G_A_wC 3_10328782_G_A_wC
3_42306446_C_T_wC 3_48632756_A_C_wC 5_150828014_G_A_wC 5_176731452_G_T_wC 5_52405693_C_A_wC
5_95163437_A_G_wL_Alt_5_95163449_G 6_127839843_A_C_wC 7_139875230_C_T_wC 7_19145326_T_C_wC_Alt_7_19145340_A_G
7_27717716_G_C_wC 7_77026460_G_A_wC 8_11720227_G_C_wC

**Kircher**
- 检验集的 17 个实验、11 个区域，和预注册 A1 一样排除开发集和 PKLR-24h。
- 共享同一区域的实验只打分一次。

**匹配与相关轨道**
- **匹配轨道**：规则同 AlphaGenome，在"生物样本名 + 轨道名"中做子串匹配。
  - Siraj：HepG2 有 4 条，K562 有 9 条；
  - Kircher：用变异表中原有的关键词。
- **相关轨道**：AlphaGenome 用的是精确名称列表，Enformer 的样本命名不同。这里用同一组关键词，在 Enformer 样本名中做子串匹配，并排除已匹配的轨道。匹配时只用了元数据。
  - HepG2：liver、hepatocyte、HuH-7，共 8 个样本名；
  - K562：CD34、common myeloid progenitor、hematopoietic multipotent progenitor、HL-60、myelocyte、bone marrow，共 43 个样本名。
  - 完整列表见本地的 `pilot/data/enformer/related_biosamples.tsv`，由本节规则生成。

## 4. 冻结的方法与运行

和 Siraj 主分析完全相同（`run_siraj.sh` 的主配置），只换预测目录、版本名和变异表：

```
proxy_select.py pilot --variants data/enformer/siraj_ref_proxy_variants.tsv --pred-dir predictions/enformer_siraj \
  --version ENFORMER --common-eval --anchor-prior 0.95 --family-prior --controls --elements <locus>.K562,<locus>.HEPG2
proxy_select.py pilot --pred-dir predictions/enformer_kircher --version ENFORMER --common-eval \
  --anchor-prior 0.95 --family-prior --controls --elements <17 个检验实验>
```

- α = 0.1，β = 0.5，π = 0.95；
- 可测变异池 50%，200 次重复，seed 0。

## 5. 假设与判定

**Siraj，40 个位点**（界值与原预注册相同）：
- **G1**（同 P1）：C1，k = 8，TrackPilot − 匹配轨道均值。单侧 95% bootstrap 下界 > −0.02。
- **G2**（同 P2）：C1，k = 8，D < −0.05 的位点比例。Clopper–Pearson 单侧 95% 上界 ≤ 15%。
- **G3**（同 P3）：C3，k = 8，TrackPilot − 全轨道均值。均值 ≥ 0.03，单侧 Wilcoxon P < 0.05。
- **G4**（同 P4）：C2，k = 16，TrackPilot − 相关轨道均值。单侧 Wilcoxon P < 0.05。
- 统计量用 `siraj_report.py` 计算，与原检验相同。

**判定**：G1、G2 成立，且 G3 或 G4 成立，即判为"在第二个模型上成立"。否则如实报告哪一条没有成立，论文中不写"模型无关"。

**Kircher，11 个区域**（次要）：
- **K1**（同 H1）：C3，k = 8，判别面板优于随机面板。单侧 Wilcoxon P < 0.05。
- **K2**（同 H2）：C3，k = 8，比全轨道均值高 ≥ 0.03，且至少 7/11 个区域为正。
- 用 `prereg_test.py` 计算。注意：原 H1/H2 用的是 v1；这里用冻结的 v2，在 C3 下没有锚点，两者只差家族先验。

**描述性报告**（不做检验）：
- 两个模型下匹配轨道均值、全轨道均值、事后最佳单轨的绝对水平；
- 在同一批位点上，TrackPilot 相对 label-free 默认值的增益。

## 6. 局限（预先写明）

- 只用了 Siraj 的 40 个位点，统计功效比 125 个低。
- Enformer 的训练数据同样覆盖这些区域，而且它没有留出折的 checkpoint，所以不做训练暴露检验。
- Enformer 没有在 MPRA 上训练，K562 和 HepG2 也都在它的训练轨道里。这次检验只回答"换模型"的问题，不回答"真正没被覆盖的细胞"的问题。

## 7. 结果：Kircher（2026-11-03）

打分 11 个区域、14,175 个变异，无报错。按第 4 节的冻结命令运行（`--tag enformer_kircher`），再用 `prereg_test.py pilot enformer_kircher` 检验。C3，k = 8，后验平均预测，按区域平均：

| 检验 | 结果 | 判定 |
|---|---|---|
| K1 判别面板 > 随机面板 | 平均差 +0.097，10/11 个区域为正，单侧 Wilcoxon P = 0.0010 | 成立 |
| K2 判别面板 − 全轨道均值 | 平均 +0.109，8/11 个区域为正 | 成立 |
| （次要）判别面板 > 高效应面板 | 平均差 +0.046，P = 0.0010 | 成立 |

- 判别面板低于全轨道均值的 3 个区域：LDLR（−0.011）、TERT（−0.031），以及 FOXE1（−0.013）。FOXE1 在 Enformer 下几乎无信号：事后最佳单轨也只有 0.085。
- 增益最大的是 HBG1（+0.463）和 SORT1（+0.315），这两个区域的全轨道均值很低。
- Siraj（G1–G4）还在打分，结果出来后补在下面。
