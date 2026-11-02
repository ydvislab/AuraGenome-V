# 论文重新定位与工具化方案（2026-11-02）

## 1. 期刊要求（来源：Bioinformatics 官方 Author guidelines）

- **范围**：计算分子生物学和基因组生物信息学。重点是能推动该领域"in a significant manner"的新算法和数据库。
- **学科类别**（投稿时选择）：Genome analysis、Sequence analysis、Gene expression、Data and text mining 等十类。我们应选 **Genome analysis**。
- **文章类型**
  - **Original Paper**：上限 7 页（约 5,000 词，不含图）。可以介绍模型、算法、软件或数据库，必须使用真实生物数据。
  - **Application Note**：上限 4 页，只介绍软件、数据库或网络服务。软件名要出现在标题里，不能是小工具。
- **软件要求**（对新算法同样适用）
  - 提供自成一体、易用的实现，附测试数据和文档（安装、运行、全部参数）；
  - 源码放在稳定网址，提交的版本存档到 Zenodo 等；
  - 说明如何复现论文结果和图；
  - **软件名尽量放进标题**。
- **机器学习要求**
  - Methods 中要有专门小节，写明训练集、交叉验证和独立测试集；
  - **要说明同源序列如何处理**，即我们第 222 行的 todo，这一项是硬性要求；
  - 不能只报告留一法的结果。
- **AI 政策**（ISCB）
  - 允许：语法检查、翻译、找不一致、写代码和文档；
  - 不允许：用 LLM 根据提示起草论文正文、图表或参考文献；
  - 使用了 AI 要在 cover letter 和 Methods 或 Acknowledgements 中声明。

## 2. 同类论文怎么写

以 GraphyloVar 为例（Bioinformatics 42(7), btag426, 2026）：

- **命名**：标题是"工具名：做什么"；摘要 Results 第一句是"We present 工具名, a …"。
- **结构**：Introduction → Materials and methods（模型、训练、评分、基线、**Implementation**）→ Results（按问题分节）→ Discussion and conclusion。
- **图**：正文 6 张、表 1 张，图 1 是方法总览示意图。
- **数据可用性**：GitHub（MIT）+ Zenodo DOI；第三方数据注明来源。

方法类论文几乎都有名字，并以工具的形式发布。所以我们不是"只能这么写"，现在的稿子反而偏离了常见写法。

## 3. 现稿和惯例的差距

| 方面 | 现稿 | 期刊惯例 | 改法 |
|---|---|---|---|
| 标题 | 问句、无名字 | 名字：做什么 | 例："X: choosing and calibrating proxy tracks of genomic foundation models with a few measurements" |
| 摘要 Results | 先讲数学，后面一串 P 值 | 第一句介绍工具，然后 2–3 个核心数字 | 先介绍 X，再给三个结论：没有匹配轨道时提升多少、有匹配轨道时不更差、一轮挑完就够 |
| 引言 | 贡献写成三条发现 | 先说缺口，再说"我们提出 X" | 补充为什么需要通用方法（任何多输出模型都会遇到） |
| 方法 | 没有 Implementation | 有 Implementation 和软件说明 | 新增：输入、输出、命令行和 API、运行时间 |
| 结果 | 按研究历史排（v1→v2→预注册→审稿追问） | 按读者关心的问题排 | 见第 5 节 |
| 术语 | disc、v2、C1/C2/C3、anchor | 面向用户的叫法 | 例如：disc → discriminating panel；C1/C2/C3 → profiled / related-only / unprofiled |
| 语气 | "requested in review""In hindsight""S3b was added after review" | 不出现审稿过程 | 删掉；"post hoc"的标注保留 |
| 图 | 1 张 tikz 流程图、1 条学习曲线、1 张 Siraj 散点 | 4–6 张，含使用示例 | 新图 1 讲问题（同一位点上不同轨道意见相反）；新增"工具输出示例"图 |
| 局限 | 一大段 | 精简，聚焦关键几条 | 保留关键几条，其余移到补充材料 |

## 4. 计算机角度的通用方法：能说到哪一步

**可以诚实声称的**
- 问题是通用的：一个多任务预测器输出很多头，用户只有少量标签，要决定用哪个头或哪种加权组合。AlphaGenome、Borzoi、Enformer 都是这种情况。
- 方法的三个部件不依赖 AlphaGenome：
  - 把每个头当作"未知仿射校准"的假设，校准参数解析积分掉；
  - 用元数据先验做默认值；
  - 按头之间的分歧挑下一个要测的样本。
- 工具的输入是"变异 × 轨道"的预测矩阵加轨道元数据，任何模型都能用。

**目前不能声称的**
- "在多个模型上都有效"：只测了 AlphaGenome，Borzoi 打分脚本有但没跑（需要 GPU）。
- 计算机层面的大创新：各部件（Box–Hill 判别设计、模型平均、共轭积分）都是经典方法。新意在于组合方式，以及"大量强相关回归头 + 未知标定 + 元数据先验"这个设定。这对 Bioinformatics 足够，对计算机顶会不够。
- 新出现的相关工作 arXiv 2609.18622（2026-09）：在分类器上证明"认证哪个模型最好"需要大量标签。它和我们的区别是：分类与回归、认证最优与预测精度。正好说明我们为什么输出加权平均而不是硬选一个，可以引用。

**如果要把"通用"写进标题或摘要**：最小补充实验是在 Siraj 125 个位点上用 Borzoi 再跑一遍冻结的方法，不改任何参数。需要 GPU，按 4 个复制模型估计约 1–2 天。不跑的话，只写"模型无关的接口，在 AlphaGenome 上验证"。

## 5. 建议的新结构（7 页）

1. **Introduction**（约 0.8 页）：缺口 → 为什么按名字挑轨道不可靠 → 我们提出 X → 主要结果一句话。
2. **Methods**（约 1.8 页）
   - 2.1 Problem
   - 2.2 Calibrated track hypotheses
   - 2.3 Metadata prior
   - 2.4 Choosing variants
   - 2.5 Implementation（新）
   - 2.6 Data
   - 2.7 Evaluation and pre-registration（含 ML 小节和同源说明）
3. **Results**（约 3 页），按问题排：
   - 3.1 没有匹配轨道时，8 次测量能恢复多少？（Kircher C3 + Siraj C3/C2）
   - 3.2 有匹配轨道时，会不会变差？（Siraj P1/P2 + Kircher v2；v1 的教训一句话带过）
   - 3.3 结论是否依赖模型见过这些位点？（held-out 折）
   - 3.4 和其他做法比（随机、回归、CODA / MODEL SELECTOR、EIG），以及为什么输出平均值
   - 3.5 实用建议：多少标签、一轮挑完、方向判断用不带截距的平均
   - 3.6 使用示例：一个位点的完整工具输出（新图）
4. **Discussion**（约 0.8 页）：适用条件（需要至少有一条大致成比例的轨道）、局限、推广到其他模型。

## 6. 工具化需要做的事

现在的 `pilot/scripts/proxy_select.py` 是评估脚本（576 行，参数里混着几十个实验开关），不是给用户用的工具。需要拆成一个包：

```
x/                         # 名字待定
  core.py        后验计算、加权平均预测、挑选变异（从 proxy_select.py 抽出，不改算法）
  priors.py      元数据匹配、轨道聚类、先验
  io.py          读预测矩阵（TSV/parquet）、轨道元数据、实测标签
  alphagenome.py 可选：用用户自己的 API key 为一个区域做 ISM 打分
  cli.py         x suggest   → 输出要测的 k 个变异（一轮）
                 x predict   → 输入实测值，输出轨道后验、每个变异的预测、HTML 报告
tests/           单元测试 + "复现论文数字"的回归测试（固定种子，对比 results/a1）
examples/        Kircher 一个元件的完整示例（数据可公开再分发，需核对许可）
docs/            安装、参数说明
```

- **发布**：PyPI 上的 `pip install x`，GitHub 打 release，Zenodo DOI。
- **工作量**：核心抽取和命令行约 1 天；测试、文档、示例约 1 天；HTML 报告图约半天。
- **复现检查**：用新包重算 Siraj 的 P1–P4，要求与现有结果逐位一致，这样论文数字不需要改。
- **可选**：交互式可视化界面。这正好是之后 VIS 论文的内容，建议现在不放进这篇。

## 7. 需要作者决定的事

1. **文章类型**：Original Paper + 命名工具（推荐），还是缩成 4 页的 Application Note？推荐前者，因为预注册检验和对照实验是我们的主要分量，压到 4 页会丢掉。
2. **工具名**：候选 `TrackPick`、`ProxyCal`、`Calibrack`。名字会进标题。2026-11-02 查过，三个名字在 PyPI 上都没被占用；GitHub 上有少量同名的个人仓库。
3. **是否补跑 Borzoi**：决定论文里能不能写"在多个模型上验证"。
4. **正文由谁执笔**：按期刊 AI 政策，正文不能由 LLM 起草。我可以做的：
   - 结构调整、逐段要点提纲；
   - 代码、工具和图；
   - 数字核对、语法检查。

   最终英文句子要由作者自己写。前几轮正文里我写的部分，也需要作者改写一遍。
