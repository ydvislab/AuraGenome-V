# 重新选题：生信 × 计算机 × 可视化的交叉问题（查新报告）

> **已被取代（2026-10-28）**：本文属于"可靠性画像"阶段，项目方向已改为 A1（代理轨道选择），见 [`docs/STATUS_zh.md`](../STATUS_zh.md)。本文仅供参考。

> 2026-10-28。起因：用户判断"分层比较环"创新力不够，要求推倒重来，找一个对生信和可视化**都是新问题**、规模稍大但切实可行的研究问题，并严谨查新。
> 本报告只做**选题与查新**，尚未开始任何新实现。所有"未找到"的结论都只代表本轮检索范围，第 6 节列出了还需要补查的地方。

---

## 1. 先承认：之前的方向为什么不够

| 之前的卖点 | 已有工作 | 结论 |
|---|---|---|
| 样本分层聚合、分组展示 | **GenomeSpy**（GigaScience 2024）已支持对上千样本做排序、筛选、分组、聚合，并记录操作溯源；你的汇总文档里也引用了它 | 核心思路已有，只是我们换成了环形 |
| 环形图上的跨轨比较 | Gleicher 的比较三策略、L'Yi 2023 的多视图与响应式设计模式、*Should I make it round?* | 可视化领域的经典问题，不新 |
| 原位展开、证据下钻 | Chromoscope、FuSViz、HiCognition 都有"概览 → 局部 → 证据" | 交互模式已有 |
| 显著性检验（TreeBH 等） | 统计方法本身是现成的 | 是工程组合，不是新问题 |

**根本原因**：之前的选题是"给已知问题换一种画法"，而不是找到一个**两个领域都尚未面对的新问题**。

---

## 2. 新选题的筛选标准

1. **生信侧**：领域内 2025–2026 年的文献明确承认这是问题，并且影响真实的分析决策；
2. **计算机 / 可视化侧**：它对应一个公认的、尚未解决的研究问题，不只是"做一个工具"；
3. **可验证**：存在客观的"事后真相"，能定量评估，不能只靠主观问卷；
4. **可行**：数据和模型可以公开获取，能衔接你已有的积累（多轨基因组可视化、LLM、AuraGenome）。

---

## 3. 候选方向与查新对比

| | **A. 基因组 AI 预测的"可校准信任"** | B. 智能体生信分析的可视化审计 | C. 大规模泛基因组图的变异可视化 |
|---|---|---|---|
| 生信侧是否承认 | **强**：见 4.2 节，多篇 2024–2026 文献 | 强：2026 年有综述呼吁"结论可重建、可质疑、可追责" | 强：综述称图基因组可视化"很快变得复杂、难以理解" |
| 可视化侧是否承认 | **强**：EuroVis 2026 短文提出"为 XAI 做可视化应以校准信任为目标"，Han & Schulz 指出"错置的信任在 VA 中基本被忽视" | 中：已有 AgentTrails、AgentLens、TraceView、AgentTrace 等通用智能体轨迹可视化 | 中：已有 PanVA 等工作 |
| 已有工具 | AlphaGenome Atlas GUI、EVEE、ESM1b 门户——**都只展示单个模型的分数，没有可信度线索** | 通用智能体可视化已经比较拥挤，生信专用的还没有 | Bandage、ODGI viz、SequenceTubeMap、PGR-TK、Panagram 等很多 |
| 有没有客观真相用于评估 | **有**：饱和基因编辑（SGE）、MPRA 饱和突变、精细定位的 eQTL，可以事后检验"该不该信" | 弱：需要人工构造带错误的分析 | 弱：多为主观任务 |
| 与你积累的衔接 | 多轨 + 大规模比较 + LLM 解释，衔接好 | LLM 智能体，衔接好 | 衔接一般 |
| 主要风险 | 需要机器学习基因组学的背景和合作者 | 容易被看成"通用智能体可视化的应用" | 工具多，难出新意 |

**推荐 A**。B 可以作为 A 中的一个子问题（LLM 解释是否会误导信任）。C 不推荐。

---

## 4. 推荐方向 A：详细论证

### 4.1 问题陈述

> 基因组 AI 模型（如 AlphaGenome、Borzoi、Evo 2）正在以前所未有的规模输出预测：每个变异在数千条功能轨道、数百种细胞类型上都有预测值，Atlas 覆盖约 90 亿个变异。**研究者如何判断某个预测在某个位点、某种细胞类型上可不可信？** 可视化如何把模型内部的可靠性信号、外部的实验证据、不同模型之间的分歧组织起来，让研究者的信任与预测的真实可靠性**相匹配**（即"校准的信任"）？

用大白话说：模型很强，但并非处处可信；现在的工具只给分数，不告诉你"这里的分数靠不靠谱"。我们要做的，是让人一眼看出"哪里能信、哪里要小心、为什么"。

### 4.2 生信侧的证据：问题真实存在，而且没有被解决

| 证据 | 来源 |
|---|---|
| AlphaGenome 对每个 1 Mb 序列预测 5930 条人类功能轨道；Atlas 为约 90 亿个变异预计算了分数，并提供交互界面。该界面展示的是预测、ISM 分数、归因和 motif，**没有提到不确定性或模型分歧** | AlphaGenome（bioRxiv 2025 / Nature 2026）；AlphaGenome Atlas（medRxiv 2026） |
| 5 个只有随机种子不同的 Basenji2 副本，在 eQTL 和个人基因组上**超过一半的情况方向不一致**；模型在参考序列上常常"错得很自信"，在带变异的序列上反而缺乏信心 | Bajwa et al., MLCB 2024（PMLR 240） |
| 集成模型的一致性与误差相关，归因的一致性能进一步区分误差，离训练数据越远误差越大；但与 PKLR 启动子 MPRA 数据的相关性只有 ρ≈0.23 | Trust-Aware Sequence-to-Function Modelling（bioRxiv 2026） |
| DNA 语言模型的 ISM 分数主要反映"序列可预测性"而非调控；不同模型的 top-100 元件**完全不重叠** | The Dark Regulome（arXiv 2026） |
| 不同可解释性方法对同一预测给出**相互矛盾**的解释，作者呼吁告别"挑好看的例子" | Position: Genomic Model Research Must Move Beyond Anecdotal Evaluation（arXiv 2026） |
| 模型能抓住位点上的大致空间模式，但预测不了个体间差异 | A modality gap in personal-genome prediction（bioRxiv 2026） |
| 变异效应预测器之间差异巨大，例如在无序区域，AlphaMissense 的灵敏度只有 29% | Deep learning tools predict variants in disordered regions with lower sensitivity（2025） |

**归纳**：机器学习基因组学社区已经找到了一批"可信度信号"（集成分歧、归因一致性、离训练分布的距离、与实验的一致性），但这些信号**只存在于论文的统计图里，没有进入研究者实际使用的工具**。EVEE 甚至用 LLM 为预测生成自然语言解释，这可能在没有校准的情况下**抬高**信任。

### 4.3 可视化侧的证据：对可视化而言也是新问题

- **校准信任是公认的开放问题**：EuroVis 2026 的一篇短文主张，可视化应让"不确定性、溯源和认知局限变得可操作"，并警告校准线索也可能变成"信任噱头"；Han & Schulz 指出可视分析中"错置的信任"基本被忽视。
- **缺少科学领域的试验场**：基因组有 SGE、MPRA 这类**事后可检验的真相**，可以客观测量"可视化是否让专家的信任更贴近真实可靠性"。（补查后修正：用 Brier 分数等方式、对照真实结果来测量依赖校准，在人机交互领域**已有先例**，例如逻辑题、临床实体抽取、皮肤镜诊断等研究。所以方法本身不新，新在把它用到"多证据可视化 + 科学专家 + 基因组 AI"这一场景。详见第 7 节。）
- **规模与结构独特**：变异 × 数千条轨道 × 数百种细胞类型 × 参考/突变等位基因 × 多个模型，再加上实验证据，是一个多维多方的比较问题，现有的 VIS4ML 工具没有处理过。
- **本轮检索中，没有找到**发表在 IEEE VIS、EuroVis 或 TVCG 上的、面向基因组序列模型（Enformer、Borzoi、AlphaGenome 等）的可视分析系统。

### 4.4 可以拆成的研究问题

- **RQ1（需求）**：变异解读和调控基因组学研究者，在决定"信不信某个模型预测"时依据什么、缺什么？通过形成性访谈加任务分析回答。
- **RQ2（表达）**：如何把三类证据组织成可导航的多尺度视图？
  - 三类证据：模型内部信号（集成分歧、归因一致性、离训练分布的距离）、外部实验证据、模型间分歧；
  - 多尺度：全基因组 → 位点 → 轨道/细胞类型 → 碱基级归因；
  - 核心可视化挑战是：在数千条轨道中，快速定位"预测与证据不一致"的地方。
- **RQ3（效果）**：在有真相的数据上，这种可视化能否让研究者的信任**更校准**？
  - 测量指标：判断与真实结果之间的 Brier 分数和区分度；
  - 对照条件：只看分数、看分数加 LLM 解释、看分数加我们的证据视图；
  - 这组对照可以直接检验"LLM 解释是否会误导信任"，与 AuraGenome 一脉相承。

### 4.5 可行性

- **模型和数据**：
  - AlphaGenome 有非商业 API 和预计算的 Atlas；Borzoi、Enformer 的权重公开，可以训练或下载多个副本来获得分歧信号；
  - 真相数据有 BRCA1 SGE（Findlay 2018）、CAGI5 调控区 MPRA 饱和突变、GTEx 精细定位的 eQTL（都是公开数据，具体访问条件待确认）。
- **你已有的可以复用**：多轨环形或线性的基因组可视化经验、原型框架、E3 的实验模式（任务库、计时、日志导出）、TreeBH 等统计严谨性的经验、LLM 系统经验。
- **需要补的**：
  - 一位机器学习基因组学或变异解读方向的合作者；
  - GPU 资源（推理或训练多个副本）；
  - 熟悉 SGE、MPRA 数据的人。

### 4.6 投稿路线（建议）

- **论文 1（工具 + 案例）**：投 *Bioinformatics* 或 *NAR*，形式是"面向基因组 AI 预测的证据可视化平台"。
- **论文 2（研究问题 + 用户实验）**：投 IEEE VIS 或 TVCG，主题是"在有可检验真相的科学领域里，可视化如何校准对 AI 的信任"。

---

## 5. 候选 B、C 为什么不作为首选

- **B（智能体生信分析的可视化审计）**：生信社区已有综述呼吁可追溯、可验证，问题是真实的。但通用的"智能体轨迹可视化"在 2025–2026 年已经很拥挤（AgentTrails、AgentLens、TraceView、AgentTrace），做成生信版本容易被认为增量不足，而且缺少客观真相。它更适合作为 A 中关于"LLM 解释"的子问题。
- **C（泛基因组图可视化）**：问题公认、规模大，但 Bandage、ODGI、SequenceTubeMap、PGR-TK、Panagram、PanVA 等工具很多，竞争激烈，而且主要是图布局和尺度问题，和你的积累衔接较弱。

---

## 6. 查新的局限与必须补查的地方（严谨起见）

1. **Trust-Aware Sequence-to-Function Modelling（bioRxiv 2026）的全文因限流没有读到**，需要确认它是否附带了交互界面。如果有，它就是直接竞品。
2. **AlphaGenome Atlas 的 GUI 只读了前 1/4 的页面**，需要读完全文，确认是否展示了不确定性或实验对照。
3. **EVEE 的界面**需要实际试用，确认它是否给出了可信度线索。
4. **IEEE VIS 2025/2026 和 EuroVis 2026 的完整论文列表**需要人工逐条翻一遍（搜索引擎对会议论文的覆盖不全）。
5. **Bioinformatics、NAR Web Server 2025–2026 期**需要检索"variant effect + visualization"类工具。
6. **人机交互/可视化的信任校准研究**里，有没有用"科学领域的真相数据"做评估的（例如医学影像），需要补查，用来界定方法论贡献。

建议用 1–2 周完成上述补查，同时与 2–3 位变异解读方向的研究者做形成性访谈，确认 RQ1。**补查完成之前，不要开始实现。**

---

## 参考（本轮检索）
- GenomeSpy（GigaScience 2024）: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11299109/
- AlphaGenome（bioRxiv）: https://www.biorxiv.org/content/10.1101/2025.06.25.661532.full.pdf ；PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12851941/
- AlphaGenome Atlas（medRxiv 2026）: https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1.full ；博客: https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/
- AlphaGenome 可视化库文档: https://www.alphagenomedocs.com/visualization_library_basics.html
- Bajwa et al. 2024（PMLR 240）: https://proceedings.mlr.press/v240/bajwa24a.html
- Trust-Aware Sequence-to-Function Modelling（bioRxiv 2026）: https://www.biorxiv.org/content/10.64898/2026.08.20.745945v1.full
- The Dark Regulome（arXiv 2026）: https://arxiv.org/abs/2606.06834
- Position: Beyond Anecdotal Evaluation of Interpretability（arXiv 2026）: https://arxiv.org/abs/2606.07607
- Personal-genome modality gap（bioRxiv 2026）: https://www.biorxiv.org/content/10.64898/2026.02.01.702969v1.full
- 无序区域变异预测器灵敏度: https://www.researchgate.net/publication/390730645
- EVEE（bioRxiv 2026）: https://www.biorxiv.org/content/10.64898/2026.04.10.717844v1
- DEGU 不确定性（npj AI 2025）: https://www.nature.com/articles/s44387-025-00053-3
- tangermeme（Nat Methods 2026）: https://www.nature.com/articles/s41592-026-03254-z
- EuroVis 2026 短文（校准信任）: https://www.diva-portal.org/smash/get/diva2:2075906/FULLTEXT01.pdf
- Han & Schulz（校准 VA 中的信任）: https://cs.au.dk/~hjschulz/pdfs/trustrex.pdf
- 可信智能体生信综述（2026）: https://pmc.ncbi.nlm.nih.gov/articles/PMC13565537
- AgentTrails（arXiv 2026）: https://arxiv.org/abs/2607.18816
- 泛基因组图综述: https://doi.org/10.3390/jgbg1010005

---

## 7. 补查结果（2026-10-28，第 6 节六项）

### 7.1 逐项结论

| # | 补查对象 | 结果 | 对新颖性的影响 |
|---|---|---|---|
| 1 | **Trust-Aware S2F / RegTrust-XAI**（bioRxiv 2026，昆士兰大学等） | 有一个 **Hugging Face 演示**（Gradio）：粘贴一条 DNA 序列，输出 K562 ATAC 的预测值、**可信度等级 A–D 徽章**（集成一致性 × 归因一致性）、离训练分布的标记，以及一张归因柱状图。只覆盖单条序列、单个细胞系、单一检测、自建的小型 CNN；**没有用户研究** | **最近的竞品**。"给单个预测显示可信度徽章"已经有人做了，不能再作为卖点 |
| 2 | **AlphaGenome Atlas**（medRxiv 2026） | 门户展示 AVI 分数、按特征分解的贡献、数千项分子预测、motif 图谱，另有一个交互式 HTML 仪表盘（285 个生物样本、17 种模态）。**没有不确定性，没有模型间分歧，界面里也没有与实验数据的对照**。作者自己建议"按区域或应用设置阈值"，并说分数"只能作为证据链的一部分" | 缺口确认：最大规模的基因组 AI 门户没有可信度线索 |
| 3 | **EVEE**（Goodfire 与 Mayo Clinic，2026） | 网页提供 420 万个 ClinVar 变异的分数、扰动谱和 **LLM 生成的机制解释**。**没有不确定性、校准或预测器分歧**。LLM 解释只用"LLM 当评委"打分（3.8/5），没有检查解释是否忠实于模型，也**没有用户研究**；作者说明解释只是"需要专家审阅的假设" | 缺口确认，并直接支持"LLM 解释可能抬高信任"这一子问题 |
| 4 | **IEEE VIS 2025**（309 篇全部筛过）/ **VIS 2026**（会议日程全文筛过）/ EuroVis 2026 | VIS 2025、2026 都**没有**面向基因组序列模型（Enformer、Borzoi、AlphaGenome、Evo 等）的工作。但出现了明显的相邻趋势：<br>- VIS 2026 **最佳论文 CellPrism**：针对"AI 虚拟细胞"扰动预测的可视分析，作者把"显示预测不确定性"列为未来工作；<br>- 同为 VIS 2026 的 Human-Guided Causal Knowledge Injection for Virtual Cells、LatentFlow（分子图神经网络）、MolExplain（工作坊）；<br>- VIS 2026 新设工作坊 *VAxAutoSci: Visual Analytics in the Age of Autonomous Science*；<br>- 信任方向：*"Trust Junk" Leads to Unjustified Support...*（VIS 2026）、*Do You "Trust" This Visualization? An Inventory to Measure Trust*（VIS 2026）、*Your Model Is Unfair...*（VIS 2025）。<br>**EuroVis 2026 的全文列表被反爬虫拦截，没有逐条核对** | 主题被可视化社区重视（最佳论文所在的领域），但**基因组序列模型 + 信任校准**这个交叉点仍是空白；同时说明竞争会很快到来 |
| 5 | **Bioinformatics / NAR 2025–2026 网络服务器** | 没有找到 2025–2026 年展示变异效应预测可信度的服务器。较早的 DeepFun（NAR 2021）有 motif 可视化，但没有可信度；Nature Genetics 2026 有一篇综述 | 缺口确认（检索覆盖有限，见 7.3） |
| 6 | **对照真实结果测量信任校准的研究** | 人机交互领域**已有**：logic puzzle 实验（AAAI，184 人）、临床实体抽取中用 Brier 分数测依赖校准（110 人，发现解释会让新手过度依赖）、医学领域的置信度辅助诊断（Tschandl 2020）、不确定性信息提升决策校准（arXiv 2309.10852） | **修正之前的说法**：用真实结果测校准的方法不新；新在专家级科学判断、多源证据可视化、基因组 AI 这一组合 |

### 7.2 修正后的新颖性边界

**已有的，不能当贡献：**
- 给单个序列预测显示可信度等级（RegTrust 演示）；
- 用真实结果、Brier 分数测量人对 AI 的依赖校准（人机交互已有方法）；
- 生物 AI 模型的可视分析这一大类（CellPrism 等）。

**目前仍未找到先例，可以作为贡献：**
1. **规模与多源性**：面向基因组基础模型的预测（变异 × 数千条轨道 × 数百种细胞类型 × 参考/突变等位基因），把三类证据放在同一个可导航的可视化里做"三角验证"，并在大规模下快速定位"预测与证据不一致"的地方：
   - 模型内部的可靠性信号；
   - 多个模型、多个副本之间的分歧；
   - 外部实验证据（SGE、MPRA、eQTL）。
2. **基因组专家的信任校准实证**：在有真实结果的变异集上，比较以下几种条件下研究者判断的校准程度。这是第一个在基因组 AI 场景下做的此类实验。
   - 只看分数；
   - 分数加单一可信度徽章（RegTrust 式）；
   - 分数加 LLM 叙述解释（EVEE 式）；
   - 分数加多证据可视化。
3. **LLM 解释与证据可视化的对比**：检验 LLM 生成的机制解释会不会让专家过度信任，可视化证据能不能纠正。目前 EVEE 等门户已经在用 LLM 解释，但没有人做过这项检验。

### 7.3 仍然存在的检索局限
- EuroVis 2026 全文列表、VIS 2026 海报的完整摘要没有逐条核对；
- Bioinformatics 和 NAR 的 2026 网络服务器专刊没有逐篇翻阅；
- AlphaGenome Atlas 门户本身没有实际试用，只读了论文和博客；
- RegTrust 的完整论文（方法和讨论部分）只读了摘要和代码仓库；
- 这些都不太可能推翻 7.2 的结论，但正式写作前应补齐。

### 7.4 结论

方向 A **可以继续**，但需要调整表述：从"首次展示可信度"改为"**在基因组基础模型的规模下做多源证据三角验证，并实证检验它能否校准专家信任**"。需要警惕的是：这个领域发展很快，RegTrust、CellPrism 的作者、DeepMind 都可能在 6–12 个月内往这里走。如果要做，应尽快完成形成性访谈和原型，抢先发出预印本。

### 7.5 补查新增的参考
- RegTrust-XAI 代码与演示: https://github.com/MujeebOnawole/RegTrust-XAI
- AlphaGenome Atlas 博客: https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/
- EVEE（Goodfire）: https://www.goodfire.com/research/evee-explaining-genetic-variants ；MCP: https://github.com/goodfire-ai/evee-mcp
- IEEE VIS 2025 论文列表（JSON）: https://ieeevis.org/year/2025/program/papers.json
- IEEE VIS 2026 日程: https://ieeevis.org/year/2026/program/content/
- CellPrism（arXiv 2608.01669）: https://arxiv.org/pdf/2608.01669
- DeepFun（NAR 2021）: https://academic.oup.com/nar/article/49/W1/W131/6287843
- 非编码变异深度学习综述（Nat Genet 2026）: https://www.nature.com/articles/s41588-026-02769-7
- Too Sure for Our Own Good（AAAI）: https://ojs.aaai.org/index.php/AAAI/article/view/38798/42760
- Novice Reliance Calibration（arXiv 2610.07800）: https://arxiv.org/pdf/2610.07800
- Using AI Uncertainty Quantification to Improve Human Decision-Making: https://arxiv.org/pdf/2309.10852
