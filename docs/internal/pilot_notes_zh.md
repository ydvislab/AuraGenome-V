# 试点：可靠性信号能否识别基因组 AI 的错误预测（研究方案第 3–6 周）

> 2026-10-28 准备完成。对应 [研究方案](research_plan_zh.md) 的 R1 和第 6 周决策点。
> **当前状态**：数据已就绪；三个脚本里只有数据准备和信号计算真正跑过，后者用的是模拟数据。**AlphaGenome 和 Borzoi 的脚本还没在真实 API 或 GPU 上运行过**（见第 4 节）。

## 1. 试点要回答的问题

> 在有实验真值的变异中，模型自身的可靠性信号（多个折版本之间的分歧、同一细胞类型不同实验之间的分歧、跨模型分歧、效应大小）能否区分"方向预测对了"和"方向预测错了"？

**决策规则（预先设定，结果出来后不再修改）**：同时满足以下两条为 **GO**，否则为 **NO-GO**，需要调整方案（例如改为以证据对照为主）：
- 组合信号的 AUROC 的 95% 置信区间下界 ≥ 0.6。置信区间按元件做 bootstrap，因为同一元件内的变异彼此不独立。
- 可靠性最高三分之一与最低三分之一的方向准确率相差 ≥ 15 个百分点。

## 2. 数据（已下载并整理）

- **来源**：Kircher et al. 2019 饱和突变 MPRA，数据取自门户官方仓库 `kircherlab/MPRA_SaturationMutagenesis` 中的 `data/elements.tsv.gz`。使用这些数据时必须引用 Kircher M, et al. *Nat Commun* 10:3583 (2019)。
- **`data/variants_grch38.tsv`**：共 39,170 个 SNV，覆盖 29 个元件，坐标为 GRCh38。
  - 每个变异附有实测效应（log2）、p 值和真值类别（+1 / −1 / 0）；
  - 每个变异附有 MPRA 所用细胞系、是否为人类细胞、在模型输出中匹配生物样本的关键词，以及该细胞系是否需要用替代细胞代表。
- **真值类别**：p < 1e-5 且 |效应| ≥ 0.1 时，按效应方向记为 +1 或 −1，其余记为 0。结果为 +1 有 2321 个，−1 有 4363 个，0 有 32,486 个。**这两个阈值沿用 CAGI5 风格，正式使用前需要对照 CAGI5 的 regulation-saturation 挑战说明核实。**
- **主分析只用"人类细胞且不需要替代细胞"的元件**，共 19,391 个 SNV、14 个实验。用替代细胞（如 HEL→K562）的元件和小鼠细胞中测得的元件留作次要分析。
- **剔除**：源文件中 11 行没有 GRCh38 坐标（MSMB 7 行、MYCrs11986220 4 行），已剔除；单碱基缺失另存为 `data/deletions_grch38.tsv`。
- `data/element_summary.tsv`：每个元件的变异数和显著变异数。

## 3. 脚本与运行顺序

| 步骤 | 脚本 | 需要什么 | 状态 |
|---|---|---|---|
| 1 | `scripts/build_variant_table.py` | 本地 Python | ✅ 已运行 |
| 2a | `scripts/alphagenome_score.py` | Python ≥ 3.10、`pip install alphagenome`、环境变量 `ALPHAGENOME_API_KEY`（非商业条款） | ⚠️ 只跑过 `--dry-run`：16 个独立区域 × 5 个模型版本，共 80 次 ISM 调用（每次再按 100 bp 分块） |
| 2b | `scripts/borzoi_score.py` | Python ≥ 3.10、GPU、`pip install borzoi-pytorch pyfaidx`、hg38 FASTA、`targets_human.txt` | ⚠️ 只做过语法检查 |
| 3 | `scripts/compute_signals.py` | 本地 Python | ✅ 用 `--simulate` 跑通，结果写在 `results_simulated/`（数字无生物学意义） |

**建议的运行顺序：**

```bash
# 0. 先在一个元件、一个模型版本上小规模试跑（核对 API 的输出列和耗时）
python scripts/alphagenome_score.py . --elements PKLR-48h --versions ALL_FOLDS
# 1. AlphaGenome 全部 5 个版本
python scripts/alphagenome_score.py .
# 2. 先只用 AlphaGenome 计算信号，此时已经可以做出 GO/NO-GO 判断
python scripts/compute_signals.py .
# 3.（可选，需要 GPU）先用 --check 核对输出形状，再补 Borzoi，加入跨模型信号
python scripts/borzoi_score.py . --fasta hg38.fa --targets targets_human.txt --check
python scripts/borzoi_score.py . --fasta hg38.fa --targets targets_human.txt
python scripts/compute_signals.py .
```

## 4. 已核实与未核实

**已经对照官方源码核实（AlphaGenome SDK 0.10.0）：**
- `dna_client.create(api_key, model_version=...)`；
- `ModelVersion` 包括 `ALL_FOLDS` 和 `FOLD_0`–`FOLD_3`。这 4 个折版本是交叉验证模型，**可以直接提供"副本分歧"信号**；
- `score_ism_variants(interval, ism_interval, variant_scorers, organism, progress_bar=...)`；
- `variant_scorers.RECOMMENDED_VARIANT_SCORERS` 的键名，如 `CAGE`、`DNASE`；
- `tidy_scores` 的输出列：`variant_id`、`output_type`、`track_name`、`biosample_name`、`raw_score` 等；
- 坐标约定：`Interval` 是 0 起始、右端开区间；`Variant` 是 1 起始。

**未核实，需要在第一次真实运行时确认：**
- API 的配额和速度。官方文档没有写明速率限制；
- `tidy_scores` 中 `variant_id` 的确切格式。脚本假定为 `chr:pos:ref>alt`，如果不同，需要调整 `compute_signals.py` 中拼接 `vkey` 的那一行；
- 匹配生物样本时用的关键词（如 `HepG2`、`HEK293`）能否在 AlphaGenome 的轨道元数据中命中。第一次运行后要检查每个元件命中了多少条轨道；
- Borzoi 移植版的输出裁剪和分箱（默认裁掉两侧各 5120 个 bin，每个 bin 32 bp），用 `--check` 核对；
- 安装环境：这台机器只有 Python 3.9，装不上 AlphaGenome SDK（要求 ≥ 3.10）。需要换一台有 Python 3.10+ 的机器，或在本机另装一个新版 Python。

## 5. 方法上需要注意的地方（写论文时也要讨论）

1. **MPRA 是质粒上的测量，不在基因组原位**。模型预测的是基因组原位的效应，两者本来就有差距，所以整体相关性偏低是预期之内的。试点问的不是"模型准不准"，而是"**可靠性信号能不能指出哪里不准**"。
2. **用哪类输出来对照**：脚本默认启动子用 CAGE、增强子用 DNASE，可以通过 `--promoter-output`、`--enhancer-output` 修改。这个选择应在看结果之前固定下来，或者把各种选择都报告出来，避免挑选对自己有利的结果。
3. **同一区域的重复实验**（PKLR-24h/48h，SORT1 的三个版本，LDLR 的两个版本）共享同一组模型预测，统计时要按元件聚类，bootstrap 已经这样做了。
4. **组合信号目前只是对各信号的排名取平均，没有拟合**。正式版本如果要学习权重，需要按元件做交叉验证。

## 6. 输出位置
- 预测结果：`predictions/alphagenome/<元件>__<版本>.tsv.gz`、`predictions/borzoi/<元件>__rep<k>.tsv.gz`
- 信号与报告：`results/signals.tsv`、`results/r1_report.json`（含 AUROC、三分位准确率、按元件的 Spearman 相关、GO/NO-GO 判断）
