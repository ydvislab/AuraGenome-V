# A1 方法 v2：Siraj SatMut 预注册（已冻结，2026-10-29）

**状态：已冻结。** 本文件与评估代码（`proxy_select.py`、`siraj_report.py`）在同一个 commit 中提交。
- 冻结时 AlphaGenome 打分还在进行，只完成了 5 个位点；
- 冻结之后不再修改方法、参数、假设、界值或判定规则；
- 如果之后发现代码错误必须修正，修正内容和理由会另行记录，原判定照样报告。

**数据接触记录**：
- 2026-10-29 09:05，经作者同意下载 `data_satmut.zip`（Zenodo 15297965，md5 0043848d…，CC-BY-NC-4.0），只存放在本地被 git 忽略的 `pilot/data/raw/siraj2026/`。
- 至今只读取过结构性的列：元件、背景、细胞类型、位置、碱基。以及只为把实测值抄进变异表而读过效应列，没有对效应做过任何统计或作图。
- 空跑：在已打分的 5 个位点上，用随机生成的假标签（`numpy` seed 123）运行了完整流程，只用于确认代码能跑通。
  - 空跑发现，`--scenarios` 的默认值漏掉了 C2，已修正。
  - 修正后再次核对，Kircher 的输出仍与改动前逐行一致。

## 1. 这次检验回答什么

按 Codex 评审的说法，Siraj 检验的是"**新位点上的外部有效性**"，不是"模型没覆盖的细胞类型"：K562 和 HepG2 都在 AlphaGenome 的训练数据里。

"模型没覆盖的细胞"这一情形，只能通过藏起轨道来近似（条件 C3），论文中会明确写出这个局限。

## 2. 数据与分析单位

- **元件**：作者最终发布的 128 个元件（`satmut_final_128elements.txt`）。每个是 200 bp 的寡核苷酸，所有单核苷酸替换都在 K562 和 HepG2 中测过。
- **主背景**：
  - 每个元件只用参考等位基因背景（R；双变异单倍型元件用 RR）。原因是这样替换才能作为普通的 GRCh38 SNV 在基因组背景下打分。
  - 替代等位基因背景（A、AA 等）只用于敏感性分析，不进入主检验。
- **坐标**：
  - 数据是 GRCh37 坐标。精细定位的变异先用 Ensembl REST 映射到 GRCh38；
  - 然后在 GRCh38 中逐字匹配 200 bp 的 m0 序列（正链或反向互补）；
  - 匹配不唯一或找不到的元件剔除，并报告剔除数。脚本见 `pilot/scripts/siraj_prepare.py`。
- **映射结果**（10/29，只用了序列，没有用效应）：
  - 128 个元件中，有 1 个的参考背景本身含缺失（6:160560880，CATG>C 加 CTGGTAAGT>C），没有单核苷酸意义上的参考序列，剔除；
  - 其余 127 个中，125 个在 GRCh38 中找到唯一的精确匹配，全部在正链；
  - 另外 2 个找不到：14:106495830（IGH 区域）和 X:148621334，剔除。
  - **最终 125 个位点**，共 74,975 个替换，每个替换在 K562 和 HepG2 中各有一个实测值。
- **实测效应**：
  - `log2Skew`，即突变序列相对背景序列的 log2 比值；标准误 `log2Skew_SE`。v2 不使用标准误（`--noise` 关闭）。
  - 已知局限：`log2Skew` 的估计在作者流程中用到了整个文库的归一化，这是一个 8 变异实验得不到的信息。不过评估只用元件内的秩相关，归一化的整体尺度不影响结果。
- **分析单位：位点（元件）。** 同一位点的两个细胞类型先平均，得到 D_i。每个细胞类型的结果作为补充报告。
- **预测**：AlphaGenome ALL_FOLDS，ISM，1 Mb 上下文，所有生物样本的 CAGE、DNase 和 ATAC 轨道，与 Kircher 相同。

## 3. 冻结的方法

v2 冻结于 commit d389586（`docs/internal/method_v2_zh.md`）：

```
proxy_select.py --common-eval --anchor-prior 0.95 --family-prior --controls
```

- α = 0.1，β = 0.5，锚点斜率 sd 系数 0.5；
- 可测变异池 50%，200 次重复，seed 0；
- 主输出：k = 8 时的 rho_bma，即池外变异上的 Spearman。

**适配 Siraj 只做三项改动，都不改变方法本身：**
1. 读取 Siraj 的变异表；
2. 同一位点的两个细胞类型共用一份预测；
3. 实现下面的条件 C2：藏起本细胞系的轨道，用相关细胞的轨道平均作锚点（场景名 `hide_exact`）。C3 对应场景 `hide_matched`，此时本细胞系和相关细胞的轨道都藏起。
   - 已实现于 `proxy_select.py` 的 `--variants`；`related_biosamples` 一列只在 Siraj 表中存在。
   - 已核对：对 Kircher，v1 和 v2 的输出都与改动前逐行一致。

## 4. 三种条件（都事先定义，不看标签）

| 条件 | 保留的轨道 | 锚点 | 对应的现实场景 |
|---|---|---|---|
| C1 保留 | 全部 | 本细胞系轨道的平均（K562 或 HepG2，与 Kircher 的匹配规则相同） | 细胞类型被模型覆盖 |
| C2 藏起本细胞系 | 去掉本细胞系的轨道 | 相关细胞轨道的平均（下表） | 只有相近细胞被覆盖 |
| C3 广泛剥夺 | 去掉本细胞系和相关细胞的轨道 | 无，此时 v2 = v1 + 家族先验 | 模型没覆盖、也没有相近细胞 |

相关细胞的名单写于看任何 Siraj 结果之前，只根据轨道元数据（生物样本名，精确匹配）确定：
- **K562**（慢性髓系白血病，红系/髓系）：common myeloid progenitor, CD34-positive；hematopoietic multipotent progenitor cell；HL-60；myelocyte；bone marrow；bone marrow cell。
- **HepG2**（肝癌）：liver；left lobe of liver；right lobe of liver；hepatocyte；HuH-7；HuH-7.5。

## 5. 预注册的假设

记 D_i = 位点 i 上 v2（disc 选点，后验加权预测）减去参照，取两个细胞类型的平均。所有比较都在同一组池外变异上进行。

| 编号 | 条件，k | 假设 | 判定标准 |
|---|---|---|---|
| **P1 非劣效** | C1，k = 8 | v2 不比原样使用本细胞系轨道平均差 | D_i 均值的单侧 95% 下界大于 −δ，δ = 0.02（按位点自助法重抽 10,000 次，seed 0） |
| **P2 损害受控** | C1，k = 8 | 明显变差的位点很少 | D_i < −τ（τ = 0.05）的位点比例，其单侧 95% Clopper–Pearson 上界 ≤ q = 0.15 |
| **P3 救回（无锚点）** | C3，k = 8 | v2 优于不用标签的全轨道平均（在新位点上复现 Kircher 的 H2） | 平均差 ≥ 0.03，且单侧 Wilcoxon 符号秩检验 p < 0.05 |
| **P4 救回（相近代理）** | C2，k = 16 | v2 优于原样使用相关细胞轨道平均 | 平均差 > 0，且单侧 Wilcoxon p < 0.05 |

**次要的，只报告效应量和置信区间，不作为判定依据：**
- **S1 选点价值**：C1 和 C3、k = 8 时，v2 的 disc 选点对 v2 的随机选点，推断方式相同。
- **S2**：v2 对同一 disc 选点上的 `meta_ridge`（双侧）。开发历史显示两者打平，预期仍打平。
- **S3**：`matched_cal`（8 个标签的仿射校正）对原样匹配平均。预期有害。
- **S4**：k = 16 时的 P1、P2。
- **S5**：家族先验对均匀先验的消融。
- **S6**：替代等位基因背景的敏感性分析。需要额外打分，不在这一次运行中，结果另行报告。
- **S7**：每个细胞类型分开报告。

停止规则（disc_stop）在 v2 下已知会因先验惯性提前停止，不纳入本次检验。

**δ、τ、q 的理由**（定于看结果之前）：
- δ = 0.02，约为 Kircher 检验集上 H1 效应（0.08）的四分之一；
- τ = 0.05，沿用评审对"明显变差"的界线；
- q = 0.15，对应"每 7 个位点里最多约 1 个明显变差"。

## 6. 判定规则

- **P1、P2 成立，且 P3 或 P4 成立**："安全网"主张成立。论文主张为：
  - 被覆盖时，不会损害现成的代理；
  - 未被覆盖或只有相近细胞时，约 8 到 16 个实测值能带来改进。
- **P1、P2 成立，但 P3 和 P4 都不成立**：只主张"带元数据锚点的推断是安全的"，以及选点有价值（S1）。放弃"救回"的说法。
- **P1 或 P2 不成立**：按评审的停止条件，放弃"安全网"主张。论文重新定位为选点策略本身，或者"少标签下的代理校正风险"（例如 S3 的翻转现象）。

## 7. 执行顺序

1. ~~映射坐标，生成变异表（不看标签）。`siraj_prepare.py`。~~ 已完成：125 个位点。
2. 作者本人在自己的终端运行 AlphaGenome 打分，共 125 次 ISM 调用。API 密钥只经过作者的终端。
   - 命令：`python pilot/scripts/alphagenome_score.py pilot --variants data/siraj/variants_grch38.tsv --out predictions/siraj_all --versions ALL_FOLDS --outputs CAGE,DNASE,ATAC --all-biosamples`
   - 输出 `pilot/predictions/siraj_all/` 被 git 忽略。
3. ~~适配评估代码，用随机生成的假标签做一次空跑，确认代码能跑通。~~ 已完成。
4. ~~把本文件改为"冻结"，与代码一起提交、推送。~~ 即本 commit。
5. 打分全部完成后运行一次，按第 6 节判定，并原样报告。命令如下，每个位点一个进程：
   ```
   python pilot/scripts/proxy_select.py pilot --variants data/siraj/proxy_variants.tsv --pred-dir predictions/siraj_all \
       --common-eval --anchor-prior 0.95 --family-prior --controls --elements <locus>.K562,<locus>.HEPG2 --tag siraj/<locus>
   python pilot/scripts/proxy_select.py ...（同上，去掉 --family-prior，加 --scenarios all，--tag siraj_uni/<locus>）   # S5
   python pilot/scripts/siraj_report.py pilot/results/proxy_select_siraj results/a1
   ```
   判定写入 `results/a1/siraj_decision.txt`，逐位点的 D_i 写入 `siraj_by_locus.tsv`，次要结果写入 `siraj_secondary.tsv`。这些都是派生的汇总表，不含 Siraj 的原始数据。
