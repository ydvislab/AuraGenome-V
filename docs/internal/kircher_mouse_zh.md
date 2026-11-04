# 第四个检验：Kircher 小鼠细胞实验（预注册，2026-11-04）

本文件、整理代码和分析代码在这些实验的任何预测打分之前写好并提交。这 6 个实验此前从未用 TrackPilot 分析过，也没有打过分。

## 1. 要回答的问题

前面"模型没覆盖的细胞"的证据有两类：
- Siraj 是把轨道藏起来模拟的；
- Kircher 的 5 个区域是事后分组分析（S13）。

这次用一组全新的数据检验：**细胞类型完全没有对应的人类轨道时，参数一个不改，TrackPilot 是否仍然有效？**

## 2. 数据

- 来源：Kircher et al. 2019，饱和突变 MPRA。人类调控元件在小鼠细胞系中测量。论文主分析一直只用人类细胞的 22 个实验，这 6 个从未使用。

| 实验 | 元件 | 细胞 | 区域 |
|---|---|---|---|
| RET | 增强子 | Neuro-2a（小鼠神经母细胞瘤） | RET |
| UC88 | 增强子 | Neuro-2a | UC88 |
| TCF7L2 | 增强子 | MIN6（小鼠胰岛 β 细胞瘤） | TCF7L2 |
| ZFAND3 | 增强子 | MIN6 | ZFAND3 |
| ZRSh-13、ZRSh-13h2 | ZRS 增强子 | NIH/3T3（小鼠成纤维细胞，共转 HOXD13 或 HOXD13 + HAND2） | ZRS（两个实验先平均） |

- 效应、坐标（GRCh38）与论文其余 Kircher 分析相同（`data/variants_grch38.tsv`）。
- 整理：`pilot/scripts/kircher_mouse_prepare.py`。
- 分析单位：**区域**，共 5 个。

## 3. 轨道划分（只用元数据）

人类轨道里没有这些小鼠细胞，所以按细胞的类型指定一个人类"替代"（proxy），在任何结果出来之前写定：

| 细胞 | AlphaGenome 替代 | Enformer 替代 |
|---|---|---|
| Neuro-2a | BE2C、SK-N-DZ、SK-N-SH（神经母细胞瘤） | 同左，另加名称含 "neuroblastoma cell line" 的 FANTOM 样本 |
| MIN6 | islet precursor cell、progenitor cell of endocrine pancreas | islet precursor cell |
| NIH/3T3 | 名称含 fibroblast 的全部样本 | 同左 |

- 匹配方式为子串匹配，不区分大小写，与论文其余部分相同。
- 不设"相关轨道"。
- 两种情景：
  - **C1（有替代）**：锚点 = 替代轨道均值；
  - **C3（无替代）**：藏起替代轨道，无锚点。

## 4. 冻结的方法与运行

```
proxy_select.py pilot --variants data/kircher_mouse/proxy_variants_<model>.tsv --pred-dir <dir> --version <ver> \
  --common-eval --anchor-prior 0.95 --family-prior --controls --elements <experiment> --tag kircher_mouse_<model>/<experiment>
```

- 参数：α = 0.1，β = 0.5，π = 0.95，家族先验 r = 0.9；可测池 50%；200 次重复；seed 0。
- 运行：`sh pilot/scripts/run_kircher_mouse.sh <alphagenome|enformer>`；汇总：`kircher_mouse_report.py`。
- 打分：
  - AlphaGenome（主模型）：用户在自己的终端用自己的 key 运行。命令为 `alphagenome_score.py pilot --elements RET,UC88,TCF7L2,ZFAND3,ZRSh-13,ZRSh-13h2 --include-mouse-cell-elements --versions ALL_FOLDS --outputs CAGE,DNASE,ATAC --all-biosamples`，打分方式与论文其余 Kircher 预测完全相同。
  - Enformer：本机运行，`enformer_score.py`，输出到 `predictions/enformer_kircher_mouse`。
- 两个模型分别检验、分别报告。

## 5. 假设（每个模型分别检验，单位为区域，k = 8）

只有 5 个区域，Wilcoxon 检验力不足。所以与 Kircher 原预注册的 H2 一样，用"平均差 + 方向一致"的标准。

- **M1（C3，主）**：TrackPilot 减全轨道均值。平均 ≥ 0.03，且 ≥ 4/5 个区域为正。
- **M2（C3）**：判别面板减随机面板（同一推断）。≥ 4/5 个区域为正。
- **M3（C1）**：TrackPilot 减替代轨道均值。平均 ≥ −0.02，且没有区域低于 −0.05。

**判定**：M1 和 M3 都成立，判为"在完全未覆盖的细胞上成立"。否则如实报告哪一条没有成立。

**描述性报告**（不做检验）：
- 各方法的绝对水平；
- 事后最佳单轨；
- 适用性检查（P > 0.5）判为"无信号"的比例。

## 6. 局限（预先写明）

- 只有 5 个区域，结论的统计强度有限。
- 人类序列在小鼠细胞里测量，反式环境是小鼠的，而轨道全部是人类的。这比"未覆盖的人类细胞"更难。
- ZRS 实验共转了转录因子（HOXD13、HAND2），成纤维细胞替代可能不贴切。
- 替代样本的选择没有经过验证，只依据细胞类型的名称。
