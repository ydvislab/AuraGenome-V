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

## 7. 结果（2026-11-04，按第 4–5 节冻结的方法运行，参数未改）

单位为区域（5 个），k = 8，200 次重复。输出：`results/a1/kircher_mouse_<model>/{decision.txt, by_region.tsv}`。

| 假设 | AlphaGenome（主模型） | Enformer |
|---|---|---|
| M1（C3）：TrackPilot − 全轨道均值，≥ 0.03 且 ≥ 4/5 为正 | −0.004，2/5 → **不成立** | +0.022，3/5 → **不成立** |
| M2（C3）：判别面板 − 随机面板，≥ 4/5 为正 | +0.086，4/5 → 成立 | +0.102，5/5 → 成立 |
| M3（C1）：TrackPilot − 替代均值，平均 ≥ −0.02 且无区域 < −0.05 | +0.025，最差 RET −0.068 → **不成立** | +0.017，最差 RET −0.044 → 成立 |
| **判定（M1 且 M3）** | **不支持** | **不支持** |

逐区域（Spearman ρ，k = 8；tp = TrackPilot）：

| 区域 | 模型 | tp C3 | 全轨道均值 | 随机 C3 | tp C1 | 替代均值 | 事后最佳单轨 |
|---|---|---|---|---|---|---|---|
| RET | AlphaGenome | 0.060 | 0.106 | 0.071 | 0.170 | 0.237 | 0.256 |
| TCF7L2 | AlphaGenome | 0.168 | 0.139 | 0.085 | 0.143 | 0.136 | 0.255 |
| UC88 | AlphaGenome | 0.318 | 0.334 | 0.189 | 0.332 | 0.327 | 0.356 |
| ZFAND3 | AlphaGenome | 0.451 | 0.426 | 0.341 | 0.407 | 0.240 | 0.513 |
| ZRS | AlphaGenome | 0.216 | 0.228 | 0.097 | 0.235 | 0.223 | 0.305 |
| RET | Enformer | 0.165 | 0.172 | 0.086 | 0.184 | 0.228 | 0.256 |
| TCF7L2 | Enformer | 0.117 | 0.105 | 0.041 | 0.115 | 0.115 | 0.204 |
| UC88 | Enformer | 0.202 | 0.185 | 0.095 | 0.214 | 0.213 | 0.283 |
| ZFAND3 | Enformer | 0.423 | 0.302 | 0.273 | 0.321 | 0.186 | 0.483 |
| ZRS | Enformer | 0.169 | 0.203 | 0.068 | 0.214 | 0.219 | 0.281 |

描述性结果（不做检验）：
- 平均水平（C3，k = 8）：AlphaGenome 的 TrackPilot 为 0.243，全轨道均值 0.246，随机面板 0.157；Enformer 分别为 0.215、0.193、0.113。
- 平均水平（C1，k = 8）：AlphaGenome 的 TrackPilot 为 0.257，替代均值 0.233；Enformer 分别为 0.210、0.192。
- 事后最佳单轨：AlphaGenome 0.337，Enformer 0.302，都远高于任何不看标签的默认。
- 标签更多时（事后）：C3 中 AlphaGenome 的 TrackPilot 在 k = 12、16 时为 0.258、0.265，超过全轨道均值；Enformer 为 0.228、0.233。

解读：
- 两个模型都没有达到预注册判定，原因相同：没有替代轨道时，8 个测量值不足以稳定地胜过全轨道均值（M1）。这组数据上全轨道均值本身接近替代轨道均值，留给方法的空间小。
- 判别面板在两个模型上都明显好于随机面板（M2），与前三个检验一致。
- 有替代轨道时，平均略好于替代（两个模型都 > 0）。替代选错的区域（ZFAND3，MIN6 用胰岛前体）提升最大（+0.17、+0.14）；替代本来就好的区域（RET）有损失，AlphaGenome 上超过了预设的 −0.05。
- 论文据此收窄"未覆盖细胞"的表述：方法能纠正选错的替代，并始终优于随机测量；但在完全没有对应轨道的细胞中，8 个标签不一定胜过全轨道均值。

### 适用性检查在小鼠实验上的表现（描述性，2026-11-04）

冻结的检查（p* = 0.5，k = 8，20 次重复，200 次置换），结果文件 `results/a1/kircher_mouse_<model>/check.txt`：

| 模型 | 情形 | 回退比例 | TrackPilot | 经检查 | 默认 |
|---|---|---|---|---|---|
| AlphaGenome | C1 | 0.05 | 0.262 | 0.266 | 0.232 |
| AlphaGenome | C3 | 0.21 | 0.245 | 0.246 | 0.246 |
| Enformer | C1 | 0.21 | 0.210 | 0.216 | 0.190 |
| Enformer | C3 | 0.29 | 0.210 | 0.229 | 0.193 |

（这里的重复数和面板与主检验不同，所以 TrackPilot 一列与上表的 200 次重复均值略有出入。）

解读：小鼠实验的测量值大多仍跟某些轨道相符，所以检查回退的比例（5–29%）远低于神经元数据（43–46%）。经检查后的结果不低于不加检查的 TrackPilot（差 +0.001 至 +0.019），但在 AlphaGenome C3 上仍只与全轨道均值持平。检查不能弥补 M1 的失败。

（第一次运行因脚本缺少 `__main__` 保护，子进程反复重启而挂起；修正后重跑，结果如上。）
