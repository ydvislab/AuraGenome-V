# C1 稳健性分析：留出折检查点（规格写于运行之前，2026-10-29）

## 要回答的问题

AlphaGenome 的 ALL_FOLDS 检查点是从全部 4 个折蒸馏得到的，训练时见过所有基因组区域的轨道数据。FOLD_0 到 FOLD_3 这 4 个检查点各自留出了一折基因组区域。

审稿人可能会问：A1 的结论是否依赖于"模型训练时见过这个位点"？如果模型在留出的位点上明显变差，少标签的判别和校正可能更有用，也可能更没用。

做法：对每个 Kircher 元件，用**留出该元件所在折**的检查点重跑 A1（v1 和 v2），与 ALL_FOLDS 对比。

- 只换预测文件（`--version`），方法、参数、seed、评估集都不变；
- MPRA 标签从来没有进入过任何检查点的训练，这里检验的只是"轨道数据的训练暴露"。

## 折的分配（来自 AlphaGenome 的折定义，打分前确定）

| 检查点 | 元件 |
|---|---|
| FOLD_0 | F9，MSMB，SORT1，SORT1-flip，SORT1.2 |
| FOLD_1 | BCL11A，FOXE1，MYCrs11986220，MYCrs6983267 |
| FOLD_2 | GP1BA，HBB，HBG1，IRF4，PKLR-48h，TERT-GAa，TERT-GBM，TERT-GSc，TERT-HEK |
| FOLD_3 | HNF4A，IRF6，LDLR，LDLR.2 |

共 22 个实验，与统一评估集复核（`proxy_select_ce_*`）相同。PKLR-24h 也已打分，但和以前一样不进入分析。

## 运行

脚本 `pilot/scripts/run_c1.sh`，每个元件在两个检查点（ALL_FOLDS 和留出折）下各跑两套配置：
- **v1**（预注册的 A1）：`--common-eval`；
- **v2**（冻结版）：`--common-eval --anchor-prior 0.95 --family-prior --controls --scenarios all,hide_matched`。

ALL_FOLDS 也在同一个脚本里重跑，保证两边代码完全一致。v1 的 ALL_FOLDS 结果应与 `proxy_select_ce_*` 逐行一致，作为一致性核对。

## 看的量（分析单位：基因组区域，同一区域的实验先平均）

主要看 k = 8，k = 16 作为补充。每个量都报告 ALL_FOLDS 下的值、留出折下的值、两者之差，以及按区域自助法得到的 95% 区间。分开报告开发集（5 个区域）、检验集（11 个区域）和全部 16 个区域。

1. **不用标签的参照**（all 场景）：匹配平均、全轨道平均、最佳单条轨道的 Spearman。用来衡量"训练暴露"本身的影响有多大。
2. **v1，hide_matched 场景**（预注册的 H1、H2、H3）：
   - disc − 随机选点；
   - disc − 全轨道平均；
   - disc − 挑预测效应最大的变异。
3. **v2，all 场景**：disc − 匹配平均，以及变差超过 0.05 的区域数。
4. **v2，hide_matched 场景**：disc − 全轨道平均。

## 事先声明的判读规则

- **结论稳健**：在检验集上，换成留出折后，第 2、3 类对比的平均差符号都不变。另外，H1（disc − 随机）和 H2（disc − 全轨道平均）的单侧 Wilcoxon 检验仍为 p < 0.05，这里只作为描述性结果报告。
- **结论依赖训练暴露**：任何一个对比的符号翻转，或者 H1、H2 不再显著。这种情况会如实写进论文的局限。
- 第 1 类的变化只用于解释，不作判定。

本分析不改动任何方法参数，也不用于选择方法。

## 结果（2026-10-29 16:21，运行一次，原样报告）

- **运行**：16:01–16:21，88 个任务全部成功。
- **一致性核对**：v1 的 ALL_FOLDS 重跑与 `proxy_select_ce_*` 逐字节一致，22 个实验都一致。
- **`results/a1/c1_decision.txt` 原文**：

```
v1 disc - random: ALL_FOLDS +0.0839 -> held-out fold +0.0929 (change +0.0090, 95% CI [-0.0154, +0.0380]) sign kept
v1 disc - all_mean: ALL_FOLDS +0.0739 -> held-out fold +0.0540 (change -0.0199, 95% CI [-0.0507, +0.0135]) sign kept
v1 disc - high_effect: ALL_FOLDS +0.0314 -> held-out fold +0.0388 (change +0.0074, 95% CI [-0.0337, +0.0499]) sign kept
v2 all: disc - matched: ALL_FOLDS +0.0006 -> held-out fold -0.0055 (change -0.0061, 95% CI [-0.0172, +0.0034]) sign FLIPPED
v2 hide_matched: disc - all_mean: ALL_FOLDS +0.0767 -> held-out fold +0.0588 (change -0.0179, 95% CI [-0.0500, +0.0210]) sign kept
H1 (v1 disc - random) under ALL_FOLDS: mean +0.0839, 9/11 positive, one-sided Wilcoxon p = 0.002441
H1 (v1 disc - random) under held-out fold: mean +0.0929, 10/11 positive, one-sided Wilcoxon p = 0.001465
H2 (v1 disc - all_mean) under ALL_FOLDS: mean +0.0739, 8/11 positive, one-sided Wilcoxon p = 0.0415
H2 (v1 disc - all_mean) under held-out fold: mean +0.0540, 7/11 positive, one-sided Wilcoxon p = 0.06152
READING: DEPENDS ON TRAINING EXPOSURE - a contrast flips sign or H1/H2 lose p < 0.05.
```

**按事先声明的规则，判定为"结论依赖训练暴露"。** 触发判定的有两处：

1. **v2 在保留匹配轨道时相对匹配平均的差值变号**：+0.001 → −0.006。
   - 这个量在两个检查点下都接近 0；变差超过 0.05 的区域都只有 1 个（FOXE1，chr9），没有增加。
   - 运行后回看，用"符号是否翻转"来判断一个接近 0 的非劣效量，规则本身过严。不过判定照样按原规则报告。
2. **H2（disc − 全轨道平均）在留出折下不再显著**：平均差 +0.074 → +0.054，11 个区域中为正的从 8 个变为 7 个，p 从 0.042 变为 0.062。

**其他结果**（描述性）：
- **H1 更强了**：disc − 随机选点为 +0.093，10/11 个区域为正，p = 0.0015。选点的价值不依赖训练暴露。
- **绝对精度普遍下降约 0.035**：检验集上最佳单条轨道从 0.465 降到 0.428，v1 disc 从 0.370 降到 0.335。这说明训练暴露确实抬高了绝对精度，但大多数对比（方法之间的差值）变化不大。
- **H2 的减弱集中在三个区域**：IRF4（−0.107）、MSMB（−0.078）和 LDLR（−0.077）。HBB 区域反而增加了 0.107。
- **v2 相对匹配平均的变化集中在 MSMB 和 LDLR**（各约 −0.03）。留出折下，这两个区域的匹配平均反而更准（MSMB 0.17 → 0.22，LDLR 0.25 → 0.35），v2 能超过它的空间变小了。
- **k = 16**：v2 保留匹配轨道时，相对匹配平均的差值为 +0.018 → +0.001，变差的区域从 0 个变为 1 个；H1、H2 的对比都保持为正。

**对论文的含义**：
- 选点优于随机（H1）不依赖训练暴露，可以放心作为主张。
- 在 11 个 Kircher 区域上，"优于全轨道平均"（H2）在留出折下只是接近显著，需要写进局限。另一方面，Siraj 的 125 个位点上 P3（同类对比）成立，但那是 ALL_FOLDS。
- 如果要把"训练暴露"这一点在大样本上讲清楚，可以用留出折检查点给 Siraj 位点重新打分。这需要先确定每个 Siraj 位点属于哪个折，打分时间约 2.5 小时。
