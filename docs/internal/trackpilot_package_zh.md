# TrackPilot 工具包（2026-11-02）

## 1. 包里有什么

| 文件 | 作用 |
|---|---|
| `src/trackpilot/core.py` | 方法本身：轨道即假设、元数据锚点先验、家族先验、判别面板、后验平均预测 |
| `src/trackpilot/io.py` | 读预测表（长格式或宽格式）、读实测值、按关键词匹配轨道 |
| `src/trackpilot/cli.py` | 命令行：`suggest`、`predict`、`score-alphagenome` |
| `src/trackpilot/report.py` | 单文件 HTML 报告（无外部资源，支持深色模式） |
| `src/trackpilot/alphagenome.py` | 可选：用用户自己的 key 为一个区域做 ISM 打分（调用方式与论文所用脚本相同） |
| `tests/test_core.py` | 8 个单元测试（合成示例） |
| `tests/test_reproduce.py` | 与论文冻结结果逐数比对（需要本地数据，没有就跳过） |
| `pilot/scripts/trackpilot_reproduce.py` | 全量比对：每个元件、每个情景 |
| `examples/` | 合成示例（论文数据不能再分发） |

依赖只有 numpy、pandas、scipy。安装：`pip install -e .`。

## 2. 两条命令

- **`trackpilot suggest`**：一轮选出 k 个要测的变异（论文的 one-round 面板：每选一个，用当前后验平均预测代替它的标签）。
  - `--measured` 传入已测的标签时，在这些标签之后再选 k 个，可以分几轮做。
- **`trackpilot predict`**：读入实测值，输出：
  - 每个变异的预测：后验平均、最可能假设的预测、无标签默认值；
  - 轨道后验表；
  - HTML 报告，首行给出领先的轨道家族和它的后验占比。

## 3. 和论文结果是否一致

核心计算从 `proxy_select.py` 逐行移植，参数是冻结值（α=0.1，β=0.5，π=0.95，家族 r=0.9）。

比对方法：对每个元件、每个情景，取第 0 次重复，用包重新计算：
- 可测变异池；
- 顺序判别面板，以及 k=2,4,8,12,16 时留出变异上的 Spearman；
- 无标签参考：锚点、全轨道均值。

再与冻结运行写下的数逐个比较。

| 数据 | 元件数 | 比对的数字 | 最大差 |
|---|---|---|---|
| Kircher（v2 冻结运行 `proxy_select_v2/p0.95_fam`） | 22 | 234 | 0 |
| Siraj（`proxy_select_siraj_ref`） | 125 个位点 × 2 个细胞 = 250 | 4,500 | 0 |

- `test_reproduce.py` 另外在 3 个 Siraj 位点、前 3 次重复上比对一轮面板（batch_kb），也完全一致。
- 敏感性检查：把 π 从 0.95 改成 0.9，差值立刻变成 1e-3 量级，说明比对能抓到改动。

## 4. 速度

1,353 个变异 × 1,018 条轨道（IRF4），M2 Pro，后台同时在跑 Enformer：
- 读预测表约 2.4 s；
- 建立先验和轨道家族约 1.1–1.4 s；
- 选 8 个变异 0.04–0.08 s；
- 预测 0.005 s。

论文 Implementation 一节已写为"under five seconds on a laptop"。

## 5. 合成示例上的一个如实记录

在合成示例上，随机测 8 个变异和判别面板差不多，随机甚至略好（0.846 对 0.826）。原因是合成数据里各家族界限分明，测哪 8 个都能认出来。

所以没有把"判别面板优于随机"写成单元测试。这一点的证据来自论文的真实数据，包已逐数复现。示例不为此调整。

## 6. 还没做的

- 发布到 PyPI，在 GitHub 打 release，存档到 Zenodo（需作者决定版本号和时间）。
- `score-alphagenome` 只做了代码移植，没有用 key 实跑过（key 只在你的终端里）。可以在你的终端里跑一个小区域试一下，例如：

  ```
  trackpilot score-alphagenome --region chr11:5227000-5227050 -o test.tsv.gz
  ```
- 论文 §3.7 的使用示例图，可以直接用 `trackpilot predict --report` 的输出来画。
