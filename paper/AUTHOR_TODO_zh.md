# 作者待办清单（2026-11-16）

`main.tex` 里还剩的红色 TODO 都需要作者本人决定或提供信息，我不能代填。下面每项都给了现状、需要作者做的事，以及可以直接改用的英文草稿。草稿只是起点：每一句都要作者确认属实，并按期刊 LLM 政策用自己的话改写。

> 篇幅提醒：正文现在正好 7 页（Bioinformatics Original Paper 的上限），最后一栏只剩约 4 行。作者信息、致谢、CRediT、资助补上之后很可能超页。到时候请告诉我，我再从正文里删几行。

---

## 1. 作者信息（`main.tex` 第 28–32 行、摘要里的 Contact）

需要提供：每位作者的姓名、单位（院系、机构、国家）、通讯作者及其邮箱，以及 ORCID（投稿系统会要求填写）。

## 2. 仓库名、PyPI 名称、Zenodo DOI（摘要 Availability 一段）

- **仓库名**：现在是 `ydvislab/AuraGenome-V`，和工具名 TrackPilot 对不上。建议在 GitHub 上把仓库改名为 `TrackPilot`，旧链接会自动跳转。改名后告诉我，我来统一替换正文、补充材料、README 和中文版里的链接。
- **PyPI**：`pyproject.toml` 里的包名是 `trackpilot`，版本 0.1.0。2026-11-16 查询时 PyPI 上还没有这个名字。发布要用作者自己的 PyPI 账号（`python -m build` 之后 `twine upload dist/*`），这一步我不替你操作。
- **Zenodo DOI**：在 Zenodo 上关联这个 GitHub 仓库，然后打一个 release（如 `v0.1.0`），Zenodo 会自动存档并分配 DOI。需要作者登录 Zenodo 授权。

拿到这些以后，Availability 一段可以写成：

> TrackPilot is a Python package and command-line tool (MIT licence), installable with `pip install trackpilot`, at https://github.com/<org>/TrackPilot; the submitted version is archived at Zenodo (doi:10.5281/zenodo.XXXXXXX).

## 3. 致谢（Acknowledgements）

英文草稿（请按实际情况增删）：

> We thank the authors of the saturation mutagenesis datasets of Kircher et al., Siraj et al. and Kosicki et al. for making their data public, and Google DeepMind for access to the AlphaGenome API. [Other people or facilities to thank.]

## 4. AI 使用声明

期刊要求披露 AI 的使用。投稿前请对照 Bioinformatics / OUP 作者指南的最新条款，确认声明应放在致谢、方法还是投稿信里。下面的草稿只写了这个项目里实际用过的工具，作者需要核对、补全，并删掉不属实的部分：

> During this work the authors used Claude (Anthropic) to help write and run analysis code, draft text and translate it, and OpenAI models through Codex to review the manuscript and analyses. The authors checked the code and results, rewrote the text, and take full responsibility for the content.

其中 "rewrote the text" 只有在作者真正重写之后才能保留。

## 5. 作者贡献（CRediT）

为每位作者从 14 个标准角色中选择：Conceptualization, Data curation, Formal analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources, Software, Supervision, Validation, Visualization, Writing – original draft, Writing – review & editing。

格式示例：

> A.B.: Conceptualization, Methodology, Software, Writing – original draft. C.D.: Supervision, Funding acquisition, Writing – review & editing.

## 6. 资助（Funding）

填写资助机构和项目编号。若没有专项资助：

> This work received no specific funding.

---

## 已完成（2026-11-16，本次不需要作者处理）

- Kircher 数据登录号：GEO GSE126550，已核对原文的 Data availability 声明。正文写明我们用的是作者 GitHub 仓库（kircherlab/MPRA_SaturationMutagenesis）发布的效应估计，下载脚本里固定了 commit。
- 图 1 重画为两栏：(a) 问题，用 IRF4/SK-MEL-28 的真实数据，由 `paper/make_fig1.py` 生成；(b) 流程。
- “实现”一节的 TODO：命令行参数已和代码（`trackpilot --help`）核对，二者一致。
