"""Build supp_tables_zh.tex, the Chinese version of supp_tables.tex, for supplement_zh.tex.

Run after make_supplement.py:  python paper/make_supplement_zh.py
Captions are replaced by label; table cells are translated by an exact-match dictionary plus a few substring
rules. Numbers are never retyped: the script checks that the multiset of numbers in every table is unchanged
(captions included), so if make_supplement.py changes a number, the matching Chinese caption must be updated here
and the script fails until it is.
"""
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "supp_tables.tex")
DST = os.path.join(HERE, "supp_tables_zh.tex")

CAPTIONS = {
    "s1": r"Siraj 检验（125 个位点）：效应相对实测的参照寡核苷酸计算（主要分析），以及按原文发表的方式以元件全部寡核苷酸中心化（原始分析）。P1--P4 为预注册假设，在两种定义下全部成立。次要分析行给出位点平均差值及其 95\% bootstrap 区间。S3b 是改用相对参照的标签时新增的。",
    "s2": r"方向准确率（事后分析）：留出的显著变异（$P<10^{-5}$，$|\text{效应}|\ge0.1$）中预测方向正确的比例。带标签的方法使用后验平均预测；“去偏移”表示去掉拟合的偏移（诊断用）。区间：按单位做 95\% bootstrap；最后一列为前者更好/更差的单位数。disc：试点变异组；random：随机组；matched：匹配轨道平均；all\_mean：全轨道平均；related：相近细胞平均；matched\_cal：对匹配轨道平均做仿射重拟合；@ 后为标签数。",
    "s3": r"为每个位点换用留出其基因组折的 AlphaGenome 检查点后，重做 Siraj 检验（标签相对参照寡核苷酸）。各列为蒸馏检查点（ALL\_FOLDS）和留出检查点下的预注册对比，以及配对变化量及其按位点的 95\% bootstrap 区间。",
    "s4": r"选点方式比较（事后分析；推断使用冻结的 v2，50 次重复）。Siraj 为场景 C1--C3；Kircher 为检验区域，分为保留匹配轨道（deploy）和没有锚点（no\_anchor）。disc：逐个选择的试点变异组；batch\_kb：测量前一轮选完；disc\_static：不随标签更新的先验分歧；eig：期望信息增益；high\_effect：预测效应最大；disagree：未加权分歧最大。判读：95\% 区间落在 $\pm0.02$ 以内为等价；整体高于或低于 0 为前者或后者更好。",
    "s5": r"CODA 式与 MODEL SELECTOR 式基线（事后分析；50 次重复）。两者都使用校准似然、不使用元数据、按信息增益选点；CODA 使用共识先验（$\kappa=10$，在开发元件上固定），MODEL SELECTOR 使用均匀先验（$\kappa=0$）。“单条轨道”：最可能的轨道，即这些方法本身的输出；“平均”：它们的后验平均。“+”：前者更好的单位数。",
    "s6a": r"信息增益（EIG）收敛检查（事后分析）之覆盖率：保留的假设所占的后验概率，按单位取中位数（括号内为最小值），分别给出第一步和 $k=8$ 变异组八步中的最低值。E0：最可能的 100 个假设，96 个积分点（表 S4 使用的设置）；E1：300 个假设；E2（未列出）每一步都保留 $\ge$99\% 的概率，使用 192 个积分点。",
    "s6b": r"信息增益收敛检查之准确度：后验平均预测的 Spearman 相关差值（相同的变异池和 50 次重复；按单位做 95\% bootstrap 区间）。",
    "s7": r"Siraj 上的标签效率（事后分析，描述性；125 个位点，200 次重复）。测了 $k$ 个试点或随机变异后后验平均预测的平均 Spearman 相关、不用标签的默认预测，以及用全部标签选出的最佳单条轨道（理想上限）。差距缩小：试点变异组挽回的“默认到理想上限”差距的比例。",
    "s8": r"提升在哪里（事后分析，描述性）：按不用标签的默认预测与理想轨道（需要全部标签）之间的差距，把 Siraj 位点分成四个四分位组；给出 $k=8$ 时试点变异组相对默认的提升，以及比默认差 0.05 以上的位点数。",
    "s9": r"跨位点先验（预注册扩展；85 个检验位点，20 次重复，$\lambda=0.75$ 在 40 个开发位点上选定）。新位点的先验是冻结先验与同一细胞类型中另外 $m$ 个位点（各自测完 8 个试点变异后）的平均后验的混合。X1（C3，$k=8$，跨位点减单位点至少 0.02）：+0.0116，65/85 个位点为正，单侧 Wilcoxon $P=4\times10^{-8}$；未达到。X2（C1，$k=8$，不劣于匹配轨道平均）：+0.0141，下界 +0.0069，3/85 个位点差了 0.05 以上（上界 0.089）；达到。",
    "s10": r"神经元 MPRA（预注册的第三个数据集；WTC11-NGN2 诱导兴奋性神经元，159 个区域，200 次重复）。测了 $k$ 个试点或随机变异后 TrackPilot 的平均 Spearman 相关，以及不用标签的默认预测（匹配轨道、相近细胞或全轨道平均）。AlphaGenome：N1 $-$0.014（未达到），N2 +0.009（$P=0.027$，成立），N3 $-$0.019，下界 $-$0.028，47/159 个区域差了 0.05 以上（未达到），N4 $-$0.015（未达到）；Enformer：N1 $-$0.035（未达到），N2 +0.005（$P=0.17$，未达到），N3 $-$0.030，下界 $-$0.038，48/159 个区域差了 0.05 以上（未达到），N4 $-$0.021（未达到）。判定（N1 和 N3）：两个模型都不支持。用全部标签选出的最佳单条轨道达到 0.32 和 0.32（AlphaGenome、Enformer），这主要是在约 1,000 条轨道中、每个片段约 34 个留出变异上做选择的结果。",
    "s11": r"适用性检查（预注册扩展；检验单位，$k=8$，20 次重复，阈值 $P>0.5$ 在开发单位上选定）。测完试点变异组后，用置换检验（把八个标签打乱 200 次）比较混合模型在真实标签和打乱标签下的对数证据；如果真实标签并没有被解释得更好，就返回不用标签的默认预测。表中给出不加检查和加检查的 TrackPilot、默认预测的平均 Spearman 相关，以及触发回退的运行比例。评估使用逐个选择的变异组。A1（神经元：加检查 $-$ 不加检查 $>0$；加检查不劣于默认，界值 0.02）和 A2（Siraj：保住相对默认的提升；加检查不劣于不加检查，界值 0.01）均成立；具体统计量见代码仓库。次要分析，神经元检验区域使用 Enformer 预测：加检查 $-$ 不加检查 $+0.019$（$P=2\times10^{-4}$），加检查 $-$ 默认 $-0.015$（下界 $-0.020$）。",
    "s12": r"适用性检查的阈值（事后分析，描述性；检验单位，$k=8$，20 次重复）。各阈值下，加检查的预测与不用标签的默认预测之间 Spearman 相关的平均差值（置换 $P$ 超过阈值时退回默认）；“无”：不做检查。冻结的阈值为 0.5；在开发单位上 0.2、0.3 和 0.5 的得分相差不超过 0.0005。在神经元实验中，置换 $P$ 值在经验上接近均匀分布，所以阈值 $p^*$ 大约让比例为 $p^*$ 的运行通过。Enf.：Enformer。",
    "s13": r"在没有自身轨道的人类细胞系中测量的 Kircher 检验区域（事后亚组分析；这些细胞系在构建变异表时、任何运行之前，就被标为未测绘，并指定一个相近生物样本作为代理）。冻结的 TrackPilot，$k=8$，200 次重复。“保留代理”：代理轨道作为匹配轨道；“隐藏代理”：代理及相近轨道均被隐藏。TP：TrackPilot。AlphaGenome，5 个区域，保留代理：TrackPilot $-$ 代理平均 = +0.002（3/5 $>0$）。AlphaGenome，5 个区域，隐藏代理：TrackPilot $-$ 全轨道平均 = +0.124（3/5 $>0$）。AlphaGenome，5 个区域，隐藏代理：试点 $-$ 随机组 = +0.108（5/5 $>0$）。Enformer，5 个区域，保留代理：TrackPilot $-$ 代理平均 = +0.010（4/5 $>0$）。Enformer，5 个区域，隐藏代理：TrackPilot $-$ 全轨道平均 = +0.160（4/5 $>0$）。Enformer，5 个区域，隐藏代理：试点 $-$ 随机组 = +0.101（5/5 $>0$）。",
    "s14": r"预注册检验：在小鼠细胞系中测量的 Kircher 实验，没有任何人类轨道测绘过这些细胞（五个区域；两个 ZRS 实验取平均）。在任何打分之前根据元数据指定了人类代理（Neuro-2a 对应神经母细胞瘤细胞系，MIN6 对应内分泌胰腺前体细胞，NIH/3T3 对应成纤维细胞）。冻结的 TrackPilot，$k=8$，200 次重复。假设：M1，隐藏代理时 TP $-$ 全轨道平均 $\ge 0.03$ 且 $\ge 4/5$ 个区域 $>0$；M2，隐藏代理时试点 $-$ 随机组在 $\ge 4/5$ 个区域 $>0$；M3，保留代理时 TP $-$ 代理平均 $\ge -0.02$ 且没有区域低于 $-0.05$；M1 和 M3 同时成立才算支持。AlphaGenome：M1 不成立，M2 成立，M3 不成立；Enformer：M1 不成立，M2 成立，M3 成立。两个模型都没有达到预注册判据。最佳轨道：用全部标签选出的最佳单条轨道（事后）。按预注册要求的描述性结果（20 次重复），适用性检查退回默认预测的运行比例：AlphaGenome 保留代理时 5\%，隐藏代理时 21\%；Enformer 保留代理时 21\%，隐藏代理时 29\%；加检查后预测的平均准确度不低于不加检查的 TrackPilot。这些描述性结果不改变对未加检查方法的判定。",
}

# Whole table cells (after stripping), translated exactly.
CELLS = {
    # headers
    "Quantity": "统计量", "Reference-relative": "相对参照", "Published centring": "原文中心化",
    "Set": "数据集", "Contrast": "对比", "First": "前者", "Second": "后者", r"Acc.\,1": "准确率 1",
    r"Acc.\,2": "准确率 2", "Diff.": "差值", "Hypothesis": "假设", "Held-out fold": "留出折",
    r"Change [95\% CI]": r"变化 [95\% CI]", "Verdict": "判定", "Data": "数据", "View": "场景", "Mean": "均值",
    "Reading": "判读", "E0 first": "E0 第一步", "E0 lowest": "E0 最低", "E1 first": "E1 第一步",
    "E1 lowest": "E1 最低", "Pilot": "试点", "Random": "随机", "Default": "默认", "Oracle": "理想上限",
    "Gap closed": "差距缩小", "Quartile": "四分位组", "Loci": "位点数", "Gap": "差距", "Gain": "提升",
    "Share closed": "缩小比例", r"Worse $>0.05$": r"差 $>0.05$", "Single locus": "单位点", "Cross-locus": "跨位点",
    "Difference": "差值", r"Loci $>0$": r"位点 $>0$", "AG pilot": "AG 试点", "AG random": "AG 随机",
    "AG default": "AG 默认", r"Enf.\ pilot": r"Enf.\ 试点", r"Enf.\ random": r"Enf.\ 随机",
    r"Enf.\ default": r"Enf.\ 默认", "Unchecked": "不加检查", "Checked": "加检查", "Fallback": "回退比例",
    "Neurons C1": "神经元 C1", "Neurons C3": "神经元 C3", "Neurons C1, Enf.": "神经元 C1，Enf.",
    "Neurons C3, Enf.": "神经元 C3，Enf.", "Model": "模型", "Region": "区域", "Cell line": "细胞系",
    "Proxy mean": "代理平均", "TP, with proxy": "TP，保留代理", "All-track mean": "全轨道平均",
    "TP, proxy hidden": "TP，隐藏代理", "Random panel": "随机组", "Best track": "最佳轨道",
    # values
    "equivalent": "等价", "first better": "前者更好", "second better": "后者更好", "undetermined": "未定",
    "holds": "成立", "none": "无", "C1 profiled": "C1 已测绘", "C2 related only": "C2 仅有相近细胞",
    "C3 unprofiled": "C3 未测绘", "Kircher (11 regions)": "Kircher（11 个区域）",
    "Neurons (109 regions)": "神经元（109 个区域）", "Siraj (85 loci)": "Siraj（85 个位点）",
    "Kircher (v1)": "Kircher（v1）", "Kircher (v2)": "Kircher（v2）",
    "Kircher (v2, no offset)": "Kircher（v2，去偏移）",
    "mean D (lower bound)": "平均 D（下界）", "harmed (upper bound)": "受损位点（上界）",
    r"v2 $-$ all-track mean ($P$)": r"v2 $-$ 全轨道平均（$P$）",
    r"v2 $-$ related mean ($P$)": r"v2 $-$ 相近细胞平均（$P$）",
    r"v2 disc panel $-$ v2 random panel": r"v2 试点组 $-$ v2 随机组",
    r"v2 $-$ meta\_ridge on the same disc panel": r"同一试点组上 v2 $-$ meta\_ridge",
    r"matched\_cal (random panel) $-$ matched": r"matched\_cal（随机组）$-$ 匹配轨道",
    r"matched\_cal (disc panel) $-$ matched": r"matched\_cal（试点组）$-$ 匹配轨道",
    r"family prior $-$ uniform prior": r"族先验 $-$ 均匀先验",
    r"v2 $-$ matched": r"v2 $-$ 匹配轨道", r"v2 $-$ all\_mean": r"v2 $-$ 全轨道平均",
    "TERT (glioblastoma)": "TERT（胶质母细胞瘤）", "glioblastoma (alt.)": "胶质母细胞瘤（替代）",
}

# Substring rules applied to every table cell after the exact match.
SUBS = [
    (" to ", " 至 "), (" harmed", " 受损"), ("(single track)", "（单条轨道）"), ("(average)", "（平均）"),
    (r"coda: average $-$ single track", r"coda：平均 $-$ 单条轨道"), (" vs random", " 对随机"),
]

NUM = re.compile(r"\d+(?:\.\d+)?")


def numbers(text):
    return Counter(NUM.findall(text))


def translate_row(line):
    body, end = line.rstrip("\n"), ""
    m = re.search(r"\s*\\\\\s*$", body)
    if m:
        body, end = body[: m.start()], " \\\\"
    cells = []
    for cell in body.split("&"):
        s = cell.strip()
        t = CELLS.get(s, s)
        for a, b in SUBS:
            t = t.replace(a, b)
        cells.append(" " + t + " " if cell.startswith(" ") else t + " ")
    return "&".join(cells).rstrip() + end + "\n"


def main():
    src = open(SRC, encoding="utf-8").read()
    out, problems = [], []
    for line in src.splitlines(keepends=True):
        if r"\caption{" in line:
            lab = re.search(r"\\label\{tab:(\w+)\}", line).group(1)
            en_cap = re.search(r"\\caption\{(.*)\}\\label", line).group(1)
            zh_cap = CAPTIONS[lab]
            if numbers(en_cap) != numbers(zh_cap):
                problems.append(f"{lab}: caption numbers differ: EN-only {numbers(en_cap) - numbers(zh_cap)}, "
                                f"ZH-only {numbers(zh_cap) - numbers(en_cap)}")
            line = line[: line.index(r"\caption{")] + r"\caption{" + zh_cap + "}" + line[line.index(r"\label{"):]
            out.append(line)
        elif "&" in line and not line.lstrip().startswith("%"):
            out.append(translate_row(line))
        else:
            out.append(line)
    zh = "".join(out)
    if numbers(src) != numbers(zh):
        problems.append(f"table numbers differ: EN-only {numbers(src) - numbers(zh)}, ZH-only {numbers(zh) - numbers(src)}")
    leftovers = sorted({w for w in re.findall(r"\b(?:better|equivalent|undetermined|profiled|holds|harmed)\b", zh)})
    if leftovers:
        problems.append(f"untranslated words: {leftovers}")
    if problems:
        sys.exit("\n".join(problems))
    open(DST, "w", encoding="utf-8").write(
        "% Generated by paper/make_supplement_zh.py from supp_tables.tex; do not edit by hand.\n" + zh)
    print("wrote", DST)


if __name__ == "__main__":
    main()
