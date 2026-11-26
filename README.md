# Tracepick · Sanger 测序引物设计

[打开在线工具](https://rei402.github.io/Tracepick-Sanger/)

把双链 DNA 原序列或一条 FASTA 记录粘贴到网页，点击 **设计测序引物**。工具为每次独立的 Sanger 反应输出一条 5′ → 3′ 引物，每次设计最多给出五条。所有计算均在浏览器本地完成。

## 使用与部署

1. 本地使用时解压完整发布包，用 Chrome、Edge、Firefox 等支持 WebAssembly 的浏览器打开 `docs/index.html`。所有页面资源均随包提供。
2. 发布到 GitHub Pages 时，将解压后的项目内容放在仓库根目录，保持 `docs/index.html` 的层级。进入 **Settings → Pages → Deploy from a branch**，选择保存文件的分支和 **/docs** 后保存。网页无需后端、密钥或构建服务。[GitHub 官方部署指南](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site)
3. 支持单条 1–20,000 bp 的 A/C/G/T 模板。空白、换行和大小写不影响坐标。遇到模糊碱基时，界面指出其位置。完整质粒选择“环状 DNA”，已切开的片段或 PCR 产物选择“线性 DNA”。若只输入插入片段，模板唯一性检查的范围也只有该片段。

网页打开时默认选择“读完指定区间或位点”、线性 DNA、正向和反向，目标起止均为 1 bp。目标坐标和预计读段范围始终显示。

## 选择哪种模式

| 目的 | 设置 | 怎样使用五条结果 |
| --- | --- | --- |
| 想知道整段模板最多能读到哪些位置 | 尽量覆盖整段模板 | 依次加入贡献最多新碱基的读段。看合并覆盖和缺口，决定需要做几次反应。 |
| 想验证某个碱基，例如第 720 位 | 读完指定区间或位点，起止均填 720 | 每条候选都可独立读到该位点。两个方向都有合格候选时，前两条分别来自两侧。 |
| 想用一条读段读过整个接头 | 输入接头两侧的目标起止坐标 | 每条候选的预计读段均包含完整区间。目标超过一次读段长度或缺少可用侧翼时，可能返回零条。 |

所有显示坐标均从 1 开始、两端均计入，正反方向都相对于输入序列。跨越环状原点时会显示两段坐标。反向引物已做反向互补，点击序列即可复制；导出 TSV 包含方向、结合区、读段区、Tm、GC、两种结构 Tm、目标和运行设置。

## 设计依据

| 原始资料 | 对设计的影响 |
| --- | --- |
| Cornell 推荐约 17–25 nt、优选 Tm 55–65 °C，强调唯一结合位点和 3′ 端精确配对。UBC 推荐 18–25 nt、Tm 55–60 °C。 | 本工具取 18–25 nt、55–65 °C，排序偏好 60 °C。不同机构使用的化学条件和 Tm 算法会影响具体数值。[Cornell](https://www.biotech.cornell.edu/core-facilities-brc/facilities/genomics-facility/qc-troubleshooting/sanger-sequencing-primer-design) · [UBC](https://sequencing.ubc.ca/our-services-equipment/sanger-sequencing/sanger-faq) |
| 芝加哥大学给出常见 GC 比例 40–60%。Monash 提醒检查连续单碱基、发卡、二聚体和额外结合位点。 | 使用 GC 40–60%，排除连续四个相同碱基及四组相同的串联二核苷酸。3′ 末位 G/C 仅作排序偏好。[芝加哥大学](https://voices.uchicago.edu/dnaseq/getting-started/order/) · [Monash 原始指南](https://www.monash.edu/__data/assets/pdf_file/0008/2014964/primer-design.pdf) |
| MD Anderson 建议引物距目标至少约 50 bp；Cornell 提醒长读段后部的可靠性可能下降。 | 网页按用户设定，默认把距引物 3′ 端第 30–750 枚新碱基计为预计读段（最多 721 bp）。质粒只计当前一圈，线性模板在端点截断。可按服务商实测表现修改读段距离。[MD Anderson](https://www.mdanderson.org/research/research-resources/core-facilities/advanced-technology-genomics-core/services-and-fees/sanger-sequencing.html) · [Cornell](https://www.biotech.cornell.edu/core-facilities-brc/facilities/genomics-facility/qc-troubleshooting/sanger-sequencing-primer-design) |

### 开源工具比较

| 工具 | 原有能力 | 本项目如何使用 |
| --- | --- | --- |
| Primer3，GPL-2.0-or-later | 内置 `pick_sequencing_primers`，支持位置、间隔和热力学约束。 | **复用代码**。将 Primer3 THAL 编译为随页面发布的 WebAssembly，实际计算候选发卡和同源二聚体 Tm。位置搜索与五条候选排序由本项目实现。[手册](https://primer3.org/manual.html) · [源码](https://github.com/primer3-org/primer3) |
| Primavera，MIT | 可从已有引物库中为装配体的重点区域选择验证引物，并绘制预计覆盖。 | **借鉴方法**。从输入序列枚举新引物，用读段并集、双向覆盖和缺口说明选择结果。[项目及示例](https://github.com/Edinburgh-Genome-Foundry/Primavera) |
| sequencing-primer-generator，MIT | 面向整质粒的批量测序引物设计，可导出订购表格。 | **借鉴交付方式**。候选以订购方向展示，一键导出 TSV；同时保留结合位置和预计读段。[项目说明](https://github.com/irahorecka/sequencing-primer-generator) |

### 筛选和排序的精确定义

- 在允许的方向上枚举全部 18–25 nt 窗口。寡核苷酸 Tm 用 SantaLucia 1998 最近邻参数及熵盐校正计算，固定 Na⁺ 50 mM、Mg²⁺ 0、dNTP 0、总链浓度 50 nM。非自互补序列浓度项为 Ct/4，自互补序列采用对称熵修正和 Ct。这是可复算的计算条件，服务商实际反应体系可能不同。[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC19045/) · [Biopython 方法文档](https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html)
- 在模板两条链上检查每一处额外位点。末端 12 nt 完全匹配，或全长无缺口比对仅有 0–1 处替换，任一成立即排除。前后两个互不重叠的 8-mer 中至少有一个必然命中全长至多一处替换的位点。环状模板的索引包含跨原点位点。这些是本工具选择的保守规则，可能排除实验上可用的引物。混合样本、模板外同源序列、插缺和更多错配超出搜索范围。
- 初筛排除 3′ 端连续 4 bp 自互补、任意连续 6 bp 自互补，或具有至少 4 bp 互补茎和 3 nt 环的发卡。按排名调用 THAL（37 °C 参考温度和上述离子条件），发卡和同源二聚体的预测 Tm 均须低于 47 °C。遇到不合格者继续复核下一条。不同引物用于独立反应。[Primer3 热力学参数](https://primer3.org/manual.html)
- 整段模式采用贪心选择。优先补足全模板未覆盖的碱基数，其次增加该方向尚未覆盖的碱基数，再按引物本身的性质择优。当五条仍留有空白时，按实际结果显示缺口。该算法没有全局最优保证，也不强制读段间重叠。
- 目标模式首先要求单条读段覆盖完整区间，再优先提供两个方向。两种模式均使用相同的引物性质代价：与 60 °C 的温差，加上 0.08 × GC 百分点偏差（参考 50%）、0.18 × 长度偏差（参考 21 nt）及 0.4 的非 GC 末端代价。同方向已选 3′ 端至少相隔 12 bp。权重、距离和 5 条上限均为产品选择，代价不表示实验成功概率。

## 复核

用 Node.js 运行 `node --test tests/engine.test.cjs`。安装 Biopython 和 primer3-py 的环境中可再运行 `python tests/audit.py`。网页日常使用无需安装这些依赖。测试模板为 ENA 登录号 L09137.2 的 pUC19c（2,686 bp），FASTA 原文保存在 tests 下，去除换行和标题后的序列 SHA-256 为 `8a205af08518dae1b0045c7b098986dca18824ef654e75233603602ab89cad04`。[ENA 原始记录](https://www.ebi.ac.uk/ena/browser/api/fasta/L09137.2)

软件测试用独立逐碱基枚举复核模板搜索、方向、环状坐标、读段并集，包含重复模板、短模板和 0–4 条结果。数值复核比较 JavaScript 与 Biopython、Primer3 原生实现，并比较浏览器 THAL 与原生 THAL 的结构 Tm。验证结果见 [VALIDATION.md](VALIDATION.md)。

实验效能仍待测量。可在有已知参考序列和历史成功引物的质粒上，将候选与历史引物分别提交至同一测序平台，保留全部候选及失败记录。以原始 AB1 峰图独立检查读段方向、目标区 Q20 碱基比例及与参考的一致性，按目标完整读通率评价；用这批结果检验默认窗口后，应在另一批模板上验证调整后的参数。

## 许可与源码

项目按 GPL-2.0-or-later 发布。原始代码、THAL 接口、附带的上游代码、第三方声明和重建命令位于 `native/`。网页的热力学源码下载包可独立保存，详见 [native/README.md](native/README.md)。
