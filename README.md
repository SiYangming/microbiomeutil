# microbiomeutil（r20110519）

本仓库归档 Broad Institute（Brian Haas）出品的 **microbiomeutil** 工具集（版本 `r20110519`，2011-05-19）的代码、文档与轻量数据；**完整发布包（127 MB，含大型参考数据库）见本仓库 Release 附件**。

## 归档说明

- 发布包 `microbiomeutil-r20110519.tgz`：133,623,556 B（127 MB），MD5 `11eaac4b0468c05297ba88ec27bd4b56`，SHA-256 `9233de80ea57bfb9e9371cbe7e3bfad2d4a51168fddaf60fa144c4046c80d823`；解压后为 `microbiomeutil-r20110519/`，共 **152** 个文件
- **入库范围**：发布包中小于 5 MB 的 **146** 个文件（代码、脚本、文档、小数据）
- **未入库、仅在 Release 附件中的 6 个大文件**（均为参考数据库/大比对文件，约 234 MB）：

| 文件 | 大小 |
|---|---|
| `WigeoN/data/WigeoN.AllVsAllReferences.txt.gz` | 105.22 MB |
| `RESOURCES/RDP/classifier_data_dir/genus_wordConditionalProbList.txt` | 38.89 MB |
| `RESOURCES/rRNA16S.gold.NAST_ALIGNED.fasta` | 38.66 MB |
| `ChimeraSlayer/__ChimeraConstructorToolkit/SeqSubstrates/rRNA16S.toChimerize.fasta.NAST` | 35.58 MB |
| `RESOURCES/rRNA16S.gold.fasta` | 8.33 MB |
| `ChimeraSlayer/__ChimeraConstructorToolkit/SeqSubstrates/rRNA16S.toChimerize.fasta` | 7.65 MB |

> 其中 `WigeoN.AllVsAllReferences.txt.gz` 达 105.22 MB，**超过 GitHub 单文件 100 MB 硬上限**，无法入库；其余同属大型参考数据，一并只在附件中提供。
>
> **要完整运行**：下载 Release 附件解压，再把 `WigeoN/`、`RESOURCES/`、`ChimeraSlayer/` 与本仓库对应目录合并（或直接用附件中的完整目录）。

- 上游：SourceForge 项目 `https://sourceforge.net/projects/microbiomeutil/`，文件页 `.../files/microbiomeutil-r20110519.tgz/download`
  - 2026-09 实测仍可下载：`Content-Length` 133,623,556 与本地一致、`last-modified` 2011-05-19、前 2 MB MD5 一致
  - 但实测下载速度仅约 **50 KB/s**，127 MB 需 40 分钟以上，这也是本地归档的原因
- 其他渠道：**bioconda 无 `microbiomeutil` recipe**（2026-09 核实 404）；GitHub 无该项目镜像
- 许可：各子目录自带 LICENSE —— `ChimeraSlayer/LICENSE`、`NAST-iEr/LICENSE`、`WigeoN/LICENSE`、`TreeChopper/LICENSE`、`AmosCmp16Spipeline/LICENSE`
- 上游原始 `README`（根目录 `README`）与各子目录 README 均按原样保留，未做改动

## 包内构成（据各子目录自带 README）

| 目录 | 说明 |
|---|---|
| `ChimeraSlayer/` | 嵌合体检测（`ChimeraSlayer.pl` + `PerlLib/`、`util/`、`sample_data/`，另有 `__ChimeraConstructorToolkit`、`__KmerChimerDetector`、`__BroadBellerophon`、`ChimeraParentSelector`、`ChimeraPhyloChecker`） |
| `NAST-iEr/` | 将单条原始核酸序列比对到一个或多个 NAST 格式序列（固定模板、无末端罚分的全局动态规划）；C 程序，`gcc NAST-iEr.c -o NAST-iEr`（或 `make`）编译 |
| `WigeoN/` | Pintail 算法的再实现（不惩罚查询中的 `N`）：比较查询与可信参考序列在 NAST 比对下的保守性，观测变异超过非异常序列变异分布 95% 分位即标记为异常；用包装脚本 `run_WigeoN.pl` 运行 |
| `TreeChopper/` | 按系统发育距离聚类树叶节点：叶间距离阈值建图 → 计算 Jaccard 相似系数剪边 → 传递闭包聚类输出 OTU |
| `AmosCmp16Spipeline/` | 用 AMOScmp 做 16S reads 的参考辅助组装：选最佳参考 → Lucy 修剪低质量末端 → 同源修剪 → AMOScmp 组装 → 按参考锚定用 `N` 补洞并给质量值 |
| `RESOURCES/` | 参考 16S 数据：`rRNA16S.gold.fasta`（近全长、相对无嵌合，取自 RDP/细菌基因组）及其 NAST 比对 `rRNA16S.gold.NAST_ALIGNED.fasta`（致谢 RDP 的 Jim Cole） |
| `docs/` | `microbiomeutil.asciidoc`、`microbiomeutil.html`（工具用法文档）、`old.microbiomeutil.asciidoc`、`images/` |

## 依赖（上游 README 与各子目录 README）

需另行安装并加入 `PATH`：

- `megablast`（NCBI BLAST）—— NAST-iEr、WigeoN
- `cdbtools`（cdbfasta/cdbfasta）—— NAST-iEr、WigeoN
- `Bioperl`、`slclust`（后者需放入 `TreeChopper/util/`）—— TreeChopper
- `Mummer`、`Amos` —— AmosCmp16Spipeline；需用包内 `AmosPatch/AMOScmp` 替换 Amos 官方的 `amos/bin/AMOScmp` 脚本

## 构建与测试

```bash
# 在包根目录
make            # 仅编译 C 写的 NAST-iEr；其余工具为 Perl，无需编译
make test       # 依次执行 testNast、testWigeon、testChimeraSlayer
# 也可单独执行：make testNast / make testWigeon / make testChimeraSlayer
```

顶层 `Makefile` 的 `DIRS = NAST-iEr ChimeraSlayer WigeoN RESOURCES`，`make all/clean` 会依次进入这些子目录执行。

## 使用

- 各工具的参数与用法：`docs/microbiomeutil.html`（源文件 `docs/microbiomeutil.asciidoc`）及各子目录的 README
- 各工具 `sample_data/` 目录内有示例输入
- 上游网站：`http://microbiomeutil.sf.net`（联系人：Brian Haas，`bhaas@broadinstitute.org`）
