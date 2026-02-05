# draw_fusions

可视化基因融合数据的R脚本工具，支持生成融合基因结构图、Circos图、蛋白质结构域图和reads覆盖图。

上游来源：[arriba](https://github.com/suhrig/arriba) - GPLv3

## 代码重构

本项目已进行代码重构，主要改进：

### 配置系统
- 统一的 `DEFAULT_CONFIG` 配置中心
- 分类管理：输出设置、颜色、尺寸、布局、阈值
- 支持命令行参数覆盖默认值

### 函数重构
- `findExons()` 使用查找表替代12层嵌套ifelse
- 转录本优先级映射表：APPRIS标注 > CCDS > RefSeq > Ensembl
- 辅助函数模块化：10+个独立函数

### 命名规范
- 变量：snake_case（如 `exon_height`、`pdf_width`）
- 常量：UPPER_SNAKE_CASE（如 `MAX_EXON_LABEL_LENGTH`）
- 函数：动词+名词（如 `calculate_coverage_region`）

## 许可证

本项目源自 [arriba](https://github.com/suhrig/arriba) 项目，采用 GPLv3 开源协议发布。

## 输入模式

本工具支持两种输入模式：

### 模式1：文件输入模式（兼容原格式）

从 fusion.tsv 文件读取融合数据，支持以下工具的输出格式：
- **Arriba**
- **STAR-Fusion**

### 模式2：单融合参数模式（新增功能）

通过命令行参数直接指定单个融合基因信息，无需准备文件：

```
--gene1=GENE1 --gene2=GENE2 --contig1=CHR1 --contig2=CHR2 --breakpoint1=POS1 --breakpoint2=POS2
```

此模式保留了原文件格式的容错处理，当未指定 `--fusions` 参数时自动启用。

## 依赖

- R (>= 3.6)
- R包：
  - `GenomicRanges` - 用于蛋白质结构域和Circos图
  - `circlize` - 用于Circos图
  - `GenomicAlignments` - 用于reads覆盖图（可选）

## 安装

### 使用Docker(推荐)

```bash
# 拉取预构建镜像
docker pull ghcr.io/schemabio/draw_fusion:latest

# 或指定版本
docker pull ghcr.io/schemabio/draw_fusion:v0.0.3

# 运行
docker run --rm -v $(pwd):/data ghcr.io/schemabio/draw_fusion --help
```

### 本地构建

```bash
# 构建镜像
docker build -t draw_fusions .

# 运行
docker run --rm -v $(pwd):/data draw_fusions --help
```

## 发布新版本

通过提交包含版本号的 commit 来触发自动构建：

```bash
# 提交时在 message 中包含版本号，如：
git commit -m "v0.0.3"

# 推送到 main 分支
git push origin main
```

GitHub Actions 会自动从 commit message 中提取版本号，构建并推送 Docker 镜像。

## 使用方法

**注意**：两种输入模式 **互斥**，二选一

### 模式1：文件输入（批量融合）

```bash
Rscript draw_fusions.R \
  --fusions=fusions.tsv \          # 必填：融合基因文件
  --annotation=annotation.gtf \    # 必填：GTF注释文件
  --output=output.pdf             # 必填：输出PDF
```

### 模式2：单融合参数（单个融合）

```bash
Rscript draw_fusions.R \
  --annotation=annotation.gtf \    # 必填：GTF注释文件
  --output=fusion.pdf \           # 必填：输出PDF
  # 单融合模式必填参数：
  --gene1=TMPRSS2 \               # 基因1名称
  --gene2=ERG \                   # 基因2名称
  --contig1=chr2 \                # 染色体1 (chr2 或 2)
  --contig2=chr2 \                # 染色体2
  --breakpoint1=42526250 \        # 断点1位置
  --breakpoint2=29447242          # 断点2位置

# 可选参数（不提供则自动推断）：
#   --transcript_id1=ENST00000439425   # 基因1的转录本ID
#   --transcript_id2=ENST00000355722   # 基因2的转录本ID
#   # direction 参数已废弃，会从 GTF 的 strand 自动推断

# fusion_type 和 direction 会根据参数自动推断：
#   - fusion_type: 根据 contig1/contig2 自动判断
#   - direction: 根据 GTF 中基因的 strand 自动推断 (+→downstream, -→upstream)
```

### 完整参数列表

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| **输入模式（二选一）** |
| `--fusions` | 文件 | 否* | - | 融合基因文件（Arriba/STAR-Fusion格式） |
| `--gene1` | 字符串 | ✅ | "" | 基因1名称 |
| `--gene2` | 字符串 | ✅ | "" | 基因2名称 |
| `--contig1` | 字符串 | ✅ | "" | 染色体1编号 |
| `--contig2` | 字符串 | ✅ | "" | 染色体2编号 |
| `--breakpoint1` | 数值 | ✅ | 0 | 断点1位置 |
| `--breakpoint2` | 数值 | ✅ | 0 | 断点2位置 |
| **必需参数** |
| `--annotation` | 文件 | ✅ | annotation.gtf | GTF注释文件 |
| `--output` | 字符串 | ✅ | output.pdf | 输出PDF文件名 |
| **可选参数** |
| `--strand1` | 字符串 | 否 | "." | 链信息 |
| `--strand2` | 字符串 | 否 | "." | 链信息 |
| `--transcript_id1` | 字符串 | 否 | "." | 基因1的转录本ID |
| `--transcript_id2` | 字符串 | 否 | "." | 基因2的转录本ID |
| `--transcriptSelection` | 字符串 | 否 | provided | 转录本选择策略 |
| `--site1` | 字符串 | 否 | exon | 断点位置类型 |
| `--site2` | 字符串 | 否 | exon | 断点位置类型 |
| `--reading_frame` | 字符串 | 否 | "." | 读码框 |
| `--split_reads` | 数值 | 否 | 0 | split reads数 |
| `--discordant_mates` | 数值 | 否 | 0 | discordant mates数 |
| `--confidence` | 字符串 | 否 | high | 可信度 |
| `--alignments` | 文件 | 否 | - | BAM比对文件 |
| `--cytobands` | 文件 | 否 | cytobands.tsv | 细胞带文件 |
| `--proteinDomains` | 文件 | 否 | protein_domains.gff3 | 蛋白质结构域 |
| `--sampleName` | 字符串 | 否 | "" | 样本名称 |
| `--plotPanels` | 字符串 | 否 | fusion,circos,domains,readcounts | 面板组合 |
| `--color1` | 字符串 | 否 | #e5a5a5 | 基因1颜色 |
| `--color2` | 字符串 | 否 | #a7c4e5 | 基因2颜色 |
| `--fontSize` | 数值 | 否 | 1 | 字体大小 |
| `--fontFamily` | 字符串 | 否 | Helvetica | 字体家族 |
| `--squishIntrons` | 布尔 | 否 | TRUE | 压缩内含子 |
| `--render3dEffect` | 布尔 | 否 | TRUE | 3D效果 |

> * 单融合模式需要 `--gene1`, `--gene2`, `--contig1`, `--contig2`, `--breakpoint1`, `--breakpoint2` 同时提供

## 配置系统

### 可配置参数类别

#### 输出设置
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `output$pdf_width` | 11.692 | PDF宽度（A4横向） |
| `output$pdf_height` | 8.267 | PDF高度（A4横向） |
| `output$font_family` | Helvetica | 字体家族 |

#### 颜色配置
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `colors$gene1` | #e5a5a5 | 基因1颜色 |
| `colors$gene2` | #a7c4e5 | 基因2颜色 |
| `colors$exon_fill` | #ffffff | 外显子填充色 |
| `colors$arc_fill` | #b6cbe3 | 弧线填充色 |

#### 尺寸配置
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `dimensions$arc_steps` | 30 | 弧线分段数 |
| `dimensions$exon_height` | 0.03 | 外显子高度 |
| `dimensions$ideogram_height` | 0.02 | 染色体图高度 |
| `dimensions$domain_height` | 0.25 | 结构域高度 |

#### 布局配置
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `layout$y_sample_name` | 1.04 | 样本名称Y位置 |
| `layout$y_coverage` | 0.72 | 覆盖图Y位置 |
| `layout$y_fusion` | 0.35 | 融合图Y位置 |
| `layout$y_ideogram` | 0.75 | 染色体图Y位置 |

#### 阈值配置
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `thresholds$split_reads` | 2 | 最小split reads |
| `thresholds$discordant_mates` | 10 | 最小discordant mates |
| `thresholds$max_domains` | 50 | 最大显示结构域数 |
| `thresholds$coverage_zoom` | 1e6 | 覆盖度缩放阈值 |

## 输出

生成包含以下面板的PDF文件：
- **fusion**: 融合基因结构图
- **circos**: Circos基因组环图
- **domains**: 蛋白质结构域图
- **readcounts**: reads覆盖图

