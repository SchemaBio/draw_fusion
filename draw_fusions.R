#!/usr/bin/env Rscript

# =============================================================================
# draw_fusions.R - 基因融合可视化工具
# =============================================================================
# 基于 https://github.com/suhrig/arriba 重构
# 许可证: GPLv3
#
# 重构特性:
#   - 配置系统 - 所有参数可配置
#   - 查找表 - 替换深层嵌套ifelse
#   - 模块化 - 拆分复杂函数
#   - 命令行参数 - 支持单融合模式
# =============================================================================
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================

# -----------------------------------------------------------------------------
# 警告设置
# -----------------------------------------------------------------------------
options(warn = 1)


# =============================================================================
# 第一部分：配置系统
# =============================================================================

DEFAULT_CONFIG <- list(
  # 输出设置
  output = list(
    pdf_width = 11.692,
    pdf_height = 8.267,
    font_family = "Helvetica",
    font_size = 1.0
  ),

  # 颜色设置
  colors = list(
    gene1 = "#e5a5a5",
    gene2 = "#a7c4e5",
    background = "#eeeeee",
    ideogram_gneg = "#ffffff",
    ideogram_acen = "#ec4f4f",
    ideogram_stalk = "#0000ff",
    fusion_types = c(
      translocation = "#000000",
      duplication = "#00bb00",
      deletion = "#ff0000",
      inversion = "#0000ff"
    )
  ),

  # 尺寸参数
  dimensions = list(
    arc_steps = 30,
    curly_brace_height = 0.03,
    ideogram_height = 0.04,
    ideogram_width = 0.4,
    exon_height = 0.03,
    gradient_steps_exon = 10,
    gradient_steps_domain = 20,
    max_resolution = 5000,
    squashed_intron_size = 200,
    strand_arrow_offset = 0.001,
    strand_arrow_spacing = 0.005,
    strand_arrow_length_scaled = 0.01,
    exon_border_padding = 0.001,
    domain_padding = 0.025,
    ideogram_arc_steps = 30
  ),

  # 布局坐标
  layout = list(
    y_sample_name = 1.04,
    y_ideograms = 0.94,
    y_breakpoint_labels = 0.86,
    y_coverage = 0.72,
    y_exons = 0.67,
    y_gene_names = 0.58,
    y_fusion = 0.50,
    y_transcript = 0.45,
    y_scale = 0.407,
    label_spacing = 0.05,
    desired_scale_size = 0.20,
    scale_whisker_size = 0.007
  ),

  # 绘图设置
  drawing = list(
    smoothness = 20,
    draw_3d_effect = TRUE,
    strand_arrow_length = 0.01,
    strand_arrow_spacing = 0.01
  ),

  # Circos设置
  circos = list(
    min_confidence = "medium",
    confidence_rank = c(low = 0, medium = 1, high = 2),
    arc_line_width = 2,
    legend_columns = 2
  ),

  # 蛋白质域设置
  domain = list(
    merge_overlap_threshold = 0.9,
    optimize_colors = FALSE,
    exon_height = 0.2,
    gene_names_y_offset = 0.05
  )
)

# 使用默认配置
config <- DEFAULT_CONFIG


# =============================================================================
# 第二部分：工具函数 - 颜色处理
# =============================================================================

changeColorBrightness <- function(color, delta) {
  rgb_vals <- col2rgb(color)
  rgb(
    pmin(255, pmax(0, rgb_vals["red", ] + delta)),
    pmin(255, pmax(0, rgb_vals["green", ] + delta)),
    pmin(255, pmax(0, rgb_vals["blue", ] + delta)),
    maxColorValue = 255
  )
}

getDarkColor <- function(color) {
  changeColorBrightness(color, -100)
}

getBrightColor <- function(color) {
  changeColorBrightness(color, +190)
}


# =============================================================================
# 第三部分：工具函数 - 染色体处理
# =============================================================================

addChr <- function(contig) {
  ifelse(contig == "MT", "chrM", paste0("chr", contig))
}

# 染色体名称标准化：统一格式，支持 chr1<->1, chrY<->Y, chrM<->MT 双向映射
# 返回不带 "chr" 前缀的格式（如 "1", "Y", "MT"）
normalizeChromosome <- function(contig) {
  # 首先处理 chrM -> MT 的特殊情况
  contig <- sub("^chrM$", "MT", contig, perl = TRUE)
  # 然后去除其他 chr 前缀
  sub("^chr", "", contig, perl = TRUE)
}

# 兼容旧函数名
removeChr <- function(contig) {
  normalizeChromosome(contig)
}

between <- function(value, start, end) {
  value >= start & value <= end
}


# =============================================================================
# 第三部分补充：工具函数 - 内联压缩逻辑
# =============================================================================

squishIntronsLogic <- function(exons_data, breakpoint, squashed_size) {
  cumulative_intron_length <- 0
  previous_exon_end <- -squashed_size

  for (exon_idx in seq_len(nrow(exons_data))) {
    if (breakpoint > previous_exon_end + 1 && breakpoint < exons_data[exon_idx, "left"]) {
      breakpoint <- (breakpoint - previous_exon_end) /
                   (exons_data[exon_idx, "left"] - previous_exon_end) *
                   squashed_size + previous_exon_end - cumulative_intron_length
    }
    if (exons_data[exon_idx, "left"] > previous_exon_end) {
      cumulative_intron_length <- cumulative_intron_length +
                                  exons_data[exon_idx, "left"] - previous_exon_end -
                                  squashed_size
      previous_exon_end <- exons_data[exon_idx, "right"]
    }
    if (breakpoint >= exons_data[exon_idx, "left"] &&
        breakpoint <= exons_data[exon_idx, "right"] + 1) {
      breakpoint <- breakpoint - cumulative_intron_length
    }
    exons_data[exon_idx, "left"] <- exons_data[exon_idx, "left"] - cumulative_intron_length
    exons_data[exon_idx, "right"] <- exons_data[exon_idx, "right"] - cumulative_intron_length
  }

  list(exons = exons_data, breakpoint = breakpoint)
}


# =============================================================================
# 第四部分：工具函数 - 转录本优先级查找表
# =============================================================================

get_transcript_priority_map <- function() {
  c(
    "appris_principal_1" = 12,
    "appris_principal_2" = 11,
    "appris_principal_3" = 10,
    "appris_principal_4" = 9,
    "appris_principal_5" = 8,
    "appris_principal" = 7,
    "appris_candidate_longest" = 6,
    "appris_candidate" = 5,
    "appris_alternative_1" = 4,
    "appris_alternative_2" = 3,
    "appris_alternative" = 2,
    "CCDS" = 1
  )
}

calculate_transcript_priority <- function(attributes) {
  priority_map <- get_transcript_priority_map()
  sapply(attributes, function(attr) {
    for (pattern in names(priority_map)) {
      if (grepl(pattern, attr)) {
        return(priority_map[[pattern]])
      }
    }
    return(0)
  })
}


# =============================================================================
# 第五部分：findExons() 重构 - 使用查找表和辅助函数
# =============================================================================

findExons <- function(exons, contig, geneID, direction, breakpoint,
                      coverage = NULL, transcriptId = ".",
                      transcriptSelection = "provided") {

  # 策略1：使用指定的转录本
  if (transcriptSelection == "provided" && transcriptId != "." && transcriptId != "") {
    candidate_exons <- exons[exons$transcript == transcriptId, ]
    if (nrow(candidate_exons) > 0) {
      return(candidate_exons)
    }
    warning(sprintf("Unknown transcript '%s', selecting alternative", transcriptId))
  }

  # 策略2：根据选择模式获取候选转录本
  candidate_exons <- switch(
    transcriptSelection,
    "canonical" = get_canonical_transcript_exons(exons, contig, geneID),
    "coverage" = get_coverage_based_transcript(
      exons, contig, geneID, direction, breakpoint, coverage
    ),
    get_fallback_transcript_exons(exons, contig, geneID, direction, breakpoint)
  )

  # 策略3：使用查找表选择共识转录本
  if (length(unique(candidate_exons$transcript)) > 1) {
    priorities <- calculate_transcript_priority(candidate_exons$attributes)
    max_priority <- max(priorities)
    candidate_exons <- candidate_exons[priorities == max_priority, ]
  }

  # 策略4：按编码序列长度选择
  if (length(unique(candidate_exons$transcript)) > 1) {
    candidate_exons <- select_by_coding_length(candidate_exons)
  }

  # 策略5：按总长度选择
  if (length(unique(candidate_exons$transcript)) > 1) {
    candidate_exons <- select_by_total_length(candidate_exons)
  }

  # 最终：返回第一个转录本
  first_transcript <- unique(candidate_exons$transcript)[1]
  return(candidate_exons[candidate_exons$transcript == first_transcript, ])
}


# --- findExons 辅助函数 ---

get_canonical_transcript_exons <- function(exons, contig, geneID) {
  exons[exons$geneID == geneID & exons$contig == contig, ]
}

get_fallback_transcript_exons <- function(exons, contig, geneID, direction, breakpoint) {
  splice_site_transcripts <- exons[
    exons$geneID == geneID &
    exons$contig == contig &
    exons$type == "exon" &
    ((direction == "downstream" & abs(exons$end - breakpoint) <= 2) |
     (direction == "upstream" & abs(exons$start - breakpoint) <= 2)),
    "transcript"
  ]

  if (length(splice_site_transcripts) > 0) {
    return(exons[exons$transcript %in% splice_site_transcripts, ])
  }
  exons[exons$geneID == geneID & exons$contig == contig, ]
}

get_coverage_based_transcript <- function(exons, contig, geneID,
                                         direction, breakpoint, coverage) {

  splice_site_transcripts <- exons[
    exons$geneID == geneID &
    exons$contig == contig &
    exons$type == "exon" &
    ((direction == "downstream" & abs(exons$end - breakpoint) <= 2) |
     (direction == "upstream" & abs(exons$start - breakpoint) <= 2)),
    "transcript"
  ]

  if (length(splice_site_transcripts) > 0) {
    return(exons[exons$transcript %in% splice_site_transcripts, ])
  }

  candidate_exons <- exons[exons$geneID == geneID & exons$contig == contig, ]

  if (!is.null(coverage) && nrow(candidate_exons) > 0) {
    coverage_by_transcript <- calculate_transcript_coverage(candidate_exons, coverage)
    length_by_transcript <- calculate_transcript_lengths(candidate_exons)
    best_transcript <- select_best_transcript(coverage_by_transcript, length_by_transcript)

    if (!is.null(best_transcript)) {
      return(candidate_exons[candidate_exons$transcript == best_transcript, ])
    }
  }

  if (length(unique(candidate_exons$transcript)) > 1) {
    encompassing <- find_encompassing_transcripts(candidate_exons, breakpoint)
    if (any(encompassing)) {
      return(candidate_exons[encompassing, ])
    }
  }

  candidate_exons
}

calculate_transcript_coverage <- function(exons, coverage) {
  transcripts <- unique(exons$transcript)
  sapply(transcripts, function(tx) {
    tx_exons <- exons[exons$transcript == tx, ]
    tx_exons$start <- sapply(tx_exons$start, max, min(start(coverage)))
    tx_exons$end <- sapply(tx_exons$end, min, max(end(coverage)))
    sum(as.numeric(coverage[IRanges(tx_exons$start, tx_exons$end)]))
  })
}

calculate_transcript_lengths <- function(exons) {
  transcripts <- unique(exons$transcript)
  sapply(transcripts, function(tx) {
    tx_exons <- exons[exons$transcript == tx, ]
    sum(tx_exons$end - tx_exons$start + 1)
  })
}

select_best_transcript <- function(coverage, lengths) {
  if (length(coverage) == 0) return(NULL)

  best_coverage <- -1
  best_transcript <- NULL
  best_length <- 0

  for (tx in names(coverage)) {
    tx_length <- lengths[tx]
    tx_coverage <- coverage[tx]
    substantial_diff <- (1 - min(tx_length, best_length) / max(tx_length, best_length)) / 10

    if (best_length == 0) {
      best_transcript <- tx
      best_coverage <- tx_coverage
      best_length <- tx_length
    } else if (tx_length > best_length &&
               tx_coverage * (1 - substantial_diff) > best_coverage) {
      best_transcript <- tx
      best_coverage <- tx_coverage
      best_length <- tx_length
    } else if (tx_length < best_length && tx_coverage > best_coverage * (1 - substantial_diff)) {
      best_transcript <- tx
      best_coverage <- tx_coverage
      best_length <- tx_length
    }
  }
  best_transcript
}

find_encompassing_transcripts <- function(exons, breakpoint) {
  transcript_starts <- aggregate(exons$start, by = list(exons$transcript), min)
  transcript_ends <- aggregate(exons$end, by = list(exons$transcript), max)
  starts <- setNames(transcript_starts$x, transcript_starts$Group.1)
  ends <- setNames(transcript_ends$x, transcript_ends$Group.1)
  starts[exons$transcript] <= breakpoint & ends[exons$transcript] >= breakpoint
}

select_by_coding_length <- function(exons) {
  coding_lengths <- aggregate(
    ifelse(exons$type == "CDS", exons$end - exons$start, 0),
    by = list(exons$transcript), sum
  )
  max_length <- max(coding_lengths$x)
  best_transcripts <- coding_lengths[coding_lengths$x == max_length, 1]
  exons[exons$transcript %in% best_transcripts, ]
}

select_by_total_length <- function(exons) {
  total_lengths <- aggregate(
    exons$end - exons$start,
    by = list(exons$transcript), sum
  )
  max_length <- max(total_lengths$x)
  best_transcripts <- total_lengths[total_lengths$x == max_length, 1]
  exons[exons$transcript %in% best_transcripts, ]
}


# =============================================================================
# 第六部分：findClosestGene()
# =============================================================================

findClosestGene <- function(exons, contig, breakpoint, extraConditions) {
  closest_exons <- exons[exons$contig == contig & extraConditions, ]
  closest_exons <- exons[
    exons$contig == contig &
    exons$geneID %in% closest_exons$geneID,
  ]

  if (length(unique(closest_exons$geneID)) > 1) {
    distance <- aggregate(
      1:nrow(closest_exons),
      by = list(closest_exons$geneID),
      function(idx) {
        min(abs(closest_exons[idx, "start"] - breakpoint),
            abs(closest_exons[idx, "end"] - breakpoint))
      }
    )
    min_dist <- min(distance$x)
    closest_gene <- distance[distance$x == min_dist, 1][1]
    closest_exons <- closest_exons[closest_exons$geneID == closest_gene, ]
  }

  if (nrow(closest_exons) == 0) {
    return(IRanges(max(1, breakpoint - 1000), breakpoint + 1000))
  }

  IRanges(min(closest_exons$start), max(closest_exons$end))
}


# =============================================================================
# 第七部分：绘图函数
# =============================================================================

drawVerticalGradient <- function(left, right, y, color, selection = NULL) {
  if (!getOption("draw_3d_effect", TRUE)) return()

  if (!is.null(selection)) {
    y <- y[selection]
    left <- left[selection]
    right <- right[selection]
  }

  for (i in seq_along(y)) {
    polygon(
      c(left[1:i], right[1:i]),
      c(y[1:i], y[i:1]),
      border = NA,
      col = rgb(
        col2rgb(color)["red", ],
        col2rgb(color)["green", ],
        col2rgb(color)["blue", ],
        col2rgb(color, alpha = TRUE)["alpha", ] * (1 / length(y)),
        maxColorValue = 255
      )
    )
  }
}

drawCurlyBrace <- function(left, right, top, bottom, tip) {
  smoothness <- 20
  x <- cumsum(exp(-seq(-2.5, 2.5, length.out = smoothness)^2))
  x <- x / max(x)
  y <- seq(top, bottom, length.out = smoothness)
  lines(left + (tip - left) + x * (left - tip), y)
  lines(tip + x * (right - tip), y)
}

drawCoverage <- function(left, right, y, coverage, start, end, color, config = NULL) {
  max_resolution <- if (!is.null(config)) config$dimensions$max_resolution else 5000
  if (is.null(coverage)) return()

  coverageData <- as.numeric(
    coverage[IRanges(
      sapply(start, max, min(start(coverage))),
      sapply(end, min, max(end(coverage)))
    )]
  )
  coverageData <- aggregate(
    coverageData,
    by = list(round(seq_along(coverageData) *
                    (right - left) * max_resolution /
                    length(coverageData))),
    mean
  )$x

  polygon(
    c(left, seq(left, right, length.out = length(coverageData)), right),
    c(y, y + coverageData * 0.1, y),
    col = color, border = NA
  )
}

drawStrand <- function(left, right, y, color, strand, config = NULL) {
  if (!(strand %in% c("+", "-"))) return()

  arrow_offset <- if (!is.null(config)) config$dimensions$strand_arrow_offset else 0.001
  arrow_spacing <- if (!is.null(config)) config$dimensions$strand_arrow_spacing else 0.005
  arrow_length <- if (!is.null(config)) config$dimensions$strand_arrow_length_scaled else 0.01

  lines(c(left + arrow_offset, right - arrow_offset), c(y, y), col = color, lwd = 2)
  lines(c(left + arrow_offset, right - arrow_offset), c(y, y), col = rgb(1, 1, 1, 0.1), lwd = 1)

  if (right - left > 0.01) {
    for (i in seq(left + arrow_spacing, right - arrow_spacing,
                   by = sign(right - left - 2 * arrow_spacing) * (2 * arrow_spacing))) {
      direction <- ifelse(strand == "+", 1, -1)
      arrows(i, y, i + arrow_length * direction, y,
             col = color, length = 0.05, lwd = 2, angle = 60)
      arrows(i, y, i + arrow_length * direction, y,
             col = rgb(1, 1, 1, 0.1), length = 0.05, lwd = 1, angle = 60)
    }
  }
}

drawExon <- function(left, right, y, color, title, type,
                    fontSize = 1, gradientSteps = 10, config = NULL) {

  exon_height <- if (!is.null(config)) config$dimensions$exon_height else 0.03
  border_padding <- if (!is.null(config)) config$dimensions$exon_border_padding else 0.001

  if (type == "CDS") {
    rect(left, y + exon_height, right, y + exon_height / 2 - border_padding,
         col = color, border = NA)
    rect(left, y - exon_height, right, y - exon_height / 2 + border_padding,
         col = color, border = NA)

    lines(c(left, left, right, right),
          c(y + exon_height / 2, y + exon_height, y + exon_height, y + exon_height / 2),
          col = getDarkColor(color), lend = 2)
    lines(c(left, left, right, right),
          c(y - exon_height / 2, y - exon_height, y - exon_height, y - exon_height / 2),
          col = getDarkColor(color), lend = 2)

    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps),
                        seq(y + exon_height, y + exon_height / 2, length.out = gradientSteps),
                        rgb(0, 0, 0, 0.2))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps),
                        seq(y - exon_height, y - exon_height / 2, length.out = gradientSteps),
                        rgb(0, 0, 0, 0.3))

  } else if (type == "exon") {
    rect(left, y + exon_height / 2, right, y - exon_height / 2,
         col = color, border = getDarkColor(color))

    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps),
                        seq(y, y + exon_height / 2, length.out = gradientSteps),
                        rgb(1, 1, 1, 0.6))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps),
                        seq(y, y - exon_height / 2, length.out = gradientSteps),
                        rgb(1, 1, 1, 0.6))

    text((left + right) / 2, y, title, cex = 0.9 * fontSize)
  }
}

drawIdeogram <- function(adjust, left, right, y, cytobands, contig, breakpoint, config = NULL) {
  arc_steps <- if (!is.null(config)) config$dimensions$ideogram_arc_steps else 30
  curly_brace_height <- if (!is.null(config)) config$dimensions$curly_brace_height else 0.03
  ideogram_height <- if (!is.null(config)) config$dimensions$ideogram_height else 0.04
  ideogram_width <- if (!is.null(config)) config$dimensions$ideogram_width else 0.4

  bandColors <- setNames(
    rgb(100:0, 100:0, 100:0, maxColorValue = 100),
    paste0("gpos", 0:100)
  )
  bandColors <- c(bandColors,
                  gneg = "#ffffff",
                  acen = "#ec4f4f",
                  stalk = "#0000ff")
  cytobands$color <- bandColors[cytobands$giemsa]

  # 标准化 contig 名称以匹配 cytobands
  contig <- normalizeChromosome(contig)
  bands <- cytobands[cytobands$contig == contig, ]
  if (nrow(bands) == 0) {
    warning(paste("Ideogram of contig", contig, "cannot be drawn"))
    return(NULL)
  }

  bands$left <- bands$start / max(cytobands$end) * ideogram_width
  bands$right <- bands$end / max(cytobands$end) * ideogram_width

  offset <- ifelse(adjust == "left", left, right - max(bands$right))
  bands$left <- bands$left + offset
  bands$right <- bands$right + offset

  tip <- min(bands$left) +
         (max(bands$right) - min(bands$left)) /
         (max(bands$end) - min(bands$start)) * breakpoint
  drawCurlyBrace(left, right, y - 0.05 + curly_brace_height, y - 0.05, tip)

  text((max(bands$right) + min(bands$left)) / 2, y + 0.07,
       paste("chromosome", contig), font = 2, cex = 1)

  bandName <- bands[between(breakpoint, bands$start, bands$end), "name"]
  text(tip, y + 0.03, bandName, cex = 1, adj = c(0.5, 0))

  band <- 1
  leftArcX <- bands[band, "left"] +
              (1 + cos(seq(pi / 2, 1.5 * pi, length.out = arc_steps))) *
              (bands[band, "right"] - bands[band, "left"])
  leftArcY <- y + sin(seq(pi / 2, 1.5 * pi, length.out = arc_steps)) *
              (ideogram_height / 2)
  polygon(leftArcX, leftArcY, col = bands[band, "color"])

  centromere_start <- NULL
  centromere_end <- NULL

  for (bandIdx in 2:(nrow(bands) - 1)) {
    if (bands[bandIdx, "giemsa"] != "acen") {
      rect(bands[bandIdx, "left"], y - ideogram_height / 2,
           bands[bandIdx, "right"], y + ideogram_height / 2,
           col = bands[bandIdx, "color"])
    } else {
      if (is.null(centromere_start)) {
        polygon(c(bands[bandIdx, "left"], bands[bandIdx, "right"], bands[bandIdx, "left"]),
                c(y - ideogram_height / 2, y, y + ideogram_height / 2),
                col = bands[bandIdx, "color"])
        centromere_start <- bands[bandIdx, "left"]
      } else {
        polygon(c(bands[bandIdx, "right"], bands[bandIdx, "left"], bands[bandIdx, "right"]),
                c(y - ideogram_height / 2, y, y + ideogram_height / 2),
                col = bands[bandIdx, "color"])
        centromere_end <- bands[bandIdx, "right"]
      }
    }
  }

  band <- nrow(bands)
  rightArcX <- bands[band, "right"] -
               (1 + cos(seq(1.5 * pi, pi / 2, length.out = arc_steps))) *
               (bands[band, "right"] - bands[band, "left"])
  rightArcY <- y + sin(seq(pi / 2, 1.5 * pi, length.out = arc_steps)) *
               ideogram_height / 2
  polygon(rightArcX, rightArcY, col = bands[band, "color"])

  if (is.null(centromere_start) || is.null(centromere_end)) {
    centromere_start <- bands[1, "right"]
    centromere_end <- bands[1, "right"]
  }

  drawVerticalGradient(leftArcX, rep(centromere_start, arc_steps), leftArcY,
                      rgb(0, 0, 0, 0.8), 1:round(arc_steps * 0.4))
  drawVerticalGradient(leftArcX, rep(centromere_start, arc_steps), leftArcY,
                      rgb(1, 1, 1, 0.7), round(arc_steps * 0.4):round(arc_steps * 0.1))
  drawVerticalGradient(leftArcX, rep(centromere_start, arc_steps), leftArcY,
                      rgb(1, 1, 1, 0.7), round(arc_steps * 0.4):round(arc_steps * 0.6))
  drawVerticalGradient(leftArcX, rep(centromere_start, arc_steps), leftArcY,
                      rgb(0, 0, 0, 0.9), arc_steps:round(arc_steps * 0.5))
  drawVerticalGradient(rightArcX, rep(centromere_end, arc_steps), rightArcY,
                      rgb(0, 0, 0, 0.8), 1:round(arc_steps * 0.4))
  drawVerticalGradient(rightArcX, rep(centromere_end, arc_steps), rightArcY,
                      rgb(1, 1, 1, 0.7), round(arc_steps * 0.4):round(arc_steps * 0.1))
  drawVerticalGradient(rightArcX, rep(centromere_end, arc_steps), rightArcY,
                      rgb(1, 1, 1, 0.7), round(arc_steps * 0.4):round(arc_steps * 0.6))
  drawVerticalGradient(rightArcX, rep(centromere_end, arc_steps), rightArcY,
                      rgb(0, 0, 0, 0.9), arc_steps:round(arc_steps * 0.5))
}


# =============================================================================
# 第八部分：蛋白质域处理函数
# =============================================================================

removeIntronsFromProteinDomains <- function(codingExons, retainedDomains) {
  if (nrow(codingExons) == 0) return(NULL)
  cumulativeIntronLength <- 0
  previousExonEnd <- 0

  for (exon in seq_len(nrow(codingExons))) {
    if (codingExons[exon, "start"] > previousExonEnd)
      cumulativeIntronLength <- cumulativeIntronLength +
                               codingExons[exon, "start"] - previousExonEnd
    domainsInExon <- which(
      between(retainedDomains$start,
              codingExons[exon, "start"],
              codingExons[exon, "end"])
    )
    retainedDomains[domainsInExon, "start"] <-
      retainedDomains[domainsInExon, "start"] - cumulativeIntronLength
    domainsInExon <- which(
      between(retainedDomains$end,
              codingExons[exon, "start"],
              codingExons[exon, "end"])
    )
    retainedDomains[domainsInExon, "end"] <-
      retainedDomains[domainsInExon, "end"] - cumulativeIntronLength
    previousExonEnd <- codingExons[exon, "end"]
  }

  if (nrow(retainedDomains) > 0) {
    retainedDomains <- do.call(rbind, lapply(
      unique(retainedDomains$proteinDomainID),
      function(x) {
        domain <- retainedDomains[retainedDomains$proteinDomainID == x, ]
        merged <- reduce(GRanges(
          domain$seqnames,
          IRanges(domain$start, domain$end),
          strand = domain$strand
        ))
        merged$proteinDomainName <- head(domain$proteinDomainName, 1)
        merged$proteinDomainID <- head(domain$proteinDomainID, 1)
        merged$color <- head(domain$color, 1)
        as.data.frame(merged)
      }
    ))
  }
  retainedDomains
}


# =============================================================================
# 第九部分：参数定义与解析
# =============================================================================

parameters <- list(
  # 文件输入模式 (二选一)
  fusions = list("fusionsFile", "file", ""),

  # 必需参数
  annotation = list("exonsFile", "file", "annotation.gtf", TRUE),
  output = list("outputFile", "string", "output.pdf", TRUE),

  # 可选参数
  alignments = list("alignmentsFile", "file", "Aligned.sortedByCoord.out.bam"),
  cytobands = list("cytobandsFile", "file", "cytobands.tsv"),
  proteinDomains = list("proteinDomainsFile", "file", "protein_domains.gff3"),
  sampleName = list("sampleName", "string", ""),
  minConfidenceForCircosPlot = list("minConfidenceForCircosPlot", "string", "medium"),
  plotPanels = list("plotPanels", "string", "fusion,circos,domains,readcounts"),
  pdfWidth = list("pdfWidth", "numeric", 11.692),
  pdfHeight = list("pdfHeight", "numeric", 8.267),
  color1 = list("color1", "string", "#e5a5a5"),
  color2 = list("color2", "string", "#a7c4e5"),
  fontSize = list("fontSize", "numeric", 1),
  fontFamily = list("fontFamily", "string", "Helvetica"),
  squishIntrons = list("squishIntrons", "bool", TRUE),
  printExonLabels = list("printExonLabels", "bool", TRUE),
  render3dEffect = list("render3dEffect", "bool", TRUE),
  showIntergenicVicinity = list("showIntergenicVicinity", "string", "0"),
  transcriptSelection = list("transcriptSelection", "string", "provided"),
  fixedScale = list("fixedScale", "numeric", 0),
  coverageRange = list("coverageRange", "string", "0"),
  mergeDomainsOverlappingBy = list("mergeDomainsOverlappingBy", "numeric", 0.9),
  optimizeDomainColors = list("optimizeDomainColors", "bool", FALSE),

  # 单融合模式参数
  gene1 = list("gene1", "string", "", TRUE),
  gene2 = list("gene2", "string", "", TRUE),
  contig1 = list("contig1", "string", "", TRUE),
  contig2 = list("contig2", "string", "", TRUE),
  breakpoint1 = list("breakpoint1", "numeric", 0, TRUE),
  breakpoint2 = list("breakpoint2", "numeric", 0, TRUE),
  direction1 = list("direction1", "string", "downstream", TRUE),
  direction2 = list("direction2", "string", "downstream", TRUE),
  strand1 = list("strand1", "string", "."),
  strand2 = list("strand2", "string", "."),
  site1 = list("site1", "string", "exon"),
  site2 = list("site2", "string", "exon"),
  transcript_id1 = list("transcript_id1", "string", "."),
  transcript_id2 = list("transcript_id2", "string", "."),
  reading_frame = list("reading_frame", "string", "."),
  fusion_transcript = list("fusion_transcript", "string", "."),
  split_reads = list("split_reads", "numeric", 0),
  discordant_mates = list("discordant_mates", "numeric", 0),
  confidence = list("confidence", "string", "high")
)


# =============================================================================
# 第十部分：主程序
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (any(grepl("^--help", args)) || length(args) == 0) {
  usage <- "Usage: draw_fusions.R"
  for (parameter in names(parameters)) {
    usage <- paste0(usage, " ")
    if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
      usage <- paste0(usage, "[")
    usage <- paste0(usage, "--", parameter, "=", parameters[[parameter]][[3]])
    if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
      usage <- paste0(usage, "]")
  }
  message(usage)
  quit("no", ifelse(length(args) == 0, 1, 0))
}

for (parameter in names(parameters))
  if (length(parameters[[parameter]]) > 3 && parameters[[parameter]][[4]])
    if (!any(grepl(paste0("^--", parameter, "="), args), perl = TRUE))
      stop(paste0("Missing mandatory argument: --", parameter))

for (parameter in names(parameters))
  assign(parameters[[parameter]][[1]],
         ifelse(parameters[[parameter]][[2]] == "file", "",
                parameters[[parameter]][[3]]))

for (arg in args) {
  argName <- sub("=.*", "", sub("^--", "", arg, perl = TRUE), perl = TRUE)
  argValue <- sub("^[^=]*=", "", arg, perl = TRUE)
  if (!(argName %in% names(parameters)) || !grepl("^--", arg, perl = TRUE))
    stop(paste("Unknown parameter:", arg))
  if (parameters[[argName]][[2]] == "bool") {
    if (argValue %in% c("TRUE", "T", "FALSE", "F")) {
      assign(parameters[[argName]][[1]], as.logical(argValue))
    } else {
      stop(paste0("Invalid argument to --", argName))
    }
  } else if (parameters[[argName]][[2]] == "string") {
    assign(parameters[[argName]][[1]], argValue)
  } else if (parameters[[argName]][[2]] == "numeric") {
    if (is.na(suppressWarnings(as.numeric(argValue))))
      stop(paste0("Invalid argument to --", argName))
    assign(parameters[[argName]][[1]], as.numeric(argValue))
  } else if (parameters[[argName]][[2]] == "file") {
    if (file.access(argValue) == -1)
      stop(paste("Cannot read file:", argValue))
    assign(parameters[[argName]][[1]], argValue)
  }
}

# 验证输入模式 (文件模式 vs 单融合模式)
is_file_mode <- fusionsFile != ""
is_single_fusion_mode <- gene1 != "" && gene2 != ""

if (!is_file_mode && !is_single_fusion_mode) {
  stop("Missing input: Either --fusions file or single fusion parameters (--gene1, --gene2, etc.) must be provided")
}

if (is_file_mode && is_single_fusion_mode) {
  stop("Conflicting input: Cannot use both --fusions file and single fusion parameters at the same time")
}

# 自动推断单融合模式的 fusion_type
if (is_single_fusion_mode) {
  if (contig1 != contig2) {
    fusion_type <- "translocation"
  } else {
    # 同一条染色体上
    if (direction1 == direction2) {
      fusion_type <- "inversion"
    } else if ((direction1 == "downstream" && breakpoint1 < breakpoint2) ||
               (direction1 == "upstream" && breakpoint1 > breakpoint2)) {
      fusion_type <- "deletion"
    } else {
      fusion_type <- "duplication"
    }
  }
}

if (cytobandsFile == "")
  warning("Missing parameter '--cytobands'. No ideograms and circos plots will be drawn.")
if (!(minConfidenceForCircosPlot %in% c("none", "low", "medium", "high")))
  stop("Invalid argument to --minConfidenceForCircosPlot")
showIntergenicVicinity <- as.list(unlist(strsplit(showIntergenicVicinity, ",", fixed = TRUE)))
if (!(length(showIntergenicVicinity) %in% c(1, 4)))
  stop(paste0("Invalid argument to --showIntergenicVicinity"))
showIntergenicVicinity <- lapply(showIntergenicVicinity, function(x) {
  if (x == "closest_gene") {
    return("exon")
  } else if (x == "closestProteinCodingGene") {
    return("CDS")
  } else if (is.na(suppressWarnings(as.numeric(x))) || as.numeric(x) < 0) {
    stop(paste0("Invalid argument to --showIntergenicVicinity"))
  } else {
    return(as.numeric(x))
  }
})
if (length(showIntergenicVicinity) == 1)
  showIntergenicVicinity <- rep(showIntergenicVicinity, 4)
if (squishIntrons)
  if (any(!is.numeric(unlist(showIntergenicVicinity))) || any(showIntergenicVicinity > 0))
    stop("--squishIntrons must be disabled, when --showIntergenicVicinity is > 0")
if (!(transcriptSelection %in% c("coverage", "provided", "canonical")))
  stop("Invalid argument to --transcriptSelection")
if (fixedScale < 0)
  stop("Invalid argument to --fixedScale")
if (!(fontFamily %in% names(pdfFonts())))
  stop(paste0("Unknown font: ", fontFamily, ". Available fonts: ", paste(names(pdfFonts()), collapse = ", ")))
coverageRange <- suppressWarnings(as.numeric(unlist(strsplit(coverageRange, ",", fixed = TRUE))))
if (!(length(coverageRange) %in% 1:2) || any(is.na(coverageRange)) || any(coverageRange < 0))
  stop("Invalid argument to --coverageRange")
plotPanels <- unlist(strsplit(plotPanels, split = ","))
invalidPlotPanels <- !(plotPanels %in% unlist(strsplit(parameters[["plotPanels"]][[3]], split = ",")))
if (any(invalidPlotPanels))
  stop(paste("Invalid argument to --plotPanels:", paste(plotPanels[invalidPlotPanels], collapse = ",")))

if (!suppressPackageStartupMessages(require(GenomicRanges)))
  warning("Package 'GenomicRanges' is not installed. No protein domains and circos plots will be drawn.")
if (!suppressPackageStartupMessages(require(circlize)))
  warning("Package 'circlize' is not installed. No circos plots will be drawn.")
if (alignmentsFile != "")
  if (!suppressPackageStartupMessages(require(GenomicAlignments)))
    stop("Package 'GenomicAlignments' must be installed when '--alignments' is used")

dark_color_gene1 <- getDarkColor(color1)
dark_color_gene2 <- getDarkColor(color2)
circosColors <- c(
  translocation = "#000000",
  duplication = "#00bb00",
  deletion = "#ff0000",
  inversion = "#0000ff"
)


# =============================================================================
# 第十一部分：读取融合数据
# =============================================================================

if (fusionsFile != "") {
  fusions <- read.table(fusionsFile, stringsAsFactors = FALSE,
                       sep = "\t", header = TRUE,
                       comment.char = "", quote = "")
  if (colnames(fusions)[1] == "X.gene1") {
    colnames(fusions)[colnames(fusions) %in%
                      c("X.gene1", "strand1.gene.fusion.", "strand2.gene.fusion.")] <-
      c("gene1", "strand1", "strand2")
    fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl = TRUE)
    fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl = TRUE)
    fusions$contig1 <- removeChr(fusions$display_contig1)
    fusions$contig2 <- removeChr(fusions$display_contig2)
    fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl = TRUE))
    fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl = TRUE))
    fusions$split_reads1 <- fusions$split_reads1
    fusions$split_reads2 <- fusions$split_reads2
    fusions$type <- sub(".*(translocation|duplication|deletion|inversion).*",
                        "\\1", fusions$type, perl = TRUE)
    fusions$fusion_transcript <- gsub("[()^$]", "", fusions$fusion_transcript)
  } else if (colnames(fusions)[1] == "X.FusionName") {
    fusions$gene1 <- sub("\\^.*", "", fusions$LeftGene, perl = TRUE)
    fusions$gene2 <- sub("\\^.*", "", fusions$RightGene, perl = TRUE)
    fusions$strand1 <- sub(".*:(.)$", "\\1/\\1", fusions$LeftBreakpoint, perl = TRUE)
    fusions$strand2 <- sub(".*:(.)$", "\\1/\\1", fusions$RightBreakpoint, perl = TRUE)
    fusions$display_contig1 <- sub(":[^:]*:[^:]*$", "", fusions$LeftBreakpoint, perl = TRUE)
    fusions$display_contig2 <- sub(":[^:]*:[^:]*$", "", fusions$RightBreakpoint, perl = TRUE)
    fusions$contig1 <- removeChr(fusions$display_contig1)
    fusions$contig2 <- removeChr(fusions$display_contig2)
    fusions$breakpoint1 <- as.numeric(
      sub(".*:([^:]*):[^:]*$", "\\1", fusions$LeftBreakpoint, perl = TRUE))
    fusions$breakpoint2 <- as.numeric(
      sub(".*:([^:]*):[^:]*$", "\\1", fusions$RightBreakpoint, perl = TRUE))
    fusions$direction1 <- ifelse(grepl(":\\+$", fusions$LeftBreakpoint, perl = TRUE),
                               "downstream", "upstream")
    fusions$direction2 <- ifelse(grepl(":\\+$", fusions$RightBreakpoint, perl = TRUE),
                               "upstream", "downstream")
    fusions$gene_id1 <- sub(".*\\^", "", fusions$LeftGene, perl = TRUE)
    fusions$gene_id2 <- sub(".*\\^", "", fusions$RightGene, perl = TRUE)
    fusions$transcript_id1 <- ifelse(
      !("CDS_LEFT_ID" %in% colnames(fusions)), ".",
      fusions$CDS_LEFT_ID)
    fusions$transcript_id2 <- ifelse(
      !("CDS_RIGHT_ID" %in% colnames(fusions)), ".",
      fusions$CDS_RIGHT_ID)
    fusions$fusion_transcript <- ifelse(
      !("FUSION_CDS" %in% colnames(fusions)), ".",
      toupper(sub("([a-z]*)", "\\1|", fusions$FUSION_CDS, perl = TRUE)))
    fusions$reading_frame <- ifelse(
      !("PROT_FUSION_TYPE" %in% colnames(fusions)), ".",
      ifelse(fusions$PROT_FUSION_TYPE == "INFRAME", "in-frame",
            ifelse(fusions$PROT_FUSION_TYPE == "FRAMESHIFT", "out-of-frame", ".")))
    fusions$split_reads <- fusions$JunctionReadCount
    fusions$discordant_mates <- fusions$SpanningFragCount
    fusions$site1 <- rep("exon", nrow(fusions))
    fusions$site2 <- rep("exon", nrow(fusions))
    fusions$confidence <- rep("high", nrow(fusions))
    fusions$type <- ifelse(
      fusions$contig1 != fusions$contig2, "translocation",
      ifelse(fusions$direction1 == fusions$direction2, "inversion",
            ifelse((fusions$direction1 == "downstream") ==
                   (fusions$breakpoint1 < fusions$breakpoint2),
                   "deletion", "duplication")))
  } else {
    stop("Unrecognized fusion file format")
  }
} else if (gene1 != "" && contig1 != "" && breakpoint1 > 0) {
  fusion <- data.frame(
    gene1 = gene1,
    gene2 = gene2,
    gene_id1 = gene1,
    gene_id2 = gene2,
    contig1 = contig1,
    contig2 = contig2,
    breakpoint1 = breakpoint1,
    breakpoint2 = breakpoint2,
    direction1 = direction1,
    direction2 = direction2,
    type = fusion_type,
    strand1 = strand1,
    strand2 = strand2,
    site1 = site1,
    site2 = site2,
    transcript_id1 = transcript_id1,
    transcript_id2 = transcript_id2,
    reading_frame = reading_frame,
    fusion_transcript = fusion_transcript,
    split_reads1 = split_reads,
    split_reads2 = split_reads,
    split_reads = split_reads,
    discordant_mates = discordant_mates,
    confidence = confidence,
    stringsAsFactors = FALSE
  )
  fusion$display_contig1 <- ifelse(grepl("^chr", contig1), contig1, paste0("chr", contig1))
  fusion$display_contig2 <- ifelse(grepl("^chr", contig2), contig2, paste0("chr", contig2))
  fusions <- fusion
} else {
  stop("Must specify either --fusions FILE or all required single fusion parameters (--gene1, --contig1, --breakpoint1, etc.)")
}


pdf(outputFile, onefile = TRUE, width = pdfWidth, height = pdfHeight,
    title = ifelse(sampleName != "", sampleName, fusionsFile))
par(family = fontFamily)

if (nrow(fusions) == 0) {
  plot(0, 0, type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  text(0, 0, "empty input file")
  warning("empty input file")
  dev.off()
  quit("no")
}


# =============================================================================
# 第十二部分：读取注释数据
# =============================================================================

cytobands <- NULL
if (cytobandsFile != "") {
  # 检测文件格式
  is_bed <- grepl("\\.bed(\\.gz)?$", cytobandsFile, perl = TRUE)

  if (is_bed) {
    # BED 格式（支持 .bed 和 .bed.gz）
    cytobands <- tryCatch({
      read.table(gzfile(cytobandsFile), header = FALSE, sep = "\t")
    }, error = function(e) {
      read.table(cytobandsFile, header = FALSE, sep = "\t")
    })

    # 统一列名（BED 是 0-based，需要 +1）
    # 标准 BED 格式至少需要 3 列，最多 12 列
    # cytoBand 通常有：contig, start, end, name, giemsa (或 strand)
    n_cols <- ncol(cytobands)
    colnames(cytobands)[1:3] <- c("contig", "start", "end")

    # 尝试自动识别 giemsa 列（包含 gneg/gpos 等关键词）
    giemsa_col <- which(sapply(cytobands, function(x) {
      any(grepl("^g(pos|neg)", x, ignore.case = TRUE))
    }))
    if (length(giemsa_col) > 0) {
      cytobands$giemsa <- cytobands[[giemsa_col[1]]]
    } else {
      cytobands$giemsa <- ifelse(grepl("^chr", cytobands$contig),
                                 sub("^chr", "", cytobands$contig), cytobands$contig)
    }

    cytobands$start <- as.numeric(cytobands$start) + 1  # BED 0-based 转 1-based
    cytobands$contig <- removeChr(cytobands$contig)
  } else {
    # TSV 格式
    cytobands <- read.table(cytobandsFile, header = TRUE,
                            colClasses = c("character", "numeric", "numeric",
                                          "character", "character"))
  }
  cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end), ]
}
if (is.null(cytobands) || !("circlize" %in% names(sessionInfo()$otherPkgs)) ||
    !("GenomicRanges" %in% names(sessionInfo()$otherPkgs)))
  plotPanels <- setdiff(plotPanels, "circos")

message("Loading annotation")

# 检测是否为 bgzip 压缩文件
is_gzipped <- grepl("\\.gz$", exonsFile, perl = TRUE)

# 读取 GTF 文件（支持 bgzip 压缩）
if (is_gzipped) {
  exons <- tryCatch({
    read.table(gzfile(exonsFile),
               header = FALSE,
               sep = "\t",
               comment.char = "#",
               quote = '"',
               stringsAsFactors = FALSE)
  }, error = function(e) {
    message("Failed to read gzipped GTF, falling back to uncompressed reading")
    read.table(exonsFile,
               header = FALSE,
               sep = "\t",
               comment.char = "#",
               quote = '"',
               stringsAsFactors = FALSE)
  })
  colnames(exons) <- c("contig", "src", "type", "start", "end", "score",
                       "strand", "frame", "attributes")
} else {
  exons <- scan(exonsFile,
                what = list(contig = "", src = "", type = "",
                           start = 0, end = 0, score = "",
                           strand = "", frame = "",
                           attributes = ""),
                sep = "\t", comment.char = "#", quote = '"',
                multi.line = FALSE)
  attr(exons, "row.names") <- .set_row_names(length(exons[[1]]))
  class(exons) <- "data.frame"
}

# 过滤外显子和 CDS
exons <- exons[exons$type %in% c("exon", "CDS"),
              c("contig", "type", "start", "end", "strand", "attributes")]
exons$contig <- removeChr(exons$contig)

parseGtfAttribute <- function(attribute, gtf) {
  parsed <- sub(
    paste0(".*", attribute, "[ =]([^;]+).*"),
    "\\1", gtf$attributes, perl = TRUE
  )
  failedToParse <- parsed == gtf$attributes
  if (any(failedToParse)) {
    warning(paste0("Failed to parse '", attribute, "' attribute of ",
                   sum(failedToParse), " record(s)."))
    parsed <- ifelse(failedToParse, "", parsed)
  }
  return(parsed)
}

exons$geneID <- parseGtfAttribute("gene_id", exons)
exons$geneName <- parseGtfAttribute("gene_name", exons)
exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
exons$transcript <- parseGtfAttribute("transcript_id", exons)
exons$exonNumber <- ifelse(rep(printExonLabels, nrow(exons)),
                          parseGtfAttribute("exon_number", exons), "")

proteinDomains <- NULL
if (proteinDomainsFile != "") {
  message("Loading protein domains")
  proteinDomains <- scan(proteinDomainsFile,
                        what = list(contig = "", src = "", type = "",
                                   start = 0, end = 0, score = "",
                                   strand = "", frame = "",
                                   attributes = ""),
                        sep = "\t", comment.char = "", quote = "",
                        multi.line = FALSE)
  attr(proteinDomains, "row.names") <- .set_row_names(length(proteinDomains[[1]]))
  class(proteinDomains) <- "data.frame"
  proteinDomains$color <- parseGtfAttribute("color", proteinDomains)
  proteinDomains$proteinDomainName <- sapply(
    parseGtfAttribute("Name", proteinDomains), URLdecode)
  proteinDomains$proteinDomainID <- parseGtfAttribute("protein_domain_id", proteinDomains)
}
if (is.null(proteinDomains) && "GenomicRanges" %in% names(sessionInfo()$otherPkgs))
  plotPanels <- setdiff(plotPanels, "domains")

if (any(fusions$site1 == "intergenic" | fusions$site2 == "intergenic")) {
  intergenicBreakpoints <- rbind(
    setNames(fusions[fusions$site1 == "intergenic",
                    c("gene1", "strand1", "contig1", "breakpoint1")],
             c("gene", "strand", "contig", "breakpoint")),
    setNames(fusions[fusions$site2 == "intergenic",
                    c("gene2", "strand2", "contig2", "breakpoint2")],
             c("gene", "strand", "contig", "breakpoint"))
  )
  exons <- rbind(exons, data.frame(
    contig = intergenicBreakpoints$contig,
    type = "intergenic",
    start = sapply(intergenicBreakpoints$breakpoint - 1000, max, 1),
    end = intergenicBreakpoints$breakpoint + 1000,
    strand = ".",
    attributes = "",
    geneName = intergenicBreakpoints$gene,
    geneID = paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
    transcript = paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
    exonNumber = "intergenic"
  ))
  fusions[fusions$site1 == "intergenic", "gene_id1"] <-
    paste0(fusions[fusions$site1 == "intergenic", "contig1"], ":",
           fusions[fusions$site1 == "intergenic", "breakpoint1"])
  fusions[fusions$site2 == "intergenic", "gene_id2"] <-
    paste0(fusions[fusions$site2 == "intergenic", "contig2"], ":",
           fusions[fusions$site2 == "intergenic", "breakpoint2"])
}


# =============================================================================
# 第十三部分：主循环
# =============================================================================

options(draw_3d_effect = render3dEffect)

for (fusion in seq_len(nrow(fusions))) {

  message(paste0("Drawing fusion #", fusion, ": ",
                 fusions[fusion, "gene1"], ":", fusions[fusion, "gene2"]))

  showVicinity <- rep(0, 4)

  if (fusions[fusion, "site1"] == "intergenic") {
    showVicinity[1] <- ifelse(
      is.numeric(showIntergenicVicinity[[1]]),
      showIntergenicVicinity[[1]],
      fusions[fusion, "breakpoint1"] -
        start(findClosestGene(exons, fusions[fusion, "contig1"],
                             fusions[fusion, "breakpoint1"],
                             exons$end < fusions[fusion, "breakpoint1"] &
                             exons$type == showIntergenicVicinity[[1]]))
    )
    showVicinity[2] <- ifelse(
      is.numeric(showIntergenicVicinity[[2]]),
      showIntergenicVicinity[[2]],
      end(findClosestGene(exons, fusions[fusion, "contig1"],
                         fusions[fusion, "breakpoint1"],
                         exons$start > fusions[fusion, "breakpoint1"] &
                         exons$type == showIntergenicVicinity[[2]])) -
        fusions[fusion, "breakpoint1"]
    )
  }
  if (fusions[fusion, "site2"] == "intergenic") {
    showVicinity[3] <- ifelse(
      is.numeric(showIntergenicVicinity[[3]]),
      showIntergenicVicinity[[3]],
      fusions[fusion, "breakpoint2"] -
        start(findClosestGene(exons, fusions[fusion, "contig2"],
                             fusions[fusion, "breakpoint2"],
                             exons$end < fusions[fusion, "breakpoint2"] &
                             exons$type == showIntergenicVicinity[[3]]))
    )
    showVicinity[4] <- ifelse(
      is.numeric(showIntergenicVicinity[[4]]),
      showIntergenicVicinity[[4]],
      end(findClosestGene(exons, fusions[fusion, "contig2"],
                         fusions[fusion, "breakpoint2"],
                         exons$start > fusions[fusion, "breakpoint2"] &
                         exons$type == showIntergenicVicinity[[4]])) -
        fusions[fusion, "breakpoint2"]
    )
  }

  coverage1 <- NULL
  coverage2 <- NULL

  if (alignmentsFile != "") {
    determineCoverageRegion <- function(exons, geneID, contig, breakpoint,
                                        showVicinityLeft, showVicinityRight) {
      closest_gene <- findClosestGene(exons, contig, breakpoint,
                                    exons$geneID == geneID)
      IRanges(min(start(closest_gene), breakpoint - showVicinityLeft),
             max(end(closest_gene), breakpoint + showVicinityRight))
    }

    coverageRegion1 <- determineCoverageRegion(
      exons, fusions[fusion, "gene_id1"], fusions[fusion, "contig1"],
      fusions[fusion, "breakpoint1"], showVicinity[1], showVicinity[2])
    coverageRegion2 <- determineCoverageRegion(
      exons, fusions[fusion, "gene_id2"], fusions[fusion, "contig2"],
      fusions[fusion, "breakpoint2"], showVicinity[3], showVicinity[4])

    readCoverage <- function(alignmentsFile, contig, coverageRegion) {
      coverageData <- tryCatch({
        alignments <- readGAlignments(
          alignmentsFile,
          param = ScanBamParam(which = GRanges(contig, coverageRegion)))
        coverage(alignments)[[contig]]
      }, error = function(e) {
        alignments <- readGAlignments(
          alignmentsFile,
          param = ScanBamParam(which = GRanges(addChr(contig), coverageRegion)))
        coverage(alignments)[[addChr(contig)]]
      })
      if (exists("alignments")) rm(alignments)
      return(coverageData)
    }

    coverage1 <- readCoverage(alignmentsFile, fusions[fusion, "contig1"],
                             coverageRegion1)
    coverage2 <- readCoverage(alignmentsFile, fusions[fusion, "contig2"],
                             coverageRegion2)

    coverageRegion1 <- IRanges(
      max(start(coverageRegion1), min(start(coverage1))),
      min(end(coverageRegion1), max(end(coverage1))))
    coverageRegion2 <- IRanges(
      max(start(coverageRegion2), min(start(coverage2))),
      min(end(coverageRegion2), max(end(coverage2))))
  }

  exons1 <- findExons(exons, fusions[fusion, "contig1"], fusions[fusion, "gene_id1"],
                     fusions[fusion, "direction1"], fusions[fusion, "breakpoint1"],
                     coverage1, fusions[fusion, "transcript_id1"], transcriptSelection)
  if (nrow(exons1) == 0) {
    par(mfrow = c(1, 1))
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(0, 0, paste0("exon coordinates of ", fusions[fusion, "gene1"],
                       " not found in\n", exonsFile))
    warning(paste("exon coordinates of", fusions[fusion, "gene1"], "not found"))
    next
  }

  exons2 <- findExons(exons, fusions[fusion, "contig2"], fusions[fusion, "gene_id2"],
                     fusions[fusion, "direction2"], fusions[fusion, "breakpoint2"],
                     coverage2, fusions[fusion, "transcript_id2"], transcriptSelection)
  if (nrow(exons2) == 0) {
    par(mfrow = c(1, 1))
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(0, 0, paste0("exon coordinates of ", fusions[fusion, "gene2"],
                       " not found in\n", exonsFile))
    warning(paste("exon coordinates of", fusions[fusion, "gene2"], "not found"))
    next
  }

  if (sum(showVicinity) > 0) {
    if (fusions[fusion, "site1"] == "intergenic") {
      for (geneID in unique(exons[
        exons$contig == fusions[fusion, "contig1"] &
        exons$exonNumber != "intergenic" &
        (between(exons$end, fusions[fusion, "breakpoint1"] - showVicinity[1],
                 fusions[fusion, "breakpoint1"] + showVicinity[2]) |
         between(exons$start, fusions[fusion, "breakpoint1"] - showVicinity[1],
                 fusions[fusion, "breakpoint1"] + showVicinity[2])),"geneID"]))
        exons1 <- rbind(exons1, findExons(
          exons, fusions[fusion, "contig1"], geneID,
          fusions[fusion, "direction1"], fusions[fusion, "breakpoint1"],
          coverage1, fusions[fusion, "transcript_id1"], transcriptSelection))
      exons1 <- exons1[
        (exons1$start >= fusions[fusion, "breakpoint1"] - showVicinity[1] &
         exons1$end <= fusions[fusion, "breakpoint1"] + showVicinity[2]) |
        exons1$exonNumber == "intergenic", ]
    }
    if (fusions[fusion, "site2"] == "intergenic") {
      for (geneID in unique(exons[
        exons$contig == fusions[fusion, "contig2"] &
        exons$exonNumber != "intergenic" &
        (between(exons$end, fusions[fusion, "breakpoint2"] - showVicinity[3],
                 fusions[fusion, "breakpoint2"] + showVicinity[4]) |
         between(exons$start, fusions[fusion, "breakpoint2"] - showVicinity[3],
                 fusions[fusion, "breakpoint2"] + showVicinity[4])),"geneID"]))
        exons2 <- rbind(exons2, findExons(
          exons, fusions[fusion, "contig2"], geneID,
          fusions[fusion, "direction2"], fusions[fusion, "breakpoint2"],
          coverage2, fusions[fusion, "transcript_id2"], transcriptSelection))
      exons2 <- exons2[
        (exons2$start >= fusions[fusion, "breakpoint2"] - showVicinity[3] &
         exons2$end <= fusions[fusion, "breakpoint2"] + showVicinity[4]) |
        exons2$exonNumber == "intergenic", ]
    }
  }

  if (alignmentsFile != "") {
    coverageNormalization <- function(coverage, coverageRegion, exons) {
      max(1, ifelse(
        squishIntrons,
        max(as.numeric(coverage[IRanges(
          sapply(exons$start, max, min(start(coverage))),
          sapply(exons$end, min, max(end(coverage))))])),
        round(quantile(coverage[coverageRegion], 0.9999))
      ))
    }
    coverageNormalization1 <- ifelse(
      head(coverageRange, 1) == 0,
      coverageNormalization(coverage1, coverageRegion1, exons1),
      head(coverageRange, 1))
    coverageNormalization2 <- ifelse(
      tail(coverageRange, 1) == 0,
      coverageNormalization(coverage2, coverageRegion2, exons2),
      tail(coverageRange, 1))
    if (length(coverageRange) == 1 && coverageRange[1] == 0) {
      coverageNormalization1 <- max(coverageNormalization1, coverageNormalization2)
      coverageNormalization2 <- max(coverageNormalization1, coverageNormalization2)
    }
    coverage1 <- coverage1 / coverageNormalization1
    coverage2 <- coverage2 / coverageNormalization2
    coverage1[coverage1 > 1] <- 1
    coverage2[coverage2 > 1] <- 1
  }

  exons1 <- exons1[order(exons1$start, -rank(exons1$type)), ]
  exons2 <- exons2[order(exons2$start, -rank(exons2$type)), ]

  breakpoint1 <- fusions[fusion, "breakpoint1"]
  breakpoint2 <- fusions[fusion, "breakpoint2"]

  if (breakpoint1 < min(exons1$start)) {
    exons1 <- rbind(
      data.frame(contig = fusions[fusion, "contig1"], type = "dummy",
                 start = max(1, breakpoint1 - 1000),
                 end = max(1, breakpoint1 - 1000),
                 strand = exons1[1, "strand"], attributes = "",
                 geneID = exons1[1, "geneID"],
                 transcript = exons1[1, "transcript"],
                 exonNumber = "dummy",
                 left = max(1, breakpoint1 - 1000),
                 right = max(1, breakpoint1 - 1000)),
      exons1)
  } else if (breakpoint1 > max(exons1$end)) {
    exons1 <- rbind(exons1,
      data.frame(contig = fusions[fusion, "contig1"], type = "dummy",
                 start = breakpoint1 + 1000, end = breakpoint1 + 1000,
                 strand = exons1[1, "strand"], attributes = "",
                 geneID = exons1[1, "geneID"],
                 transcript = exons1[1, "transcript"],
                 exonNumber = "dummy",
                 left = breakpoint1 + 1000, right = breakpoint1 + 1000))
  }
  if (breakpoint2 < min(exons2$start)) {
    exons2 <- rbind(
      data.frame(contig = fusions[fusion, "contig2"], type = "dummy",
                 start = max(1, breakpoint2 - 1000),
                 end = max(1, breakpoint2 - 1000),
                 strand = exons2[1, "strand"], attributes = "",
                 geneID = exons2[1, "geneID"],
                 transcript = exons2[1, "transcript"],
                 exonNumber = "dummy",
                 left = max(1, breakpoint2 - 1000),
                 right = max(1, breakpoint2 - 1000)),
      exons2)
  } else if (breakpoint2 > max(exons2$end)) {
    exons2 <- rbind(exons2,
      data.frame(contig = fusions[fusion, "contig2"], type = "dummy",
                 start = breakpoint2 + 1000, end = breakpoint2 + 1000,
                 strand = exons2[1, "strand"], attributes = "",
                 geneID = exons2[1, "geneID"],
                 transcript = exons2[1, "transcript"],
                 exonNumber = "dummy",
                 left = breakpoint2 + 1000, right = breakpoint2 + 1000))
  }

  exons1$start <- as.integer(exons1$start)
  exons1$end <- as.integer(exons1$end)
  exons2$start <- as.integer(exons2$start)
  exons2$end <- as.integer(exons2$end)

  exons1$left <- exons1$start
  exons1$right <- exons1$end
  exons2$left <- exons2$start
  exons2$right <- exons2$end

  squashed_intron_size <- if (fixedScale > 0) fixedScale / 100 else 200

  if (squishIntrons) {
    result1 <- squishIntronsLogic(exons1, breakpoint1, squashed_intron_size)
    exons1 <- result1$exons
    breakpoint1 <- result1$breakpoint

    result2 <- squishIntronsLogic(exons2, breakpoint2, squashed_intron_size)
    exons2 <- result2$exons
    breakpoint2 <- result2$breakpoint
  } else {
    exons1$right <- exons1$right - min(exons1$left)
    breakpoint1 <- breakpoint1 - min(exons1$left)
    exons1$left <- exons1$left - min(exons1$left)
    exons2$right <- exons2$right - min(exons2$left)
    breakpoint2 <- breakpoint2 - min(exons2$left)
    exons2$left <- exons2$left - min(exons2$left)
  }

  scalingFactor <- max(exons1$right) + max(exons2$right)
  if (fixedScale > 0) {
    if (fixedScale >= scalingFactor) {
      scalingFactor <- fixedScale
    } else {
      warning(paste("fallback to automatic scaling, because value for",
                    "--fixedScale is too small (increase it to", scalingFactor,
                    "to avoid this)"))
    }
  }
  exons1$left <- exons1$left / scalingFactor
  exons1$right <- exons1$right / scalingFactor
  exons2$left <- exons2$left / scalingFactor
  exons2$right <- exons2$right / scalingFactor
  breakpoint1 <- breakpoint1 / scalingFactor
  breakpoint2 <- breakpoint2 / scalingFactor

  second_gene_horizontal_offset <- 1 + 0.05 - max(exons2$right)

  fusion_offset_gene1 <- (max(exons1$right) + second_gene_horizontal_offset) / 2 -
                   ifelse(fusions[fusion, "direction1"] == "downstream",
                         breakpoint1, max(exons1$right) - breakpoint1)
  fusion_offset_gene2 <- fusion_offset_gene1 +
                   ifelse(fusions[fusion, "direction1"] == "downstream",
                         breakpoint1, max(exons1$right) - breakpoint1)

  topRowPresent <- "fusion" %in% plotPanels
  bottomRowPresent <- any(c("circos", "domains", "readcounts") %in% plotPanels)
  layoutWidths <- c(
    ifelse("circos" %in% plotPanels, 1.1, 0.01),
    ifelse("domains" %in% plotPanels, 1.2, 0.01),
    ifelse("readcounts" %in% plotPanels, 0.7, 0.01)
  )
  layoutHeights <- c(
    ifelse(topRowPresent, 1.55, 0.3),
    ifelse(bottomRowPresent, 1.2, 0.01),
    ifelse("circos" %in% plotPanels, 0.25, 0.01)
  )
  layout(matrix(c(1, 1, 1, 2, 4, 5, 3, 4, 5), 3, 3, byrow = TRUE),
         widths = layoutWidths, heights = layoutHeights)
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "l", xlim = c(-0.12, 1.12),
       ylim = ifelse(rep(bottomRowPresent, 2), c(0.4, 1.1), c(0.2, 1.3)),
       bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

  y_sample_name <- ifelse(topRowPresent, config$layout$y_sample_name, 0.5)
  y_ideograms <- ifelse(alignmentsFile != "", config$layout$y_ideograms, 0.84)
  y_breakpoint_labels <- ifelse(alignmentsFile != "", config$layout$y_breakpoint_labels, 0.76)
  y_coverage <- config$layout$y_coverage
  y_exons <- config$layout$y_exons
  y_gene_names <- config$layout$y_gene_names
  y_fusion <- config$layout$y_fusion
  y_transcript <- config$layout$y_transcript
  y_scale <- config$layout$y_scale
  y_trajectory_breakpoint_labels <- y_breakpoint_labels - 0.035
  y_trajectory_exon_top <- y_exons + config$dimensions$exon_height
  y_trajectory_exon_bottom <- y_exons - 0.055
  y_trajectory_fusion <- y_fusion + config$dimensions$exon_height

  text(0.5, y_sample_name, sampleName, font = 2, cex = fontSize * 1.5,
       adj = c(0.5, 0))

  if ("fusion" %in% plotPanels) {

    if (!is.null(cytobands)) {
      drawIdeogram("left", min(exons1$left), max(exons1$right), y_ideograms,
                  cytobands, fusions[fusion, "contig1"],
                  fusions[fusion, "breakpoint1"])
      drawIdeogram("right", second_gene_horizontal_offset, second_gene_horizontal_offset + max(exons2$right),
                  y_ideograms, cytobands, fusions[fusion, "contig2"],
                  fusions[fusion, "breakpoint2"])
    }

    if (fusions[fusion, "gene1"] != ".")
      text(max(exons1$right) / 2, y_gene_names, fusions[fusion, "gene1"],
           font = 2, cex = fontSize, adj = c(0.5, 0))
    if (fusions[fusion, "site1"] != "intergenic")
      text(max(exons1$right) / 2, y_gene_names - 0.01, head(exons1$transcript, 1),
           cex = 0.9 * fontSize, adj = c(0.5, 1))
    if (fusions[fusion, "gene2"] != ".")
      text(second_gene_horizontal_offset + max(exons2$right) / 2, y_gene_names,
           fusions[fusion, "gene2"], font = 2, cex = fontSize, adj = c(0.5, 0))
    if (fusions[fusion, "site2"] != "intergenic")
      text(second_gene_horizontal_offset + max(exons2$right) / 2, y_gene_names - 0.01,
           head(exons2$transcript, 1), cex = 0.9 * fontSize, adj = c(0.5, 1))

    if (fusions[fusion, "site1"] == "intergenic")
      for (gene in unique(exons1$geneName)) {
        exonsOfGene <- exons1[exons1$geneName == gene & exons1$type != "dummy", ]
        if (any(exonsOfGene$type == "exon"))
          text(mean(c(min(exonsOfGene$left), max(exonsOfGene$right))),
               y_exons - 0.04, gene, cex = 0.9 * fontSize, adj = c(0.5, 1))
      }
    if (fusions[fusion, "site2"] == "intergenic")
      for (gene in unique(exons2$geneName)) {
        exonsOfGene <- exons2[exons2$geneName == gene & exons2$type != "dummy", ]
        if (any(exonsOfGene$type == "exon"))
          text(second_gene_horizontal_offset +
               mean(c(min(exonsOfGene$left), max(exonsOfGene$right))),
               y_exons - 0.04, gene, cex = 0.9 * fontSize, adj = c(0.5, 1))
      }

    text(breakpoint1 + 0.01, y_breakpoint_labels - 0.03,
         paste0("breakpoint1\n", fusions[fusion, "display_contig1"], ":",
               fusions[fusion, "breakpoint1"]),
         adj = c(1, 0), cex = fontSize)
    text(second_gene_horizontal_offset + breakpoint2 - 0.01, y_breakpoint_labels - 0.03,
         paste0("breakpoint2\n", fusions[fusion, "display_contig2"], ":",
               fusions[fusion, "breakpoint2"]),
         adj = c(0, 0), cex = fontSize)

    if (alignmentsFile != "") {
      lines(c(-0.02, -0.01, -0.01, -0.02),
            c(y_coverage, y_coverage, y_coverage + 0.1, y_coverage + 0.1))
      text(-0.025, y_coverage, "0", adj = c(1, 0.5), cex = 0.9 * fontSize)
      text(-0.025, y_coverage + 0.1, coverageNormalization1,
           adj = c(1, 0.5), cex = 0.9 * fontSize)
      text(-0.05, y_coverage + 0.08, "Coverage", srt = 90,
           cex = 0.9 * fontSize, adj = c(1, 0.5))

      if (length(coverageRange) == 2) {
        rightCoverageAxisX <- second_gene_horizontal_offset + max(exons2$right)
        lines(c(rightCoverageAxisX + 0.02, rightCoverageAxisX + 0.01,
                rightCoverageAxisX + 0.01, rightCoverageAxisX + 0.02),
              c(y_coverage, y_coverage, y_coverage + 0.1, y_coverage + 0.1))
        text(rightCoverageAxisX + 0.025, y_coverage, "0",
             adj = c(0, 0.5), cex = 0.9 * fontSize)
        text(rightCoverageAxisX + 0.025, y_coverage + 0.1, coverageNormalization2,
             adj = c(0, 0.5), cex = 0.9 * fontSize)
        text(rightCoverageAxisX + 0.05, y_coverage + 0.08, "Coverage",
             srt = 90, cex = 0.9 * fontSize, adj = c(1, 0.5))
      }

      rect(min(exons1$left), y_coverage, max(exons1$right), y_coverage + 0.1,
           col = "#eeeeee", border = NA)
      if (squishIntrons) {
        for (exon in seq_len(nrow(exons1)))
          if (exons1[exon, "type"] != "CDS")
            drawCoverage(exons1[exon, "left"], exons1[exon, "right"], y_coverage,
                        coverage1, exons1[exon, "start"], exons1[exon, "end"],
                        color1)
      } else {
        drawCoverage(min(exons1$left), max(exons1$right), y_coverage, coverage1,
                    min(exons1$start), max(exons1$start), color1)
      }

      rect(second_gene_horizontal_offset + min(exons2$left), y_coverage,
           second_gene_horizontal_offset + max(exons2$right), y_coverage + 0.1,
           col = "#eeeeee", border = NA)
      if (squishIntrons) {
        for (exon in seq_len(nrow(exons2)))
          if (exons2[exon, "type"] != "CDS")
            drawCoverage(second_gene_horizontal_offset + exons2[exon, "left"],
                        second_gene_horizontal_offset + exons2[exon, "right"], y_coverage,
                        coverage2, exons2[exon, "start"], exons2[exon, "end"],
                        color2)
      } else {
        drawCoverage(second_gene_horizontal_offset + min(exons2$left),
                    second_gene_horizontal_offset + max(exons2$right), y_coverage, coverage2,
                    min(exons2$start), max(exons2$start), color2)
      }
    }

    lines(c(min(exons1$left), max(exons1$right)), c(y_exons, y_exons),
          col = dark_color_gene1)
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene, ]
      drawStrand(min(exonsOfGene$left), max(exonsOfGene$left),
                y_exons, dark_color_gene1, head(exonsOfGene$strand, 1))
    }
    for (exon in seq_len(nrow(exons1)))
      drawExon(exons1[exon, "left"], exons1[exon, "right"], y_exons, color1,
              exons1[exon, "exonNumber"], exons1[exon, "type"], fontSize)

    lines(c(second_gene_horizontal_offset, second_gene_horizontal_offset + max(exons2$right)),
          c(y_exons, y_exons), col = dark_color_gene2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene, ]
      drawStrand(second_gene_horizontal_offset + min(exonsOfGene$left),
                second_gene_horizontal_offset + max(exonsOfGene$right),
                y_exons, dark_color_gene2, head(exonsOfGene$strand, 1))
    }
    for (exon in seq_len(nrow(exons2)))
      drawExon(second_gene_horizontal_offset + exons2[exon, "left"],
              second_gene_horizontal_offset + exons2[exon, "right"], y_exons, color2,
              exons2[exon, "exonNumber"], exons2[exon, "type"], fontSize)

    if (fusions[fusion, "direction1"] == "downstream") {
      lines(c(fusion_offset_gene1, fusion_offset_gene1 + breakpoint1),
            c(y_fusion, y_fusion), col = dark_color_gene1)
      for (gene in unique(exons1$geneName)) {
        exonsOfGene <- exons1[exons1$geneName == gene, ]
        if (min(exonsOfGene$start) <= fusions[fusion, "breakpoint1"])
          drawStrand(fusion_offset_gene1 + min(exonsOfGene$left),
                    fusion_offset_gene1 + min(breakpoint1, max(exonsOfGene$right)),
                    y_fusion, col = dark_color_gene1, exonsOfGene$strand[1])
      }
      for (exon in seq_len(nrow(exons1)))
        if (exons1[exon, "start"] <= fusions[fusion, "breakpoint1"])
          drawExon(fusion_offset_gene1 + exons1[exon, "left"],
                  fusion_offset_gene1 + min(breakpoint1, exons1[exon, "right"]),
                  y_fusion, color1, exons1[exon, "exonNumber"],
                  exons1[exon, "type"], fontSize)
      lines(c(0, 0, fusion_offset_gene1), c(y_trajectory_exon_top, y_trajectory_exon_bottom,
                                       y_trajectory_fusion), col = "red", lty = 2)
      lines(c(breakpoint1, breakpoint1, fusion_offset_gene1 + breakpoint1),
            c(y_trajectory_breakpoint_labels, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
    } else if (fusions[fusion, "direction1"] == "upstream") {
      lines(c(fusion_offset_gene1, fusion_offset_gene2), c(y_fusion, y_fusion), col = dark_color_gene1)
      for (gene in unique(exons1$geneName)) {
        exonsOfGene <- exons1[exons1$geneName == gene, ]
        if (max(exonsOfGene$end + 1) >= fusions[fusion, "breakpoint1"])
          drawStrand(fusion_offset_gene2 - max(exonsOfGene$right) + breakpoint1,
                    min(fusion_offset_gene2, fusion_offset_gene2 - min(exonsOfGene$left) + breakpoint1),
                    y_fusion, col = dark_color_gene1,
                    chartr("+-", "-+", exonsOfGene$strand[1]))
      }
      for (exon in seq_len(nrow(exons1)))
        if (exons1[exon, "end"] + 1 >= fusions[fusion, "breakpoint1"])
          drawExon(fusion_offset_gene1 + max(exons1$right) - exons1[exon, "right"],
                  min(fusion_offset_gene2,
                      fusion_offset_gene1 + max(exons1$right) - exons1[exon, "left"]),
                  y_fusion, color1, exons1[exon, "exonNumber"],
                  exons1[exon, "type"], fontSize)
      lines(c(max(exons1$right), max(exons1$right), fusion_offset_gene1),
            c(y_trajectory_exon_top, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
      lines(c(breakpoint1, breakpoint1,
              fusion_offset_gene1 + max(exons1$right) - breakpoint1),
            c(y_trajectory_breakpoint_labels, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
    }

    if (fusions[fusion, "direction2"] == "downstream") {
      lines(c(fusion_offset_gene2, fusion_offset_gene2 + breakpoint2),
            c(y_fusion, y_fusion), col = dark_color_gene2)
      for (gene in unique(exons2$geneName)) {
        exonsOfGene <- exons2[exons2$geneName == gene, ]
        if (min(exonsOfGene$start) <= fusions[fusion, "breakpoint2"])
          drawStrand(max(fusion_offset_gene2,
                        fusion_offset_gene2 + breakpoint2 - max(exonsOfGene$right)),
                    fusion_offset_gene2 + breakpoint2 - min(exonsOfGene$left),
                    y_fusion, col = dark_color_gene2,
                    chartr("+-", "-+", exonsOfGene$strand[1]))
      }
      for (exon in seq_len(nrow(exons2)))
        if (exons2[exon, "start"] <= fusions[fusion, "breakpoint2"])
          drawExon(max(fusion_offset_gene2,
                       fusion_offset_gene2 + breakpoint2 - exons2[exon, "right"]),
                  fusion_offset_gene2 + breakpoint2 - exons2[exon, "left"],
                  y_fusion, color2, exons2[exon, "exonNumber"],
                  exons2[exon, "type"], fontSize)
      lines(c(second_gene_horizontal_offset, second_gene_horizontal_offset, fusion_offset_gene2 + breakpoint2),
            c(y_trajectory_exon_top, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
      lines(c(second_gene_horizontal_offset + breakpoint2, second_gene_horizontal_offset + breakpoint2,
              fusion_offset_gene2),
            c(y_trajectory_breakpoint_labels, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
    } else if (fusions[fusion, "direction2"] == "upstream") {
      lines(c(fusion_offset_gene2, fusion_offset_gene2 + max(exons2$right) - breakpoint2),
            c(y_fusion, y_fusion), col = dark_color_gene2)
      for (gene in unique(exons2$geneName)) {
        exonsOfGene <- exons2[exons2$geneName == gene, ]
        if (max(exonsOfGene$end + 1) >= fusions[fusion, "breakpoint2"])
          drawStrand(max(fusion_offset_gene2,
                        fusion_offset_gene2 + min(exonsOfGene$left) - breakpoint2),
                    fusion_offset_gene2 + max(exonsOfGene$right) - breakpoint2,
                    y_fusion, col = dark_color_gene2, exonsOfGene$strand[1])
      }
      for (exon in seq_len(nrow(exons2)))
        if (exons2[exon, "end"] + 1 >= fusions[fusion, "breakpoint2"])
          drawExon(max(fusion_offset_gene2,
                       fusion_offset_gene2 + exons2[exon, "left"] - breakpoint2),
                  fusion_offset_gene2 + exons2[exon, "right"] - breakpoint2,
                  y_fusion, color2, exons2[exon, "exonNumber"],
                  exons2[exon, "type"], fontSize)
      lines(c(second_gene_horizontal_offset + max(exons2$right),
              second_gene_horizontal_offset + max(exons2$right),
              fusion_offset_gene2 + max(exons2$right) - breakpoint2),
            c(y_trajectory_exon_top, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
      lines(c(second_gene_horizontal_offset + breakpoint2, second_gene_horizontal_offset + breakpoint2,
              fusion_offset_gene2),
            c(y_trajectory_breakpoint_labels, y_trajectory_exon_bottom, y_trajectory_fusion),
            col = "red", lty = 2)
    }

    if (fusions[fusion, "fusion_transcript"] != ".") {
      fusion_transcript1 <- gsub("\\|.*", "", fusions[fusion, "fusion_transcript"],
                                 perl = TRUE)
      fusion_transcript1 <- substr(fusion_transcript1,
                                   max(1, nchar(fusion_transcript1) - 30),
                                   nchar(fusion_transcript1))
      fusion_transcript2 <- gsub(".*\\|", "", fusions[fusion, "fusion_transcript"],
                                 perl = TRUE)
      fusion_transcript2 <- substr(fusion_transcript2, 1,
                                   min(nchar(fusion_transcript2), 30))
      non_template_bases <- gsub(".*\\|([^|]*)\\|.*", "\\1",
                                  fusions[fusion, "fusion_transcript"], perl = TRUE)
      if (non_template_bases == fusions[fusion, "fusion_transcript"])
        non_template_bases <- ""
      non_template_bases1 <- substr(non_template_bases, 1,
                                    floor(nchar(non_template_bases) / 2))
      non_template_bases2 <- substr(non_template_bases,
                                    ceiling(nchar(non_template_bases) / 2 + 0.5),
                                    nchar(non_template_bases))
      text(fusion_offset_gene2, y_transcript,
           bquote(.(fusion_transcript1) * phantom(.(non_template_bases1))),
           col = dark_color_gene1, adj = c(1, 0.5), cex = fontSize)
      text(fusion_offset_gene2, y_transcript,
           bquote(phantom(.(non_template_bases2)) * .(fusion_transcript2)),
           col = dark_color_gene2, adj = c(0, 0.5), cex = fontSize)
      text(fusion_offset_gene2, y_transcript, non_template_bases1,
           adj = c(1, 0.5), cex = fontSize)
      text(fusion_offset_gene2, y_transcript, non_template_bases2,
           adj = c(0, 0.5), cex = fontSize)
    }

    realScale <- max(exons1$end - exons1$start, exons2$end - exons2$start)
    mapScale <- max(exons1$right - exons1$left, exons2$right - exons2$left)
    desired_scale_size <- config$layout$desired_scale_size
    whisker_size <- config$layout$scale_whisker_size
    realScale <- desired_scale_size / mapScale * realScale
    mapScale <- desired_scale_size
    realScaleOptimalFit <- signif(realScale, 1)
    mapScaleOptimalFit <- realScaleOptimalFit / realScale * mapScale
    lines(c(1 - mapScaleOptimalFit, 1), c(y_scale, y_scale))
    lines(c(1 - mapScaleOptimalFit, 1 - mapScaleOptimalFit),
          c(y_scale - whisker_size, y_scale + whisker_size))
    lines(c(1, 1), c(y_scale - whisker_size, y_scale + whisker_size))
    realScaleThousands <- max(0, min(3, floor(log10(realScaleOptimalFit) / 3)))
    scaleUnits <- c("bp", "kbp", "Mbp", "Gbp")
    scaleLabel <- paste(realScaleOptimalFit / max(1, 1000 ^ realScaleThousands),
                       scaleUnits[realScaleThousands + 1])
    text(1 - mapScaleOptimalFit / 2, y_scale + 0.005, scaleLabel,
         adj = c(0.5, 0), cex = fontSize * 0.9)
    if (squishIntrons)
      text(1 - mapScaleOptimalFit / 2, y_scale - 0.005, "introns not to scale",
           adj = c(0.5, 1), cex = fontSize * 0.9, font = 3)

  }

  if (!("circos" %in% plotPanels)) {
    plot(0, 0, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    plot(0, 0, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  } else {
    par(mar = c(2, 4, 0, 0), xpd = NA)
    circos.clear()
    circos.initializeWithIdeogram(cytoband = cytobands, plotType = NULL)
    geneLabels <- data.frame(
      contig = c(fusions[fusion, "contig1"], fusions[fusion, "contig2"]),
      start = c(fusions[fusion, "breakpoint1"], fusions[fusion, "breakpoint2"])
    )
    # 确保 contig 名称与 cytobands 一致（去除 chr 前缀以匹配）
    geneLabels$contig <- removeChr(geneLabels$contig)
    geneLabels$end <- geneLabels$start + 1
    geneLabels$gene <- c(fusions[fusion, "gene1"], fusions[fusion, "gene2"])
    geneLabels$gene <- ifelse(
      c(fusions[fusion, "site1"], fusions[fusion, "site2"]) == "intergenic",
      paste0(c(fusions[fusion, "display_contig1"],
               fusions[fusion, "display_contig2"]), ":",
             geneLabels$start),
      geneLabels$gene)
    circos.genomicLabels(geneLabels, labels.column = 4, side = "outside",
                         cex = fontSize, labels_height = 0.27)
    for (contig in unique(cytobands$contig)) {
      if (contig %in% get.all.sector.index()) {
        set.current.cell(track.index = 2, sector.index = contig)
        circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex = 0.85)
      }
    }
    circos.genomicIdeogram(cytobands)
    confidenceRank <- c(low = 0, medium = 1, high = 2)
    for (i in c(setdiff(seq_len(nrow(fusions)), fusion), fusion)) {
      f <- fusions[i, ]
      # 确保 contig 名称一致性后再比较
      f_contig1 <- removeChr(f$contig1)
      f_contig2 <- removeChr(f$contig2)
      if (any(cytobands$contig == f_contig1) &&
          any(cytobands$contig == f_contig2))
        if (minConfidenceForCircosPlot != "none" &&
            confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot] ||
            i == fusion)
          circos.link(f_contig1, f$breakpoint1, f_contig2, f$breakpoint2,
                     lwd = 2,
                     col = ifelse(i == fusion, circosColors[f$type],
                                getBrightColor(circosColors[f$type])))
    }
    plot(0, 0, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")
    legend("top", legend = names(circosColors),
           col = sapply(circosColors, getBrightColor), lwd = 3,
           ncol = 2, box.lty = 0)
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
  }

  plot(0, 0, type = "l", xlim = c(-0.1, 1.1), ylim = c(0, 1),
       bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  par(xpd = NA)
  if ("domains" %in% plotPanels) {

    exon_height <- config$domain$exon_height
    y_exons_domain <- 0.5
    y_gene_names_domain <- y_exons_domain - exon_height / 2 - config$domain$gene_names_y_offset

    coding_exons_gene1 <- exons1[exons1$type == "CDS" &
                           fusions[fusion, "site1"] != "intergenic", ]
    coding_exons_gene2 <- exons2[exons2$type == "CDS" &
                           fusions[fusion, "site2"] != "intergenic", ]

    if (fusions[fusion, "direction1"] == "upstream") {
      coding_exons_gene1 <- coding_exons_gene1[coding_exons_gene1$end >= fusions[fusion, "breakpoint1"], ]
      coding_exons_gene1$start <- ifelse(coding_exons_gene1$start < fusions[fusion, "breakpoint1"],
                                  fusions[fusion, "breakpoint1"],
                                  coding_exons_gene1$start)
    } else {
      coding_exons_gene1 <- coding_exons_gene1[coding_exons_gene1$start <= fusions[fusion, "breakpoint1"], ]
      coding_exons_gene1$end <- ifelse(coding_exons_gene1$end > fusions[fusion, "breakpoint1"],
                                fusions[fusion, "breakpoint1"],
                                coding_exons_gene1$end)
    }
    if (fusions[fusion, "direction2"] == "upstream") {
      coding_exons_gene2 <- coding_exons_gene2[coding_exons_gene2$end >= fusions[fusion, "breakpoint2"], ]
      coding_exons_gene2$start <- ifelse(coding_exons_gene2$start < fusions[fusion, "breakpoint2"],
                                  fusions[fusion, "breakpoint2"],
                                  coding_exons_gene2$start)
    } else {
      coding_exons_gene2 <- coding_exons_gene2[coding_exons_gene2$start <= fusions[fusion, "breakpoint2"], ]
      coding_exons_gene2$end <- ifelse(coding_exons_gene2$end > fusions[fusion, "breakpoint2"],
                                fusions[fusion, "breakpoint2"],
                                coding_exons_gene2$end)
    }

    exonsGRanges1 <- GRanges(coding_exons_gene1$contig,
                             IRanges(coding_exons_gene1$start, coding_exons_gene1$end),
                             strand = coding_exons_gene1$strand)
    exonsGRanges2 <- GRanges(coding_exons_gene2$contig,
                             IRanges(coding_exons_gene2$start, coding_exons_gene2$end),
                             strand = coding_exons_gene2$strand)
    domainsGRanges <- GRanges(proteinDomains$contig,
                             IRanges(proteinDomains$start, proteinDomains$end),
                             strand = proteinDomains$strand)
    domainsGRanges$proteinDomainName <- proteinDomains$proteinDomainName
    domainsGRanges$proteinDomainID <- proteinDomains$proteinDomainID
    domainsGRanges$color <- proteinDomains$color
    domainsGRanges <- domainsGRanges[
      suppressWarnings(unique(queryHits(
        findOverlaps(domainsGRanges, union(exonsGRanges1, exonsGRanges2))
      )))
    ]

    domainsGRangesList <- GRangesList(lapply(
      unique(domainsGRanges$proteinDomainID),
      function(x) domainsGRanges[domainsGRanges$proteinDomainID == x]
    ))

    trimDomains <- function(domainsGRangesList, exonsGRanges) {
      do.call(rbind, lapply(domainsGRangesList, function(x) {
        intersected <- as.data.frame(reduce(
          suppressWarnings(intersect(x, exonsGRanges))
        ))
        if (nrow(intersected) > 0) {
          intersected$proteinDomainName <- head(x$proteinDomainName, 1)
          intersected$proteinDomainID <- head(x$proteinDomainID, 1)
          intersected$color <- head(x$color, 1)
        } else {
          intersected$proteinDomainName <- character()
          intersected$proteinDomainID <- character()
          intersected$color <- character()
        }
        return(intersected)
      }))
    }
    coding_domains1 <- trimDomains(domainsGRangesList, exonsGRanges1)
    coding_domains2 <- trimDomains(domainsGRangesList, exonsGRanges2)

    coding_exons_gene1$length <- coding_exons_gene1$end - coding_exons_gene1$start + 1
    coding_exons_gene2$length <- coding_exons_gene2$end - coding_exons_gene2$start + 1

    if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
      text(0.5, 0.5, "Genes are not protein-coding.")
    } else {
      codingLength1 <- sum(coding_exons_gene1$length)
      codingLength2 <- sum(coding_exons_gene2$length)
      if (codingLength1 + codingLength2 > 0) {
        if (!((codingLength1 == 0 ||
               grepl("\\.$", fusions[fusion, "strand1"])) &&
              (codingLength2 == 0 ||
               grepl("\\.$", fusions[fusion, "strand2"])))) {

          antisenseTranscription1 <- sub("/.*", "", fusions[fusion, "strand1"]) !=
            sub(".*/", "", fusions[fusion, "strand1"])
          antisenseTranscription2 <- sub("/.*", "", fusions[fusion, "strand2"]) !=
            sub(".*/", "", fusions[fusion, "strand2"])
          if (!((codingLength1 == 0 || antisenseTranscription1) &&
                (codingLength2 == 0 || antisenseTranscription2))) {

            coding_domains1 <- removeIntronsFromProteinDomains(coding_exons_gene1,
                                                               coding_domains1)
            coding_domains2 <- removeIntronsFromProteinDomains(coding_exons_gene2,
                                                               coding_domains2)

            if (!is.null(coding_domains1) || !is.null(coding_domains2)) {

              mergeSimilarDomains <- function(domains, mergeDomainsOverlappingBy) {
                if (is.null(domains)) return(domains)
                merged <- domains[FALSE, ]
                domains <- domains[order(domains$end - domains$start,
                                         decreasing = TRUE), ]
                for (domain in rownames(domains)) {
                  if (!any((abs(merged$start - domains[domain, "start"]) +
                            abs(merged$end - domains[domain, "end"])) /
                           (domains[domain, "end"] - domains[domain, "start"] + 1) <=
                           1 - mergeDomainsOverlappingBy))
                    merged <- rbind(merged, domains[domain, ])
                }
                return(merged)
              }
              coding_domains1 <- mergeSimilarDomains(coding_domains1,
                                                      mergeDomainsOverlappingBy)
              coding_domains2 <- mergeSimilarDomains(coding_domains2,
                                                      mergeDomainsOverlappingBy)

              if (optimizeDomainColors) {
                uniqueDomains <- unique(c(coding_domains1$proteinDomainID,
                                           coding_domains2$proteinDomainID))
                colors <- rainbow(length(uniqueDomains))
                colors <- apply(col2rgb(colors), 2,
                               function(x) 0.3 + y / 255 * 0.7)
                colors <- apply(colors, 2,
                               function(x) rgb(x["red"], x["green"], x["blue"]))
                names(colors) <- uniqueDomains
                coding_domains1$color <- colors[coding_domains1$proteinDomainID]
                coding_domains2$color <- colors[coding_domains2$proteinDomainID]
              }

              if (any(coding_exons_gene1$strand == "-")) {
                coding_exons_gene1$length <- rev(coding_exons_gene1$length)
                temp <- coding_domains1$end
                coding_domains1$end <- codingLength1 - coding_domains1$start
                coding_domains1$start <- codingLength1 - temp
              }
              if (any(coding_exons_gene2$strand == "-")) {
                coding_exons_gene2$length <- rev(coding_exons_gene2$length)
                temp <- coding_domains2$end
                coding_domains2$end <- codingLength2 - coding_domains2$start
                coding_domains2$start <- codingLength2 - temp
              }

              coding_exons_gene1$length <- coding_exons_gene1$length /
                                     (codingLength1 + codingLength2)
              coding_exons_gene2$length <- coding_exons_gene2$length /
                                     (codingLength1 + codingLength2)
              coding_domains1$start <- coding_domains1$start /
                                        (codingLength1 + codingLength2)
              coding_domains1$end <- coding_domains1$end /
                                        (codingLength1 + codingLength2)
              coding_domains2$start <- coding_domains2$start /
                                        (codingLength1 + codingLength2)
              coding_domains2$end <- coding_domains2$end /
                                        (codingLength1 + codingLength2)

              rect(0, y_exons_domain - exon_height / 2,
                   sum(coding_exons_gene1$length), y_exons_domain + exon_height / 2,
                   col = color1, border = NA)
              rect(sum(coding_exons_gene1$length), y_exons_domain - exon_height / 2,
                   sum(coding_exons_gene1$length) + sum(coding_exons_gene2$length),
                   y_exons_domain + exon_height / 2, col = color2, border = NA)

              exonBoundaries <- cumsum(c(coding_exons_gene1$length,
                                         coding_exons_gene2$length))
              if (length(exonBoundaries) > 1) {
                exonBoundaries <- exonBoundaries[-length(exonBoundaries)]
                for (exonBoundary in exonBoundaries)
                  lines(c(exonBoundary, exonBoundary),
                        c(y_exons_domain - exon_height, y_exons_domain + exon_height),
                        col = "white", lty = 3)
              }

              nestDomains <- function(domains) {
                if (length(unlist(domains)) == 0) return(domains)
                domains <- domains[order(domains$end - domains$start,
                                         decreasing = TRUE), ]
                rownames(domains) <- seq_len(nrow(domains))
                domains$parent <- 0
                for (domain in rownames(domains))
                  domains[domains$start >= domains[domain, "start"] &
                          domains$end <= domains[domain, "end"] &
                          rownames(domains) != domain, "parent"] <- domain
                maxOverlappingDomains <- max(1, as.integer(coverage(
                  IRanges(domains$start * 10e6, domains$end * 10e6))))
                padding <- 1 / maxOverlappingDomains * 0.4
                domains$y <- 0
                domains$height <- 0
                adjustPositionAndHeight <- function(parentDomain, y, height,
                                                    padding, e) {
                  for (domain in which(e$domains$parent == parentDomain)) {
                    overlappingDomains <- which(
                      (between(e$domains$start,
                               e$domains[domain, "start"],
                               e$domains[domain, "end"]) |
                       between(e$domains$end, e$domains[domain, "start"],
                              e$domains[domain, "end"])) &
                       e$domains$parent == parentDomain)
                    e$domains[domain, "height"] <-
                      height / length(overlappingDomains) -
                      padding * (length(overlappingDomains) - 1) /
                      length(overlappingDomains)
                    e$domains[domain, "y"] <-
                      y + (which(domain == overlappingDomains) - 1) *
                      (e$domains[domain, "height"] + padding)
                    adjustPositionAndHeight(domain,
                                            e$domains[domain, "y"] + padding,
                                            e$domains[domain, "height"] -
                                            2 * padding, padding, e)
                  }
                }
                adjustPositionAndHeight(0, 0, 1, padding, environment())
                domains <- domains[order(domains$height, decreasing = TRUE), ]
                return(domains)
              }
              coding_domains1 <- nestDomains(coding_domains1)
              coding_domains2 <- nestDomains(coding_domains2)
              coding_domains1$y <-
                y_exons_domain - exon_height / 2 + 0.025 +
                (exon_height - 2 * 0.025) * coding_domains1$y
              coding_domains2$y <-
                y_exons_domain - exon_height / 2 + 0.025 +
                (exon_height - 2 * 0.025) * coding_domains2$y
              coding_domains1$height <-
                coding_domains1$height * (exon_height - 2 * 0.025)
              coding_domains2$height <-
                coding_domains2$height * (exon_height - 2 * 0.025)

              drawProteinDomainRect <- function(left, bottom, right, top, color) {
                rect(left, bottom, right, top, col = color,
                     border = getDarkColor(color))
                gradientSteps <- 20
                drawVerticalGradient(rep(left, gradientSteps),
                                    rep(right, gradientSteps),
                                    seq(top, bottom, length.out = gradientSteps),
                                    rgb(1, 1, 1, 0.7))
                drawVerticalGradient(rep(left, gradientSteps),
                                    rep(right, gradientSteps),
                                    seq(bottom, bottom + (top - bottom) * 0.4,
                                        length.out = gradientSteps),
                                    rgb(0, 0, 0, 0.1))
              }
              if (length(unlist(coding_domains1)) > 0)
                for (domain in seq_len(nrow(coding_domains1)))
                  drawProteinDomainRect(
                    coding_domains1[domain, "start"],
                    coding_domains1[domain, "y"],
                    coding_domains1[domain, "end"],
                    coding_domains1[domain, "y"] +
                    coding_domains1[domain, "height"],
                    coding_domains1[domain, "color"])
              if (length(unlist(coding_domains2)) > 0)
                for (domain in seq_len(nrow(coding_domains2)))
                  drawProteinDomainRect(
                    sum(coding_exons_gene1$length) +
                    coding_domains2[domain, "start"],
                    coding_domains2[domain, "y"],
                    sum(coding_exons_gene1$length) +
                    coding_domains2[domain, "end"],
                    coding_domains2[domain, "y"] +
                    coding_domains2[domain, "height"],
                    coding_domains2[domain, "color"])

              if (codingLength1 > 0)
                text(sum(coding_exons_gene1$length) / 2, y_gene_names_domain,
                     fusions[fusion, "gene1"], font = 2, cex = fontSize)
              if (codingLength2 > 0)
                text(sum(coding_exons_gene1$length) +
                     sum(coding_exons_gene2$length) / 2, y_gene_names_domain,
                     fusions[fusion, "gene2"], font = 2, cex = fontSize)

              countUniqueDomains <- function(domains) {
                uniqueDomains <- 0
                if (length(unlist(domains)) > 0) {
                  uniqueDomains <- 1
                  if (nrow(domains) > 1) {
                    previousDomain <- domains[1, "proteinDomainID"]
                    for (domain in 2:nrow(domains)) {
                      if (previousDomain != domains[domain, "proteinDomainID"])
                        uniqueDomains <- uniqueDomains + 1
                      previousDomain <- domains[domain, "proteinDomainID"]
                    }
                  }
                }
                return(uniqueDomains)
              }
              if (length(unlist(coding_domains1)) > 0)
                coding_domains1 <- coding_domains1[
                  order(coding_domains1$start), ]
              uniqueDomains1 <- countUniqueDomains(coding_domains1)
              if (length(unlist(coding_domains2)) > 0)
                coding_domains2 <- coding_domains2[
                  order(coding_domains2$end, decreasing = TRUE), ]
              uniqueDomains2 <- countUniqueDomains(coding_domains2)

              titleY <- y_exons_domain + exon_height / 2 + (uniqueDomains1 + 2) * 0.05
              text(0.5, titleY + 0.01, "RETAINED PROTEIN DOMAINS",
                   adj = c(0.5, 0), font = 2, cex = fontSize)
              text(0.5, titleY,
                   ifelse(fusions[fusion, "reading_frame"] %in%
                          c("in-frame", "out-of-frame"),
                         paste(fusions[fusion, "reading_frame"], "fusion"),
                         ifelse(fusions[fusion, "reading_frame"] ==
                                "stop-codon",
                                "stop codon before fusion junction",
                                "reading frame unclear")),
                   adj = c(0.5, 1), cex = fontSize)

              if (length(unlist(coding_domains1)) > 0) {
                previousConnectorX <- -1
                previousLabelX <- -1
                labelY <- y_exons_domain + exon_height / 2 + uniqueDomains1 * 0.05
                for (domain in seq_len(nrow(coding_domains1))) {
                  connectorX <- min(
                    coding_domains1[domain, "start"] + 0.01,
                    (coding_domains1[domain, "start"] +
                     coding_domains1[domain, "end"]) / 2)
                  if (connectorX - previousConnectorX < 0.01 &&
                      coding_domains1[domain, "end"] >
                      previousConnectorX + 0.01)
                    connectorX <- previousConnectorX + 0.01
                  labelX <- max(connectorX, previousLabelX) + 0.02
                  adjacentDomainsOfSameType <-
                    domain + 1 <= nrow(coding_domains1) &&
                    coding_domains1[domain + 1, "proteinDomainID"] ==
                    coding_domains1[domain, "proteinDomainID"]
                  if (adjacentDomainsOfSameType) {
                    labelX <- coding_domains1[domain + 1, "start"] + 0.015
                  } else {
                    text(labelX, labelY,
                         coding_domains1[domain, "proteinDomainName"],
                         adj = c(0, 0.5),
                         col = getDarkColor(coding_domains1[domain, "color"]),
                         cex = fontSize)
                  }
                  lines(c(labelX - 0.005, connectorX, connectorX),
                        c(labelY, labelY,
                          coding_domains1[domain, "y"] +
                          coding_domains1[domain, "height"]),
                        col = getDarkColor(coding_domains1[domain, "color"]))
                  if (!adjacentDomainsOfSameType)
                    labelY <- labelY - 0.05
                  previousConnectorX <- connectorX
                  previousLabelX <- labelX
                }
              }

              if (length(unlist(coding_domains2)) > 0) {
                previousConnectorX <- 100
                previousLabelX <- 100
                labelY <- y_exons_domain - exon_height / 2 - (uniqueDomains2 + 1) * 0.05
                for (domain in seq_len(nrow(coding_domains2))) {
                  connectorX <- sum(coding_exons_gene1$length) +
                    max(coding_domains2[domain, "end"] - 0.01,
                        (coding_domains2[domain, "start"] +
                         coding_domains2[domain, "end"]) / 2)
                  if (previousConnectorX - connectorX < 0.01 &&
                      sum(coding_exons_gene1$length) +
                      coding_domains2[domain, "start"] <
                      previousConnectorX - 0.01)
                    connectorX <- previousConnectorX - 0.01
                  labelX <- min(connectorX, previousLabelX) - 0.02
                  adjacentDomainsOfSameType <-
                    domain + 1 <= nrow(coding_domains2) &&
                    coding_domains2[domain + 1, "proteinDomainID"] ==
                    coding_domains2[domain, "proteinDomainID"]
                  if (adjacentDomainsOfSameType) {
                    labelX <- sum(coding_exons_gene1$length) +
                      coding_domains2[domain + 1, "end"] - 0.015
                  } else {
                    text(labelX, labelY,
                         coding_domains2[domain, "proteinDomainName"],
                         adj = c(1, 0.5),
                         col = getDarkColor(coding_domains2[domain, "color"]),
                         cex = fontSize)
                  }
                  lines(c(labelX + 0.005, connectorX, connectorX),
                        c(labelY, labelY, coding_domains2[domain, "y"]),
                        col = getDarkColor(coding_domains2[domain, "color"]))
                  if (!adjacentDomainsOfSameType)
                    labelY <- labelY + 0.05
                  previousConnectorX <- connectorX
                  previousLabelX <- labelX
                }
              }
            }
          }
        }
      }
    }
  }
  par(xpd = FALSE)

  if ("readcounts" %in% plotPanels) {
    plot(0, 0, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(0, 0.575, "SUPPORTING READ COUNT", font = 2,
         adj = c(0, 0), cex = fontSize)
    if ("split_reads" %in% colnames(fusions)) {
      text(0, 0.525,
           paste0("Split reads = ", fusions[fusion, "split_reads"], "\n",
                 "Discordant mates = ", fusions[fusion, "discordant_mates"]),
           adj = c(0, 1), cex = fontSize)
    } else {
      text(0, 0.525,
           paste0("Split reads at breakpoint1 = ",
                 fusions[fusion, "split_reads1"], "\n",
                 "Split reads at breakpoint2 = ",
                 fusions[fusion, "split_reads2"], "\n",
                 "Discordant mates = ",
                 fusions[fusion, "discordant_mates"]),
           adj = c(0, 1), cex = fontSize)
    }
  }
}

devNull <- dev.off()
message("Done")
