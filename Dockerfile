# draw_fusions.R Dockerfile
# 基于 R 官方镜像 (Debian slim 版本，体积较小)

FROM r-base:4.3.1

# 设置工作目录
WORKDIR /app

# 安装系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libcairo2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg62-turbo-dev \
    && rm -rf /var/lib/apt/lists/*

# 安装 R 依赖包
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'GenomicAlignments', 'circlize'), ask=FALSE, update=TRUE)"

# 复制 R 脚本
COPY draw_fusions.R /app/draw_fusions.R

# 设置执行权限
RUN chmod +x /app/draw_fusions.R

# 设置默认参数（如果需要可以覆盖）
ENV SAMPLE_NAME="fusion_sample"
ENV PLOT_PANELS="fusion,circos,domains,readcounts"

# 入口点
ENTRYPOINT ["Rscript", "draw_fusions.R"]

# 默认命令（帮助信息）
CMD ["--help"]
