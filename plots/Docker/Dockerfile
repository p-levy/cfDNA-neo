# Use the official rocker tidyverse image (includes R + tidyverse)
FROM rocker/tidyverse:4.3.2

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update system packages and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libncurses5-dev \
    libncursesw5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and remaining CRAN packages (tidyverse already included in base image)
RUN R -e "install.packages(c('BiocManager', 'data.table', 'ggridges', 'ggrepel', 'patchwork', 'scales'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'circlize', 'GenomicRanges', 'EnrichedHeatmap'))"

# Default command
CMD ["/bin/bash"]