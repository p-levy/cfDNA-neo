# Use the official rocker tidyverse image (includes R + tidyverse)
FROM rocker/tidyverse:4.3.2

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update system packages and install dependencies for additional packages
RUN apt-get update && apt-get install -y \
    bedtools \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install devtools for GitHub installations
RUN R -e "install.packages('devtools', repos='https://cran.rstudio.com/')"

# Install CRAN packages
RUN R -e "install.packages(c('data.table', 'argparse'), repos='https://cran.rstudio.com/')"

# Install bedtoolsr from GitHub
RUN R -e "devtools::install_github('PhanstielLab/bedtoolsr')"

# Default command
CMD ["/bin/bash"]