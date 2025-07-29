FROM ubuntu:25.04

RUN apt-get -y update \
      && apt-get -y install r-base

RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('deepSNV'))"
RUN R -e "install.packages(c('tidyverse','foreach','doParallel'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
