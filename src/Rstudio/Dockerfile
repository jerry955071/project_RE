FROM rocker/rstudio:4.2.2

# install packages
RUN Rscript -e "install.packages(c('foreach', 'doParallel', 'magrittr', 'rjson', 'data.table'));"

# install tabix
RUN apt-get update && \
    apt-get install -y \
        libncurses5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        tabix && \
    apt-get clean
