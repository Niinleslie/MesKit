FROM debian:latest
COPY sources.list   /etc/apt/sources.list

#FROM debian:latest
RUN apt-get clean
RUN apt-get -y update && apt-get install -y -f --no-install-recommends \
    build-essential \
    wget \
    sudo \
    libgfortran5 \
    gdebi-core \
    # r-needed:
    r-base \
    # for devtools: https://stackoverflow.com/questions/31114991/installation-of-package-devtools-had-non-zero-exit-status-in-a-powerpc
    apt-utils \
    libgfortran5 \
    libcurl4-gnutls-dev \
    libxt-dev \
    libssl-dev \
    libxml2 \
    libxml2-dev \
    libpng-dev

# Download and install shiny server
#RUN R -e "options(repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN')"
RUN sudo R -e "install.packages(c('shiny', 'devtools', 'BiocManager'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN')"
RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb \
    && gdebi -n shiny-server-1.5.9.923-amd64.deb \
    && rm -f shiny-server-1.5.9.923-amd64.deb

# MesKit part:
RUN R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor')"
RUN R -e "remotes::install_github('Niinleslie/Meskit')"
RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))"

# shiny server application & configuration
COPY inst/shiny /srv/shiny-server/

EXPOSE 8888

CMD ["/bin/bash"]
