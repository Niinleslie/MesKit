FROM debian:latest
COPY sources.list   /etc/apt/sources.list

#FROM debian:latest
RUN apt-get clean
RUN apt-get -y update && apt-get install -y -f --no-install-recommends \
    #build-essential \
    #gfortran \
    wget \
    sudo \
    libgfortran5 \
    gdebi-core \
    # r-needed:
    r-base \
    r-base-dev \
    libgfortran5 \
    libcurl4-gnutls-dev \
    libxt-dev \
    libssl-dev \
    libxml2 \
    libxml2-dev \
    libpng-dev

# Download and install shiny server
#RUN R -e "options(repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN')"
RUN R -e "install.packages(c('DT', 'devtools', 'BiocManager'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN')"
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyWidgets', 'shinyBS', 'shinycssloaders', 'shinyjs'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN')"
RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb \
    && gdebi -n shiny-server-1.5.9.923-amd64.deb \
    && rm -f shiny-server-1.5.9.923-amd64.deb


# MesKit part:
#RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))"
RUN R -e "BiocManager::install(c('BSgenome', 'GenomeInfoDb', 'org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))" 

#RUN R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor')"
RUN R -e "devtools::install_github('Ninomoriaty/MesKit', ref = 'master')"


# shiny server application & configuration
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY inst/shiny /srv/shiny-server/

EXPOSE 3838

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

#CMD ["/bin/bash"]
CMD ["/usr/bin/shiny-server.sh"]
