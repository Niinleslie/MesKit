FROM debian:testing
#ENV DEBIAN_FRONTEND noninteractive

# Set the 
RUN apt-get clean
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y locales && \
    apt-get install -y libterm-readkey-perl

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8 

RUN apt-get install --assume-yes apt-utils &&\
    apt-get install -y software-properties-common &&\
    apt-get install -y gnupg2 

## Now install R and littler, and create a link for littler in /usr/local/bin
ENV R_BASE_VERSION 4.0.0

RUN apt-get install -y --no-install-recommends \
        #gcc-9-base \
        libopenblas0-pthread \
		littler \
        r-cran-littler \
		r-base=${R_BASE_VERSION}-* \
		r-base-dev=${R_BASE_VERSION}-* \
		r-recommended=${R_BASE_VERSION}-* \
	&& ln -s /usr/lib/R/site-library/littler/examples/build.r /usr/local/bin/build.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/check.r /usr/local/bin/check.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get -y update && apt-get install -y -f --no-install-recommends \
    apt-utils \
    software-properties-common \
    wget \
    sudo \
    libgfortran5 \
    gdebi-core \
    libgfortran5 \
    libcurl4-gnutls-dev \
    libxt-dev \
    libssl-dev \
    libxml2 \
    libxml2-dev \
    libpng-dev

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
#RUN apt-get clean && apt-get -y update && apt-get install -y locales && locale-gen en_US.UTF-8
#ENV LANG='en_US.UTF-8' LANGUAGE='en_US.UTF-8' LC_ALL='en_US.UTF-8'


# Download and install shiny server
RUN R -e "install.packages(c('DT', 'devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyWidgets', 'shinyBS', 'shinycssloaders', 'shinyjs', 'rjson'), repos='http://cran.rstudio.com/')"
RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb \
    && gdebi -n shiny-server-1.5.9.923-amd64.deb \
    && rm -f shiny-server-1.5.9.923-amd64.deb

# MesKit part:
RUN R -e "BiocManager::install(c('BSgenome', 'GenomeInfoDb', 'org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19'))" 
RUN R -e "devtools::install_github('Niinleslie/MesKit', ref = 'master')"

# shiny server application & configuration
#COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY inst/shiny /srv/shiny-server/

EXPOSE 3838

# Copy further configuration files into the Docker image
#COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/bin/bash"]
#CMD ["/usr/bin/shiny-server.sh"]
