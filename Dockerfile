FROM debian:testing
ENV DEBIAN_FRONTEND noninteractive


RUN apt-get autoclean && \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends locales libterm-readkey-perl 
    #rm -rf /var/lib/apt/lists/*

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8 

RUN apt-get update && \
    apt-get install -y --assume-yes apt-utils &&\
    apt-get install -y software-properties-common &&\
    apt-get install -y gnupg2 


# RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" > /etc/apt/sources.list.d/cran.list
# RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
# RUN echo "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" > /etc/apt/sources.list

#ENV R_BASE_VERSION 4.0.2

RUN apt update && apt-get update && \
    apt-get install -y --no-install-recommends libopenblas0-pthread \
	littler \
    r-cran-littler \
    r-base \
    r-base-dev \
    r-recommended \
	# r-base=${R_BASE_VERSION}* \
	# r-base-dev=${R_BASE_VERSION}* \
	# r-recommended=${R_BASE_VERSION}* \
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
    libpng-dev \
    libcairo2-dev

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
#RUN apt-get clean && apt-get -y update && apt-get install -y locales && locale-gen en_US.UTF-8
#ENV LANG='en_US.UTF-8' LANGUAGE='en_US.UTF-8' LC_ALL='en_US.UTF-8'

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('DT', 'devtools', 'BiocManager', 'shiny', 'rmarkdown',  'shinydashboard',
        'shinyWidgets', 'shinyBS', 'shinycssloaders', 'shinyjs', 'rjson'), 
        repos='http://cran.rstudio.com/')" && \
    chown shiny:shiny /var/lib/shiny-server

#RUN R -e "install.packages(c('DT', 'devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyWidgets', 'shinyBS', 'shinycssloaders', 'shinyjs', 'rjson'), repos='http://cran.rstudio.com/')"
#RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.9.923-amd64.deb \
#    && gdebi -n shiny-server-1.5.9.923-amd64.deb \
#    && rm -f shiny-server-1.5.9.923-amd64.deb \
#    && chown shiny:shiny /var/lib/shiny-server

# MesKit part:
RUN R -e "BiocManager::install(c('BSgenome', 'GenomeInfoDb', 'org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19'))" 
RUN R -e "devtools::install_github('jokergoo/ComplexHeatmap')"
RUN R -e "devtools::install_github('Niinleslie/MesKit', ref = 'master')"

# Copy the app to the image
COPY inst/shiny /srv/shiny-server/

# Copy further configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf

EXPOSE 3838

# CMD ["/bin/bash"]
CMD ["/usr/bin/shiny-server.sh"]