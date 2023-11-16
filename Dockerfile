FROM rocker/shiny-verse:4.3

# Install system dependencies (if any)
RUN apt-get update && apt-get install -y libglpk40>=5

WORKDIR /srv/shiny-server

# Install your package
RUN R -e 'install.packages(c("ggplot2", "bsicons", "DT", "shinyjs", "shinyFiles", "shinyDisconnect"))'
RUN R -e 'install.packages(c("Seurat"))'
RUN R -e 'remotes::install_github("daattali/shinycssloaders")'
RUN R -e 'install.packages(c("R.utils"))'

# RUN rm -rf /srv/shiny-server/*
RUN rm -rf /srv/shiny-server/[0-9][0-9]_*
# RUN /bin/bash -c "find /srv/shiny-server -mindepth 1 ! -name 'app_cache' -exec rm -rf {} \;"
# RUN rm -rf /srv/shiny-server/sample-apps && \
# 	rm /srv/shiny-server/index.html

RUN mkdir -p /srv/shiny-server/SOM && \
	echo "SeuratObjectMaker::run_SOM(launch.browser = FALSE, running_locally = FALSE)" > /srv/shiny-server/app.R && \
	echo "SeuratObjectMaker::run_SOM(launch.browser = FALSE, running_locally = FALSE)" > /srv/shiny-server/SOM/app.R && \
	chmod u=rwx,go=u-w /srv/shiny-server/app.R && \ 
	chmod u=rwx,go=u-w -R /srv/shiny-server/SOM

# Copy the current directory contents into the container at /package
# note: this will always be run, thus any lines below will not be cached
# when building the Docker image
COPY . /package
# COPY . /srv/shiny-server
RUN cp /package/shiny-server.conf -t /etc/shiny-server/
# RUN cp ./shiny-server.conf -t /etc/shiny-server/

RUN R -e "devtools::install('/package')"
# RUN R -e "devtools::install('/srv/shiny-server')"
# CMD R -e 'SeuratObjectMaker::SeuratObjectMaker(port = 3838, launch.browser = F)'
