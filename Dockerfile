FROM rocker/shiny-verse:4.3

# Install system dependencies (if any)
RUN apt-get update && apt-get install -y libglpk40>=5

# Set the working directory to /package
WORKDIR /srv/shiny-server

# Install your package
RUN R -e 'install.packages(c("ggplot2", "bsicons", "DT", "shinyjs", "shinyFiles", "shinyDisconnect"))'
RUN R -e 'install.packages(c("Seurat"))'
RUN R -e 'remotes::install_github("daattali/shinycssloaders")'
RUN R -e 'install.packages(c("R.utils"))'

RUN rm -rf /srv/shiny-server/sample-apps && \
	rm /srv/shiny-server/index.html

RUN mkdir -p /srv/shiny-server/SOM && \
	echo "SeuratObjectMaker::run_SOM(launch.browser = FALSE, running_locally = FALSE)" > /srv/shiny-server/SOM/app.R && \
	chmod +777 -R /srv/shiny-server/SOM # TODO: change this to smth better

# Copy the current directory contents into the container at /package
# note: this will always be run, thus any lines below will not be cached
# when building the Docker image
COPY . /srv/shiny-server
RUN cp ./shiny-server.conf -t /etc/shiny-server/

RUN R -e "devtools::install('/srv/shiny-server')"
# CMD R -e 'SeuratObjectMaker::SeuratObjectMaker(port = 3838, launch.browser = F)'
