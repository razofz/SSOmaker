FROM rocker/shiny-verse:4.3

# Install system dependencies (if any)
RUN apt-get update && apt-get install -y libglpk40>=5

# Set the working directory to /package
WORKDIR /package

# Install your package
RUN R -e 'install.packages(c("ggplot2", "bsicons", "DT", "shinyjs", "shinyFiles", "shinyDisconnect"))'
RUN R -e 'install.packages(c("Seurat"))'
RUN R -e 'remotes::install_github("daattali/shinycssloaders")'

# Copy the current directory contents into the container at /package
COPY . /package

RUN R -e "devtools::install('/package')"
# RUN R -e 'SeuratObjectMaker(port = 3838, launch.browser = F)'
