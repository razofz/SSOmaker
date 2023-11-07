# Use an official R base image
FROM r-base:latest

# Install system dependencies (if any)
RUN apt update && apt install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev \
    libfreetype-dev libpng-dev libtiff5-dev libjpeg-dev

# # Set the working directory to /package
# WORKDIR /package

# Copy the current directory contents into the container at /package
ADD . /package

#, repos='http://cran.rstudio.com/')"
# Install your package
RUN R -e "install.packages('devtools')"
# RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install('/package')"
