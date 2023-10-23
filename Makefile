EXAMPLE_DATA_DIR = pbmc_small
PROJECT_NAME = SOM
PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

MAMBA := $(shell which mamba)
CONDA := $(shell which conda)

ifdef MAMBA
	CONDA_COMMAND := mamba
else ifdef CONDA
	CONDA_COMMAND := conda
else
	$(error Neither conda nor mamba is installed. Please install at least conda.)
endif

## Set up python interpreter environment
environment:
ifdef CONDA_COMMAND
	$(CONDA_COMMAND) create --name $(PROJECT_NAME) --file no-build-env.yaml -y
	@echo ">>> New conda env created. Activate with:\nconda activate $(PROJECT_NAME)"
else
	@echo ">>> Neither mamba nor conda found, please install at least conda."
endif

.PHONY: delete_environment
CONDA_ENV_LOCATION := $(shell conda info --envs | grep "^$(PROJECT_NAME)" | tr -s ' ' ' ' | cut -d ' ' -f 2)
## Delete the virtual (conda) environment for this project
delete_environment: clean
	@echo ">>> Removing environment $(PROJECT_NAME)"
	@echo ">>> in location $(CONDA_ENV_LOCATION)"
	$(CONDA_COMMAND) env remove --name $(PROJECT_NAME)

clean:
	rm -rf ${EXAMPLE_DATA_DIR}

lint:
	R --quiet -e "styler::style_dir(path=\".\", dry = \"on\")"

format_code:
	R --quiet -e "styler::style_dir(path=\".\")"

download_example_data:
	Rscript pbmc_small_download.R
	gzip --recursive ${EXAMPLE_DATA_DIR}

install_extra_r_packages:
	R --quiet -e 'remotes::install_github("daattali/shinycssloaders")'
	R --quiet -e 'suppressWarnings(if (!require("bsicons")) install.packages("bsicons", repos="https://ftp.acc.umu.se/mirror/CRAN/"))'
