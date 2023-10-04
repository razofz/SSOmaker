EXAMPLE_DATA_DIR = pbmc_small

clean:
	rm -rf ${EXAMPLE_DATA_DIR}

lint:
	R --quiet -e "styler::style_dir(path=\".\", dry = \"on\")"

format_code:
	R --quiet -e "styler::style_dir(path=\".\")"

download_example_data:
	Rscript download_pbmc_small.R
	gzip --recursive ${EXAMPLE_DATA_DIR}
