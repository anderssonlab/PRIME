.PHONY: extdata

extdata:
	Rscript --vanilla data-raw/make-ctss-extdata.R
