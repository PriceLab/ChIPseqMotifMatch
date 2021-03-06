all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes ChIPseqMotifMatch)

install:
	(cd ..; R CMD INSTALL ChIPseqMotifMatch)

check:
	(cd ..; R CMD check --ignore-vignettes --no-manual `ls -t ChIPseqMotifMatch_* | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

