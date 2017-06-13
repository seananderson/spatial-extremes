TEXT = mvt-fields
SOM = supp
REFS = spatial-extremes

all: $(TEXT).pdf

%.pdf: %.tex $(REFS).bib
	latexmk -pdf $<

continuous:
	while true; do make --silent; sleep 1; done

dropbox:
	cp mvt-fields.pdf ~/Dropbox/
