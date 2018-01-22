TEXT = mvt-fields
SOM = supp
REFS = spatial-extremes

all: $(TEXT).pdf response.tex

response.tex: response.md
	pandoc response.md -o response.tex

%.pdf: %.tex $(REFS).bib response.tex
	latexmk -pdf $<

continuous:
	while true; do make --silent; sleep 1; done

dropbox:
	cp mvt-fields.pdf ~/Dropbox/Public/
