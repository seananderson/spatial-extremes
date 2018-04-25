TEXT = mvt-fields
SOM = supp
REFS = spatial-extremes

all: $(TEXT).pdf Appendix-S1.pdf

response.tex: response.md
	pandoc response.md -o response.tex

Appendix-S1.pdf: Appendix-S1.tex
	latexmk -pdf $<

%.pdf: %.tex $(REFS).bib response.tex
	latexmk -pdf $<

continuous:
	while true; do make --silent; sleep 1; done

clean:
	latexmk -C Appendix-S1.tex
	latexmk -C mvt-fields.tex


dropbox:
	cp mvt-fields.pdf ~/Dropbox/Public/
