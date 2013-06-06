BIB=Antarctica2013.bib

PDF=Antarctica2013-all-control.pdf 

TEX=$(patsubst  %.pdf, %.tex,  $(PDF) )
TOC=$(patsubst  %.pdf, %.toc,  $(PDF) )
BBL=$(patsubst  %.pdf, %.bbl,  $(PDF) )
LOG=$(patsubst  %.pdf, %.log,  $(PDF) )
AUX=$(patsubst  %.pdf, %.aux,  $(PDF) )
BLG=$(patsubst  %.pdf, %.blg,  $(PDF) )

all: $(PDF)

%.aux: %.tex
	pdflatex $<

%.bbl: %.aux $(BIB)
	bibtex $(patsubst  %.tex, %,  $< ) 

%.pdf: %.tex %.aux %.bbl
	pdflatex $< && pdflatex $< 

clean: 
	$(RM) $(AUX) $(PDF) $(LOG) $(BBL) $(TOC) $(BLG)
