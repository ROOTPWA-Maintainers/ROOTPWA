SHELL:=$(shell which bash)

DOCUMENT_NAME=main

SOURCES=likelihood.tex

STYLE_FILES=packages.sty \
            macros.sty

LATEX_CMD=pdflatex -synctex=1
BIBTEX_CMD=bibtex


default: $(DOCUMENT_NAME).pdf

view:	$(DOCUMENT_NAME).pdf
	acroread $(DOCUMENT_NAME).pdf &

$(DOCUMENT_NAME).pdf: $(DOCUMENT_NAME).tex $(SOURCES) $(STYLE_FILES) #$(DOCUMENT_NAME).bbl
	$(LATEX_CMD) $(DOCUMENT_NAME)
	$(LATEX_CMD) $(DOCUMENT_NAME)

#$(DOCUMENT_NAME).bbl: $(DOCUMENT_NAME).tex $(SOURCES) $(STYLE_FILES) ../spectroscopy.bib
#	$(LATEX_CMD) $(DOCUMENT_NAME)
#	$(BIBTEX_CMD) $(DOCUMENT_NAME)

clean:
	rm -fv $(DOCUMENT_NAME).{pdf,synctex.gz,aux,log,out,soc,toc,tdo,bbl,blg}
	for i in $(SOURCES); do \
		rm -fv $${i%.tex}.aux; \
	done

emacs:
	emacs $(DOCUMENT_NAME).tex  $(SOURCES) &
