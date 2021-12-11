DOC_BASENAME :=  "HERVfinder"
DOC :=  "$(DOC_BASENAME).tex"
LATEX_OPTS := -shell-escape -interaction=nonstopmode

all: pdf

dvi:
	rm -rf _minted*
	latex $(LATEX_OPTS) $(DOC)
	latex $(LATEX_OPTS) $(DOC)

ps:dvi
	dvips "$(DOC_BASENAME).dvi"

pdf:
	rm -rf _minted*
	pdflatex $(LATEX_OPTS) $(DOC)
	pdflatex $(LATEX_OPTS) $(DOC)

answer_len_distrib.png:
	Rscript get_len_distrib.R

clean:
	rm -rf *.aux *.fdb_latexmk *.fls *.log *.out *.pdf *.synctex* *.toc
