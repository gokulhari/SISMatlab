.PHONY: feedback_svdfr_v24.pdf all clean

all: SISMatlab_v7.pdf

%.tex: %.raw
	./raw2tex $< > $@

%.tex: %.dat
	./dat2tex $< > $@


SISMatlab_v7.pdf: SISMatlab_v7.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make SISMatlab_v7.tex

clean:
	latexmk -CA

