.PHONY: SISMatlab_v10.pdf all clean

all: SISMatlab_v10.pdf

%.tex: %.raw
	./raw2tex $< > $@

%.tex: %.dat
	./dat2tex $< > $@


SISMatlab_v10.pdf: SISMatlab_v10.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make SISMatlab_v10.tex

clean:
	latexmk -CA

