.PHONY: SISMatlab_v11.pdf all clean

all: SISMatlab_v11.pdf

%.tex: %.raw
	./raw2tex $< > $@

%.tex: %.dat
	./dat2tex $< > $@


SISMatlab_v11.pdf: SISMatlab_v11.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make SISMatlab_v11.tex

clean:
	latexmk -CA

