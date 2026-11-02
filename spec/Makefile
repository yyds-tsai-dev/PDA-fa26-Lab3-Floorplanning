# Builds Lab3.pdf from lab3.tex with xelatex (run twice for references).
FIGS = figs/toy.tex figs/toy_partial.tex

all: Lab3.pdf

Lab3.pdf: lab3.tex $(FIGS) ../cases/toy/toy.in ../cases/toy/toy.out ../cases/toy/partial.out
	xelatex -jobname=Lab3 -interaction=nonstopmode -halt-on-error lab3.tex
	xelatex -jobname=Lab3 -interaction=nonstopmode -halt-on-error lab3.tex

# Regenerate the toy figures (needs python3); the generated files are kept in the tree.
figs:
	python3 ../tools/plot_tikz.py ../cases/toy/toy.in ../cases/toy/toy.out figs/toy.tex --nets
	python3 ../tools/plot_tikz.py ../cases/toy/toy.in ../cases/toy/partial.out figs/toy_partial.tex --width 0.45

clean:
	rm -f Lab3.aux Lab3.log Lab3.out Lab3.toc

distclean: clean
	rm -f Lab3.pdf

.PHONY: all figs clean distclean
