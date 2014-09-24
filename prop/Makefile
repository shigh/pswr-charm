
TEX_FILE = prop.tex
PDF_FILE = prop.pdf

all : $(PDF_FILE)

$(PDF_FILE) : convg_plot.png $(TEX_FILE)
	pdflatex $(TEX_FILE)

force:
	pdflatex $(TEX_FILE)
	pdflatex $(TEX_FILE)

convg_plot.png:
	python convg_plot.py

clean:
	rm -vf *.png
	rm -vf prop.pdf
	rm -vf *.log *.blg *.bbl *.aux
