pdflatex --shell-escape $1
bibtex $1
pdflatex --shell-escape $1
pdflatex --shell-escape $1

rm $1.aux
rm $1.bbl
rm $1.blg
rm $1.log
rm $1.out

mv $1.pdf ../$1.pdf
