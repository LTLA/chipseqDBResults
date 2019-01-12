targets: all

RCMD=R

# Obtaining the vignette sources.
src:
	 git clone https://github.com/LTLA/chipseqDB src

update: src
	cd src && git pull

# Moving the scripts to avoid directory chaos.
%.Rmd: src/vignettes/%.Rmd
	cp $< .

ref.bib: src/vignettes/ref.bib
	cp src/vignettes/ref.bib .

# Defining the HTML outputs.
%.knit.md : %.Rmd ref.bib
	${RCMD} --no-save --slave -e "rmarkdown::render('$<', clean=FALSE)"

all: intro.knit.md cbp.knit.md h3k9ac.knit.md

# Cleaning commands.
clean:
	rm -rf *.html *_files *.Rmd 
	git checkout *.knit.md
	cd src && git reset --hard HEAD
