PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGDATE := $(shell sed -n "s/Date: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

GITDATE=$(shell (git log -1 --date=short --pretty=format:"%ad"))
GITVERS=$(shell (date -d `git log -1 --date=short --pretty=format:"%ad"` +%Y%m%d))

R_FILES := $(wildcard $(PKGSRC)/R/*.R)
HELP_FILES := $(wildcard $(PKGSRC)/man/*.Rd)

all: build

.PHONY: all release roxygen

README.md: DESCRIPTION
	sed -i 's/Version: *\([^ ]*\)/Version: $(PKGVERS)/' README.md
	sed -i 's/Date: *\([^ ]*\)/Date: $(PKGDATE)/' README.md

NEWS: NEWS.md
	sed 's/^# / /' NEWS.md > NEWS
	sed -i 's/^##//' NEWS

docs: $(HELP_FILES) README.md
	R --vanilla --silent -e "options(repos='http://cran.r-project.org'); pkgdown::build_site(preview=FALSE)"

roxygen: $(HELP_FILES)

$(HELP_FILES): $(R_FILES)
	R --vanilla --silent -e "library(devtools);" \
		-e "document(roclets='rd')"

update:
	sed -i 's/Date: *\([^ ]*\)/Date: $(GITDATE)/' DESCRIPTION

release: build docs

build:
	cd ..;\
	R --vanilla CMD build $(PKGSRC) --compact-vignettes="gs+qpdf" --resave-data

buildNV: README.md NEWS
	cd ..;\
	R --vanilla CMD build $(PKGSRC) --no-build-vignettes

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R --vanilla CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

clean:
	cd ..;\
	rm -rf $(PKGNAME).Rcheck $(PKGNAME)_$(PKGVERS).tar.gz
