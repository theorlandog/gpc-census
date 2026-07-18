UV ?= uv
VERSION := $(shell sed -n 's/^version = "\([^"]*\)"$$/\1/p' pyproject.toml)
SDIST := dist/gpc_census-$(VERSION).tar.gz
SPEC := gpc-census.spec
GEN_SPEC := build/gpc-census.spec
RPM_TOPDIR := $(CURDIR)/build/rpm
REPORT_TEX := results/report/main.tex
REPORT_PDF := $(REPORT_TEX:.tex=.pdf)
REPORT_MD := $(REPORT_TEX:.tex=.md)

# pandoc-crossref binaries are version-locked to the pandoc they were
# built against; v0.3.19 pairs with pandoc 3.6.4 (Fedora 43). Bump the
# two together.
PANDOC_CROSSREF_VERSION := v0.3.19
PANDOC_CROSSREF_SHA256 := d6bdac44dbe9209e0bca0b35ea377cf0a0fd7433cc9a31b9a603512c24317d60
PANDOC_CROSSREF := build/tools/pandoc-crossref
MD_CHECK := build/report-md-check.html

.PHONY: sync test lint build sdist wheel srpm rpm report report-md upgrade clean

sync:
	$(UV) sync

test:
	$(UV) run pytest

lint:
	$(UV) run ruff check

# House rule: only accept dependency releases at least 14 days old.
upgrade:
	$(UV) lock --upgrade --exclude-newer "$$(date -u -d '14 days ago' +%Y-%m-%dT%H:%M:%SZ)"

build:
	$(UV) build

sdist:
	$(UV) build --sdist

wheel:
	$(UV) build --wheel

# The spec's Version line is rewritten to match pyproject.toml so CI-stamped
# versions (e.g. 0.1.0+main.<sha>) flow into the RPM.
$(GEN_SPEC): $(SPEC) pyproject.toml
	mkdir -p build
	sed 's/^Version:.*/Version:        $(VERSION)/' $(SPEC) > $(GEN_SPEC)

srpm: sdist $(GEN_SPEC)
	rpmbuild -bs $(GEN_SPEC) --define "_topdir $(RPM_TOPDIR)" --define "_sourcedir $(CURDIR)/dist"

rpm: sdist $(GEN_SPEC)
	rpmbuild -ba $(GEN_SPEC) --define "_topdir $(RPM_TOPDIR)" --define "_sourcedir $(CURDIR)/dist" \
	  || { command -v dnf >/dev/null 2>&1 \
	       && dnf -y builddep $(RPM_TOPDIR)/SRPMS/*.buildreqs.nosrc.rpm \
	       && rpmbuild -ba $(GEN_SPEC) --define "_topdir $(RPM_TOPDIR)" --define "_sourcedir $(CURDIR)/dist"; }
	@echo "==> RPMs under $(RPM_TOPDIR)/RPMS, SRPM under $(RPM_TOPDIR)/SRPMS"

report: $(REPORT_PDF)

$(REPORT_PDF): $(REPORT_TEX)
	latexmk -pdf -interaction=nonstopmode -cd $(REPORT_TEX)

report-md: $(REPORT_MD)

$(PANDOC_CROSSREF):
	mkdir -p build/tools
	curl -sSfL -o build/tools/pandoc-crossref.tar.xz \
	  https://github.com/lierdakil/pandoc-crossref/releases/download/$(PANDOC_CROSSREF_VERSION)/pandoc-crossref-Linux-X64.tar.xz
	echo "$(PANDOC_CROSSREF_SHA256)  build/tools/pandoc-crossref.tar.xz" | sha256sum -c
	tar -xJf build/tools/pandoc-crossref.tar.xz -C build/tools
	touch $(PANDOC_CROSSREF)

# The render check proves every crossref reference resolves; unresolved
# ones survive as literal [-@sec:x] citations or crossref's ?? marker.
$(REPORT_MD): $(REPORT_TEX) scripts/tex2md.py $(PANDOC_CROSSREF)
	python3 scripts/tex2md.py $(REPORT_TEX) $(REPORT_MD)
	pandoc $(REPORT_MD) -F $(PANDOC_CROSSREF) --number-sections -s --mathjax -o $(MD_CHECK)
	@if grep -qE 'reference-type|@(sec|eq|tbl):|⁇' $(MD_CHECK); then \
	  echo "report-md: unresolved references, inspect $(MD_CHECK)"; exit 1; fi
	@echo "==> $(REPORT_MD) (render check: $(MD_CHECK))"

clean:
	rm -rf dist build data-output data-output.zip
	rm -f $(REPORT_MD)
	-latexmk -C -cd $(REPORT_TEX) 2>/dev/null
