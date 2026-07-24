UV ?= uv
VERSION := $(shell sed -n 's/^version = "\([^"]*\)"$$/\1/p' pyproject.toml)
SDIST := dist/gpc_census-$(VERSION).tar.gz
SPEC := gpc-census.spec
GEN_SPEC := build/gpc-census.spec
RPM_TOPDIR := $(CURDIR)/build/rpm
REPORT_MD := results/report/main.md
REPORT_PDF := results/report/main.pdf
REPORT_BIB := results/report/references.bib
REPORT_CSL := results/report/american-physics-society.csl

# pandoc/extra bundles pandoc, a matched pandoc-crossref, and a TeX
# engine, so one pinned image covers the whole markdown-to-PDF build.
# CI runs the same recipe inside the image with PANDOC_RUN=pandoc.
CONTAINER ?= podman
PANDOC_IMAGE := docker.io/pandoc/extra:3.6.4@sha256:6a53f5ac29999b2084691b133546f57a80464a4a3991c15cd1a373133b97e7a7
PANDOC_RUN ?= $(CONTAINER) run --rm -v $(CURDIR):/data:Z -w /data $(PANDOC_IMAGE)
PANDOC_FLAGS := -F pandoc-crossref --citeproc --number-sections

.PHONY: sync test lint build sdist wheel srpm rpm report upgrade clean

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

# pandoc-crossref runs before citeproc so it consumes the [-@sec:x]
# style citations that are cross-references, not bibliography keys;
# anything either filter leaves unresolved surfaces as a citeproc
# "not found" warning or a LaTeX "Reference undefined" warning, and
# pandoc-crossref renders missing targets as literal ?id? marks in the
# PDF, so the guard checks all three and fails the build on any.
$(REPORT_PDF): $(REPORT_MD) $(REPORT_BIB) $(REPORT_CSL)
	mkdir -p build
	$(PANDOC_RUN) $(REPORT_MD) $(PANDOC_FLAGS) -o $(REPORT_PDF) 2> build/report.log \
	  || { cat build/report.log; exit 1; }
	@if grep -Ei 'not found|undefined' build/report.log; then \
	  echo "report: unresolved references or citations, see build/report.log"; exit 1; fi
	@if command -v pdftotext >/dev/null 2>&1 \
	  && pdftotext $(REPORT_PDF) - 2>/dev/null | grep -q '¿'; then \
	  echo "report: unresolved crossref marks in PDF output"; exit 1; fi
	@echo "==> $(REPORT_PDF)"

clean:
	rm -rf dist build data-output data-output.zip
	rm -f $(REPORT_PDF)
