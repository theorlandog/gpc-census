UV ?= uv
VERSION := $(shell sed -n 's/^version = "\([^"]*\)"$$/\1/p' pyproject.toml)
SDIST := dist/gpc_census-$(VERSION).tar.gz
SPEC := gpc-census.spec
GEN_SPEC := build/gpc-census.spec
RPM_TOPDIR := $(CURDIR)/build/rpm
REPORT_TEX := results/report/main.tex
REPORT_PDF := $(REPORT_TEX:.tex=.pdf)
REPORT_MD := $(REPORT_TEX:.tex=.md)

# pandoc/extra bundles pandoc with a matched pandoc-crossref, so one
# pinned image replaces pairing binary versions by hand.
CONTAINER ?= podman
PANDOC_IMAGE := docker.io/pandoc/extra:3.6.4@sha256:6a53f5ac29999b2084691b133546f57a80464a4a3991c15cd1a373133b97e7a7
PANDOC_RUN := $(CONTAINER) run --rm -v $(CURDIR):/data:Z -w /data $(PANDOC_IMAGE)
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

# The render check proves every crossref reference resolves; unresolved
# ones survive as literal [-@sec:x] citations or crossref's ?? marker.
$(REPORT_MD): $(REPORT_TEX) scripts/tex2md.py
	mkdir -p build
	PANDOC="$(PANDOC_RUN)" python3 scripts/tex2md.py $(REPORT_TEX) $(REPORT_MD)
	$(PANDOC_RUN) $(REPORT_MD) -F pandoc-crossref --number-sections -s --mathjax -o $(MD_CHECK)
	@if grep -qE 'reference-type|@(sec|eq|tbl):|⁇' $(MD_CHECK); then \
	  echo "report-md: unresolved references, inspect $(MD_CHECK)"; exit 1; fi
	@echo "==> $(REPORT_MD) (render check: $(MD_CHECK))"

clean:
	rm -rf dist build data-output data-output.zip
	rm -f $(REPORT_MD)
	-latexmk -C -cd $(REPORT_TEX) 2>/dev/null
