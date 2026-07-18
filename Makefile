UV ?= uv
VERSION := $(shell sed -n 's/^version = "\([^"]*\)"$$/\1/p' pyproject.toml)
SDIST := dist/gpc_census-$(VERSION).tar.gz
SPEC := gpc-census.spec
GEN_SPEC := build/gpc-census.spec
RPM_TOPDIR := $(CURDIR)/build/rpm
REPORT_TEX := report/main.tex
REPORT_PDF := $(REPORT_TEX:.tex=.pdf)

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
	rpmbuild -ba $(GEN_SPEC) --define "_topdir $(RPM_TOPDIR)" --define "_sourcedir $(CURDIR)/dist"
	@echo "==> RPMs under $(RPM_TOPDIR)/RPMS, SRPM under $(RPM_TOPDIR)/SRPMS"

report: $(REPORT_PDF)

$(REPORT_PDF): $(REPORT_TEX)
	latexmk -pdf -interaction=nonstopmode -cd $(REPORT_TEX)

clean:
	rm -rf dist build
	-latexmk -C -cd $(REPORT_TEX) 2>/dev/null
