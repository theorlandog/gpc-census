UV ?= uv
VERSION := $(shell sed -n 's/^version = "\([^"]*\)"$$/\1/p' pyproject.toml)
SDIST := dist/gpc_census-$(VERSION).tar.gz
SPEC := gpc-census.spec
GEN_SPEC := build/gpc-census.spec
RPM_TOPDIR := $(CURDIR)/build/rpm
REPORT_MD := results/report/main.md
REPORT_PDF := results/report/main.pdf
REPORT_TEX := results/report/main.tex
REPORT_BIB := results/report/references.bib
REPORT_CSL := results/report/physics-numeric.csl
ANON_DIR := results/report/anonymized
ANON_MD := $(ANON_DIR)/main.md
ANON_PDF := $(ANON_DIR)/main.pdf
ANON_TEX := $(ANON_DIR)/main.tex
STATE_VERIFIER := scripts/verify_states_standalone.py
NO_DESIGN_VERIFIER := scripts/verify_interference_certificates_standalone.py
VERTEX_VERIFIER := scripts/verify_vertex_exhaustion_standalone.py
SUPPLEMENT_BUILDER := scripts/build_paper1_supplement.py
SUPPLEMENT_ZIP := build/gpc-census-paper1-supplement.zip
SUPPLEMENT_PYTHON := 3.12.13
SUPPLEMENTARY_MD := results/report/supplementary_material.md
SUPPLEMENTARY_PDF := results/report/supplementary_material.pdf

# pandoc/extra bundles pandoc, citeproc support, and a TeX engine, so one
# pinned image covers the whole markdown-to-PDF build.
# CI runs the same recipe inside the image with PANDOC_RUN=pandoc.
CONTAINER ?= podman
PANDOC_IMAGE := docker.io/pandoc/extra:3.6.4@sha256:6a53f5ac29999b2084691b133546f57a80464a4a3991c15cd1a373133b97e7a7
PANDOC_RUN ?= $(CONTAINER) run --rm -v $(CURDIR):/data:Z -w /data $(PANDOC_IMAGE)
PANDOC_FLAGS := --citeproc --number-sections

.PHONY: sync test checksums checksums-check verify-paper verify-data emit note paper1-supplement paper1-supplement-locked supplementary-material lint build sdist wheel srpm rpm report report-tex anonymize report-anon report-anon-tex upgrade clean

sync:
	$(UV) sync

test:
	$(UV) run pytest

# Regenerate the data manifest. This is the resolution whenever SHA256SUMS is
# reported stale: never edit the file, rebuild it.
checksums:
	$(UV) run python scripts/update_data_checksums.py

# Fail if the manifest is stale, without rewriting it. Cheap enough to run
# before pushing, and the fast way to catch a merge that resolved the manifest
# textually while the files it describes merged independently.
checksums-check:
	$(UV) run python scripts/update_data_checksums.py --check

verify-paper:
	$(UV) run python $(STATE_VERIFIER) results/data/states.jsonl
	python3 $(NO_DESIGN_VERIFIER) results/data/interference_certificates.json
	python3 $(VERTEX_VERIFIER)
	python3 scripts/check_manuscript_counts.py

# The manuscript checker is scoped to Paper 1, so the gauge, cascade, orbit,
# holonomy, natural-orbital, v103 and fiber-symmetry artifacts it released
# keep their own gate here.
verify-data:
	$(UV) run python scripts/check_data_consistency.py
	$(UV) run python scripts/audit_data_completeness.py
	$(UV) run python scripts/build_census_master.py --check
	$(UV) run python scripts/emit_doc_tables.py --check

# Regenerate everything the sync invariant owns: the census summary and every
# emitted doc table. Run this, never hand-edit a block between sync markers.
emit:
	$(UV) run python scripts/build_census_master.py
	$(UV) run python scripts/emit_doc_tables.py --write
	$(UV) run python scripts/update_data_checksums.py

paper1-supplement:
	$(UV) run --python $(SUPPLEMENT_PYTHON) python $(SUPPLEMENT_BUILDER) --output $(SUPPLEMENT_ZIP)

paper1-supplement-locked:
	$(UV) run --python $(SUPPLEMENT_PYTHON) python $(SUPPLEMENT_BUILDER) \
	  --output $(SUPPLEMENT_ZIP) --test-locked-environment

supplementary-material: $(SUPPLEMENTARY_PDF)

$(SUPPLEMENTARY_PDF): $(SUPPLEMENTARY_MD)
	mkdir -p build
	$(PANDOC_RUN) $(SUPPLEMENTARY_MD) --number-sections -o $(SUPPLEMENTARY_PDF) \
	  2> build/supplementary-material.log \
	  || { cat build/supplementary-material.log; exit 1; }

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

# Citeproc resolves the bibliography. Internal theorem links are ordinary
# Markdown anchors, while the two tables are referred to by their stable order.
# Any unresolved citation surfaces as a citeproc warning and fails the build.
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

# The level-5 companion note: same pinned image, same citation discipline,
# same failure-on-unresolved-citation rule as the main report.
NOTE_MD := results/report/level5/main.md
NOTE_PDF := results/report/level5/main.pdf
NOTE_BIB := results/report/level5/references.bib

note: $(NOTE_PDF)

$(NOTE_PDF): $(NOTE_MD) $(NOTE_BIB) $(REPORT_CSL)
	mkdir -p build
	$(PANDOC_RUN) $(NOTE_MD) $(PANDOC_FLAGS) -o $(NOTE_PDF) 2> build/note.log \
	  || { cat build/note.log; exit 1; }
	@if grep -Ei 'not found|undefined' build/note.log; then \
	  echo "note: unresolved references or citations, see build/note.log"; exit 1; fi
	@echo "==> $(NOTE_PDF)"

report-tex: $(REPORT_TEX)

# Standalone LaTeX for journal upload, from the same pinned image and
# citation chain as the PDF so both outputs stay in lockstep. citeproc
# bakes the formatted bibliography into the .tex, so the file is
# self-contained (no .bbl/.bib needed at submission). For journals that
# require natbib \citep markup and a .bib instead, swap --citeproc for
# --natbib in PANDOC_FLAGS for this target and upload $(REPORT_BIB)
# alongside. Unlike PDF output, .tex is not standalone by default,
# hence -s. The guard mirrors the PDF one and rejects unresolved markers.
$(REPORT_TEX): $(REPORT_MD) $(REPORT_BIB) $(REPORT_CSL)
	mkdir -p build
	$(PANDOC_RUN) $(REPORT_MD) $(PANDOC_FLAGS) -s -o $(REPORT_TEX) 2> build/report-tex.log \
	  || { cat build/report-tex.log; exit 1; }
	@if grep -Ei 'not found|undefined' build/report-tex.log; then \
	  echo "report-tex: unresolved references or citations, see build/report-tex.log"; exit 1; fi
	@if grep -q '?[A-Za-z:-]*?}' $(REPORT_TEX); then \
	  echo "report-tex: unresolved crossref marks in TeX output"; exit 1; fi
	@echo "==> $(REPORT_TEX)"

# Double-anonymous review copy for journal submission. The folder is
# generated from the master manuscript, never edited by hand; the
# generator strips the author list, affiliation footnote, preprint DOI,
# and personal acknowledgments, and fails if any identifying string
# survives. tests/test_anonymized_report.py guards against drift.
anonymize:
	python3 scripts/make_anonymized_report.py

report-anon: anonymize
	mkdir -p build
	$(PANDOC_RUN) $(ANON_MD) $(PANDOC_FLAGS) -o $(ANON_PDF) 2> build/report-anon.log \
	  || { cat build/report-anon.log; exit 1; }
	@if grep -Ei 'not found|undefined' build/report-anon.log; then \
	  echo "report-anon: unresolved references or citations, see build/report-anon.log"; exit 1; fi
	@if command -v pdftotext >/dev/null 2>&1 \
	  && pdftotext $(ANON_PDF) - 2>/dev/null | grep -q '¿'; then \
	  echo "report-anon: unresolved crossref marks in PDF output"; exit 1; fi
	@echo "==> $(ANON_PDF)"

# Standalone LaTeX of the anonymized copy, for journals that ask for
# source at submission. Same citation chain and guards as report-tex;
# citeproc bakes the formatted bibliography into the .tex, and the
# folder's references.bib ships alongside for editors who want the raw
# database.
report-anon-tex: anonymize
	mkdir -p build
	$(PANDOC_RUN) $(ANON_MD) $(PANDOC_FLAGS) -s -o $(ANON_TEX) 2> build/report-anon-tex.log \
	  || { cat build/report-anon-tex.log; exit 1; }
	@if grep -Ei 'not found|undefined' build/report-anon-tex.log; then \
	  echo "report-anon-tex: unresolved references or citations, see build/report-anon-tex.log"; exit 1; fi
	@if grep -q '?[A-Za-z:-]*?}' $(ANON_TEX); then \
	  echo "report-anon-tex: unresolved crossref marks in TeX output"; exit 1; fi
	@echo "==> $(ANON_TEX)"

clean:
	rm -rf dist build data-output data-output.zip
	rm -f $(REPORT_PDF) $(REPORT_TEX) $(ANON_PDF) $(ANON_TEX) $(SUPPLEMENTARY_PDF)
