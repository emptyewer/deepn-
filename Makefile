CMAKE       ?= cmake
ifeq ($(PRESET),release)
    BUILD_DIR ?= build-release
else ifeq ($(PRESET),debug)
    BUILD_DIR ?= build-debug
else
    BUILD_DIR ?= build-debug
endif
PRESET      ?= default
JOBS        ?= $(shell sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)
BLAST_DB    ?= data/blast_db
MAKEBLASTDB  = ncbi/build/bin/makeblastdb
MAKEMBINDEX  = ncbi/build/bin/makembindex

CODESIGN_IDENTITY ?= Developer ID Application: Proteverse LLC (4D867DTJWY)

VERSION     ?= $(shell git describe --tags --always --dirty 2>/dev/null || echo dev)
ARCHIVE     ?= DEEPN++-$(VERSION)-macOS.tar.gz

.PHONY: help submodules configure compile build build-quick blast-db blast-db-copy codesign package clean rebuild distclean

.DEFAULT_GOAL := help

help:
	@echo "DEEPN++ Build System"
	@echo ""
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@echo "  build        Build all + BLAST DBs + codesign + package ($(JOBS) jobs)"
	@echo "  build-quick  Build all, use prebuilt BLAST databases"
	@echo "  codesign     Code-sign DEEPN++.app with Developer ID"
	@echo "  package      Create $(ARCHIVE) from signed app"
	@echo "  blast-db     Generate BLAST databases from FASTA files in $(BLAST_DB)/"
	@echo "  configure    Run CMake configure only"
	@echo "  clean        Clean build artifacts"
	@echo "  rebuild      Clean + build"
	@echo "  distclean    Remove entire build directory"
	@echo ""
	@echo "Options:"
	@echo "  JOBS=N       Parallel jobs (default: $(JOBS))"
	@echo "  PRESET=x     CMake preset (default: $(PRESET))"
	@echo "                 debug   - Qt 6.11 debug build"
	@echo "                 release - Qt 6.11 static release build"
	@echo "  BLAST_DB=dir FASTA source directory (default: $(BLAST_DB))"

submodules:
	git submodule update --init zlib ncbi blat xlsx

configure: submodules
	$(CMAKE) --preset $(PRESET)

compile: configure
	$(CMAKE) --build $(BUILD_DIR) -j$(JOBS)

blast-db:
	@mkdir -p $(BLAST_DB)
	@found=0; \
	for db in $(BLAST_DB)/*.fasta $(BLAST_DB)/*.fa; do \
		[ -f "$$db" ] || continue; \
		found=1; \
		echo "Generating BLAST database for $$db"; \
		$(MAKEBLASTDB) -in "$$db" -dbtype nucl; \
		$(MAKEMBINDEX) -input "$$db"; \
	done; \
	[ $$found -eq 0 ] && echo "No FASTA files found in $(BLAST_DB)/" || true

blast-db-copy:
	@cp $(BLAST_DB)/*.* $(BUILD_DIR)/junction_dice/JunctionDice++.app/Contents/Data/ 2>/dev/null || true
	@cp -R $(BUILD_DIR)/junction_dice/JunctionDice++.app $(BUILD_DIR)/deepn/DEEPN++.app/Contents/Resources/ 2>/dev/null || true

codesign:
	$(CMAKE) \
		"-DAPP_PATH=$(BUILD_DIR)/deepn/DEEPN++.app" \
		"-DIDENTITY=$(CODESIGN_IDENTITY)" \
		"-DENTITLEMENTS=$(CURDIR)/cmake/entitlements.plist" \
		-P cmake/codesign.cmake
	codesign --verify --deep --strict $(BUILD_DIR)/deepn/DEEPN++.app

package:
	@echo "Packaging DEEPN++.app → $(ARCHIVE)"
	tar -czf $(BUILD_DIR)/$(ARCHIVE) -C $(BUILD_DIR)/deepn DEEPN++.app
	@echo "Created $(BUILD_DIR)/$(ARCHIVE) ($(shell du -h $(BUILD_DIR)/$(ARCHIVE) 2>/dev/null | cut -f1))"

build: compile blast-db blast-db-copy codesign package

build-quick: compile

clean:
	$(CMAKE) --build $(BUILD_DIR) --target clean 2>/dev/null || true

rebuild: clean build

distclean:
	rm -rf $(BUILD_DIR)
