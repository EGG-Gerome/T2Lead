.PHONY: help install install-conda-stack install-pip-stack install-torch \
	install-crem-db install-reinvent4-full install-env \
	install-reinvent4 download-reinvent4-prior run run-docking clean clean-all clean-logs

# Default /bin/sh often lacks conda on PATH (non-login, no profile). Use bash -l so
# system profile (~/.profile, /etc/profile) runs; prepend conda bin as fallback.
SHELL := /bin/bash
.SHELLFLAGS := -lc
CONDA_ROOT ?= $(HOME)/miniconda3
export PATH := $(CONDA_ROOT)/bin:$(PATH)

CONDA_ENV  := t2lead
# Without --no-capture-output, conda run holds logs until exit (looks hung on long runs).
PYTHON     := conda run --no-capture-output -n $(CONDA_ENV) python
# Conda env is isolated; suppress pip's generic "running as root" warning in containers.
export PIP_ROOT_USER_ACTION := ignore
DISEASE    ?= breast cancer
# PyTorch CUDA wheels: override for your GPU (Blackwell cu128 ⚠ OpenMM MD falls back to CPU — see README, Ampere cu118, CPU: omit index).
TORCH_INDEX_URL ?= https://download.pytorch.org/whl/cu124
REINVENT4_DIR ?= $(CURDIR)/REINVENT4
REINVENT4_PRIOR_URL ?= https://zenodo.org/api/records/15641297/files/reinvent.prior/content
CREM_DB_FILE ?= $(CURDIR)/data/crem_db/chembl33_sa25_f5.db

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-26s\033[0m %s\n", $$1, $$2}'

install: install-conda-stack install-pip-stack install-torch install-crem-db install-reinvent4-full install-env  ## Full install (RDKit, OpenMM stack, docking, CReM, PyTorch, CReM DB, REINVENT4, .env)
	@echo ""
	@echo "Full install finished. Run: make run DISEASE=\"your disease\""

install-conda-stack:  ## Create conda env (if missing) and conda-forge scientific stack
	@if ! conda env list | grep -qE '^$(CONDA_ENV)[[:space:]]'; then \
		echo "Creating conda env $(CONDA_ENV) ..."; \
		conda create -n $(CONDA_ENV) python=3.11 -y; \
	else \
		echo "Conda env $(CONDA_ENV) already exists; skipping create."; \
	fi
	conda install -n $(CONDA_ENV) -c conda-forge \
		rdkit openmm pdbfixer mdtraj openmmforcefields openff-toolkit -y

install-pip-stack:  ## Editable install + docking + CReM (crem, vina, meeko, gemmi)
	$(PYTHON) -m pip install -U pip
	$(PYTHON) -m pip install -e ".[docking,h2l]"

install-torch:  ## PyTorch + torchvision (set TORCH_INDEX_URL for your CUDA / CPU)
	$(PYTHON) -m pip install torch torchvision --index-url $(TORCH_INDEX_URL)

install-crem-db:  ## Download and extract default CReM fragment database (Zenodo)
	bash scripts/download_crem_db.sh

install-reinvent4-full:  ## Clone REINVENT4 (if needed), pip install, download prior, optional metadata fix
	@if [ ! -f "$(REINVENT4_DIR)/pyproject.toml" ]; then \
		if [ -d "$(REINVENT4_DIR)" ]; then \
			rmdir "$(REINVENT4_DIR)" 2>/dev/null || { echo "ERROR: $(REINVENT4_DIR) exists and is not empty (or not a REINVENT4 tree); remove it or set REINVENT4_DIR." >&2; exit 1; }; \
		fi; \
		git clone https://github.com/MolecularAI/REINVENT4.git "$(REINVENT4_DIR)"; \
	fi
	$(PYTHON) -m pip install -e "$(REINVENT4_DIR)"
	mkdir -p "$(REINVENT4_DIR)/priors"
	@if [ ! -s "$(REINVENT4_DIR)/priors/reinvent.prior" ]; then \
		curl -L "$(REINVENT4_PRIOR_URL)" -o "$(REINVENT4_DIR)/priors/reinvent.prior"; \
	else \
		echo "REINVENT4 prior already present; skipping download."; \
	fi
	@$(PYTHON) scripts/fix_reinvent_prior_metadata.py "$(REINVENT4_DIR)/priors/reinvent.prior" || \
		echo "Note: prior metadata fix skipped or failed (often harmless)."

install-env:  ## Create .env from .env.example if needed and write DP_* paths
	@test -f .env || cp .env.example .env
	@CONDA_P=$$(conda run --no-capture-output -n $(CONDA_ENV) bash -c 'echo $$CONDA_PREFIX'); \
	PRIOR=$$(realpath "$(REINVENT4_DIR)/priors/reinvent.prior"); \
	CREM=$$(realpath "$(CREM_DB_FILE)"); \
	$(PYTHON) scripts/bootstrap_env.py --env-file .env --conda-prefix "$$CONDA_P" \
		--prior-path "$$PRIOR" --crem-db "$$CREM"

install-reinvent4: install-reinvent4-full  ## (Partial) REINVENT4 only — alias of install-reinvent4-full

download-reinvent4-prior:  ## (Partial) Download prior into REINVENT4_DIR/priors
	mkdir -p "$(REINVENT4_DIR)/priors"
	curl -L "$(REINVENT4_PRIOR_URL)" -o "$(REINVENT4_DIR)/priors/reinvent.prior"
	@echo "Downloaded: $(REINVENT4_DIR)/priors/reinvent.prior"

run:  ## Run full pipeline (set DISEASE="lung cancer" to override)
	$(PYTHON) scripts/run_pipeline.py --disease "$(DISEASE)" -v

run-docking:  ## Run in docking-only mode (set TARGET=CHEMBLXXXX)
	$(PYTHON) scripts/run_pipeline.py --target $(TARGET) --docking-only -v

clean:  ## Remove cached data (fingerprints, models, scored candidates)
	rm -rf data/fp_cache data/*/model_cache data/*/scored_candidates.csv
	rm -rf data/*/docking_poses data/*/md_trajectories
	@echo "Caches cleared. Raw crawl data (molecules/activities CSVs) preserved."

clean-all:  ## Remove ALL output data (including crawled ChEMBL data)
	rm -rf data/*
	@echo "All output data removed."

clean-logs:  ## Remove log files
	rm -rf data/logs
	@echo "Logs removed."
