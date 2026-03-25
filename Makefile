.PHONY: help install install-reinvent4 download-reinvent4-prior run run-docking clean clean-all clean-logs

# Default /bin/sh often lacks conda on PATH (non-login, no profile). Use bash -l so
# system profile (~/.profile, /etc/profile) runs; prepend conda bin as fallback.
SHELL := /bin/bash
.SHELLFLAGS := -lc
CONDA_ROOT ?= $(HOME)/miniconda3
export PATH := $(CONDA_ROOT)/bin:$(PATH)

CONDA_ENV  := t2lead
# Without --no-capture-output, conda run holds logs until exit (looks hung on long runs).
PYTHON     := conda run --no-capture-output -n $(CONDA_ENV) python
DISEASE    ?= breast cancer
REINVENT4_DIR ?= /root/REINVENT4
REINVENT4_PRIOR_URL ?= https://zenodo.org/api/records/15641297/files/reinvent.prior/content

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

install:  ## Create conda env and install core dependencies
	conda create -n $(CONDA_ENV) python=3.11 -y
	conda install -n $(CONDA_ENV) -c conda-forge rdkit -y
	$(PYTHON) -m pip install -r requirements.txt
	@echo "\n--- Optional: GPU support ---"
	@echo "  pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128"
	@echo "--- Optional: docking ---"
	@echo "  pip install -e \".[docking]\"   # or: vina meeko gemmi"
	@echo "--- Optional: MD simulation ---"
	@echo "  conda install -c conda-forge openmm pdbfixer mdtraj openmmforcefields openff-toolkit -y"
	@echo "--- Optional: REINVENT4 ---"
	@echo "  make install-reinvent4"

install-reinvent4:  ## Install REINVENT4 into current conda env
	git clone https://github.com/MolecularAI/REINVENT4.git $(REINVENT4_DIR) || true
	$(PYTHON) -m pip install -e $(REINVENT4_DIR)
	@echo "REINVENT4 CLI should be at: /root/miniconda3/envs/$(CONDA_ENV)/bin/reinvent"
	@echo "Next: make download-reinvent4-prior"

download-reinvent4-prior:  ## Download REINVENT4 prior model
	mkdir -p $(REINVENT4_DIR)/priors
	curl -L "$(REINVENT4_PRIOR_URL)" -o $(REINVENT4_DIR)/priors/reinvent.prior
	@echo "Downloaded: $(REINVENT4_DIR)/priors/reinvent.prior"
	@echo "Optional hash-metadata fix:"
	@echo "  $(PYTHON) scripts/fix_reinvent_prior_metadata.py $(REINVENT4_DIR)/priors/reinvent.prior"

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
