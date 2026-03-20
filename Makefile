.PHONY: install run run-docking clean help

CONDA_ENV  := t2lead
# Without --no-capture-output, conda run holds logs until exit (looks hung on long runs).
PYTHON     := conda run --no-capture-output -n $(CONDA_ENV) python
DISEASE    ?= breast cancer

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

install:  ## Create conda env and install all dependencies
	conda create -n $(CONDA_ENV) python=3.11 -y
	conda install -n $(CONDA_ENV) -c conda-forge rdkit -y
	$(PYTHON) -m pip install -r requirements.txt
	$(PYTHON) -m pip install -e .
	@echo "\n--- Optional: GPU support ---"
	@echo "  pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128"
	@echo "--- Optional: docking ---"
	@echo "  pip install vina meeko gemmi"
	@echo "--- Optional: MD simulation ---"
	@echo "  conda install -c conda-forge openmm pdbfixer mdtraj openmmforcefields openff-toolkit -y"

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
