# T2Lead — reproducible runtime (CPU PyTorch by default; use host CUDA stack for GPU).
FROM continuumio/miniconda3:24.11.1-0

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_ROOT_USER_ACTION=ignore

RUN apt-get update && apt-get install -y --no-install-recommends \
    git curl ca-certificates openjdk-21-jre-headless \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

RUN conda create -n t2lead python=3.11 -y \
    && conda install -n t2lead -c conda-forge \
        rdkit openmm pdbfixer mdtraj openmmforcefields openff-toolkit ambertools -y \
    && conda run -n t2lead pip install --no-cache-dir -U pip \
    && conda run -n t2lead pip install --no-cache-dir -e ".[docking,h2l]" \
    && conda run -n t2lead pip install --no-cache-dir \
        torch torchvision --index-url https://download.pytorch.org/whl/cpu

# Nextflow for nf-core/sarek (optional variant-calling upstream)
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

ENV PATH="/opt/conda/envs/t2lead/bin:${PATH}" \
    CONDA_DEFAULT_ENV=t2lead

CMD ["python", "scripts/run_pipeline.py", "--help"]
