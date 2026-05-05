# UPDD Dockerfile (CPU-only, validation + small-unit-test path)
# Generated: 2026-05-05
# Anchor: commit fd637a8 (post Phase 1+2 unification)
#
# Purpose:
#   Reproducible CPU-only environment for code-review, unit-test, and
#   schema-regression workflows. Does NOT include gpu4pyscf or cupy
#   (GPU stack is too large for general-purpose Docker distribution and
#   requires NVIDIA-Container-Toolkit at host-side; see Dockerfile.gpu
#   for the GPU production image).
#
#   Suitable for:
#     - GitHub Actions CI (no GPU runner needed)
#     - External reviewers who want to inspect schema / R-15-R-16 charge
#       guards / R-17 cofactor preservation / branched_ddg statistics
#     - bioRxiv / Zenodo reproducibility audit
#
#   NOT suitable for:
#     - Production MD / QM-MM SCF runs (use Dockerfile.gpu or native conda env)
#     - Full pytest suite (~ 30+ tests require openmm + mdtraj real compute paths
#       that complete on CPU but are slower than GPU); pytest baseline 284 pass /
#       9 skip / 1 fail is achievable but takes ~10 min on CPU vs ~1 min on GPU.
#
# Build:
#   docker build -t updd:cpu -f Dockerfile .
#
# Run pytest:
#   docker run --rm -v $(pwd):/workspace updd:cpu \
#       /opt/conda/envs/qmmm/bin/python -m pytest tests/

FROM continuumio/miniconda3:24.5.0-0

LABEL maintainer="UPDD project"
LABEL description="UPDD CPU-only validation environment (qmmm conda env, no GPU stack)"
LABEL version="0.7.1"

# System deps for ambertools + matplotlib + general scientific stack
RUN apt-get update && apt-get install -y \
    build-essential \
    libxrender1 \
    libxext6 \
    libsm6 \
    libgl1-mesa-glx \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Copy environment specification first (for Docker layer caching)
COPY environment.yml /tmp/environment.yml

# Build the qmmm conda env without the GPU pip packages.
# We strip gpu4pyscf-* and cupy-cuda12x from environment.yml on the fly,
# preserving everything else (including pyscf which is the CPU-fallback path).
RUN sed -i \
    -e '/gpu4pyscf-cuda12x/d' \
    -e '/gpu4pyscf-libxc-cuda12x/d' \
    -e '/cupy-cuda12x/d' \
    /tmp/environment.yml && \
    conda env create -f /tmp/environment.yml && \
    conda clean -afy

# Activate qmmm by default
SHELL ["conda", "run", "-n", "qmmm", "/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/qmmm/bin:$PATH
ENV CONDA_DEFAULT_ENV=qmmm

# Smoke check at build time (fail-fast on env breakage)
RUN python -c "import openmm; import mdtraj; import pyscf; import rdkit; print('UPDD CPU env smoke OK')"

# Default command: run pytest with --tb=no -q (CI-friendly summary output).
# Override with `docker run ... <cmd>` for interactive use.
CMD ["python", "-m", "pytest", "tests/", "--tb=no", "-q"]
