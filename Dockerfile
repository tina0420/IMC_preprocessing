FROM python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    libhdf5-dev \
    libatlas-base-dev \
    liblapack-dev \
    libblas-dev \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    scanpy==1.10.1 \
    anndata==0.10.5 \
    numpy==1.26.4 \
    scipy==1.12.0 \
    pandas==2.2.1 \
    matplotlib==3.8.0 \
    scikit-learn==1.4.0 \
    seaborn \
    tqdm



