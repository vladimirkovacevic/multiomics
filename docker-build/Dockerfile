# eQTL Analysis Tutorial - Python/Jupyter Notebook Version
# Based on the original R tutorial from https://jknightlab.github.io/eqtl-intro/
#
# This Docker image provides a complete environment for learning eQTL analysis
# using Python, pandas, and statistical libraries.

FROM python:3.11-slim

LABEL maintainer="Vladimir Kovacevic <vladimirkovacevic@example.com>"
LABEL description="eQTL Analysis Tutorial with Python and Jupyter Notebook"
LABEL version="1.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /home/eqtl

# Copy requirements and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Create directory structure
RUN mkdir -p /home/eqtl/data/simulated \
    && mkdir -p /home/eqtl/notebooks

# Copy notebook (data will be mounted at runtime)
COPY notebooks/eQTL_tutorial.ipynb /home/eqtl/notebooks/

# Create a welcome message
RUN echo '#!/bin/bash\n\
echo ""\n\
echo "========================================================================"\n\
echo "  eQTL Analysis Tutorial - Python/Jupyter Notebook"\n\
echo "========================================================================"\n\
echo ""\n\
echo "Welcome to the eQTL Analysis Tutorial!"\n\
echo ""\n\
echo "This tutorial teaches expression Quantitative Trait Locus (eQTL) analysis"\n\
echo "using Python, covering:"\n\
echo ""\n\
echo "  1. Simple SNP-Gene Expression Associations"\n\
echo "  2. Handling Confounding Variation"\n\
echo "  3. Principal Components Analysis"\n\
echo "  4. Genome-Wide eQTL Analysis"\n\
echo ""\n\
echo "------------------------------------------------------------------------"\n\
echo "Quick Start:"\n\
echo "------------------------------------------------------------------------"\n\
echo ""\n\
echo "Your Jupyter Notebook server is starting..."\n\
echo ""\n\
echo "Access the notebook at: http://localhost:8887"\n\
echo ""\n\
echo "Open: notebooks/eQTL_tutorial.ipynb"\n\
echo ""\n\
echo "Data location: /home/eqtl/data/simulated/"\n\
echo ""\n\
echo "========================================================================"\n\
echo ""\n' > /home/eqtl/welcome.sh && chmod +x /home/eqtl/welcome.sh

# Expose Jupyter port
EXPOSE 8887

# Set environment variables
ENV JUPYTER_ENABLE_LAB=yes

# Start Jupyter Notebook
CMD ["/bin/bash", "-c", "/home/eqtl/welcome.sh && jupyter notebook --ip=0.0.0.0 --port=8887 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password='' --notebook-dir=/home/eqtl"]
