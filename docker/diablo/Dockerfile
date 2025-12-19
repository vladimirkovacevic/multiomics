# DIABLO Analysis - Multi-Omics Integration using mixOmics
# Docker image for Data Integration Analysis for Biomarker Discovery using
# Latent Variable Approaches for Omics Studies

FROM rocker/r-ver:4.3.0

LABEL maintainer="Vladimir Kovacevic <vladimirkovacevic@example.com>"
LABEL description="DIABLO multi-omics integration with mixOmics package"
LABEL version="1.0"

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /usr/local/lib/R/etc/Rprofile.site

# Install BiocManager first
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"

# Install required R packages
RUN R -e "install.packages(c('tidyverse', 'ggplot2', 'RColorBrewer', 'pheatmap', 'gridExtra', 'reshape2', 'ellipse', 'igraph'), dependencies=TRUE)"

# Install mixOmics from Bioconductor
RUN R -e "BiocManager::install('mixOmics', dependencies=TRUE, update=FALSE)"

# Create scripts directory
RUN mkdir -p /scripts

# Copy R analysis script
COPY diablo_analysis.R /scripts/

# Set working directory
WORKDIR /data

# Create welcome message
RUN echo '#!/bin/bash\n\
echo ""\n\
echo "========================================================================"\n\
echo "  DIABLO Multi-Omics Integration with mixOmics"\n\
echo "========================================================================"\n\
echo ""\n\
echo "This container provides DIABLO analysis for multi-omics data integration."\n\
echo ""\n\
echo "Usage:"\n\
echo "  docker run --rm -v /path/to/data:/data vladimirkovacevic/diablo:latest \\\\\n\
echo "    Rscript /scripts/diablo_analysis.R \\\\\n\
echo "      expression.csv genotypes.csv covariates.csv output_prefix [n_comp] [keepX_expr] [keepX_geno]"\n\
echo ""\n\
echo "mixOmics version: $(R --slave -e \"cat(as.character(packageVersion('\'mixOmics\'')))\")" \n\
echo ""\n\
echo "========================================================================"\n\
echo ""\n' > /usr/local/bin/welcome.sh && chmod +x /usr/local/bin/welcome.sh

# Set entrypoint to show welcome message if no command provided
CMD ["/usr/local/bin/welcome.sh"]
