FROM rocker/r-ver:4.5.1

# Install system dependencies + Python 3
RUN apt-get update && apt-get install -y \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    curl \
    python3 \
    python3-pip \
    apt-transport-https \
    ca-certificates \
    gnupg \
    && rm -rf /var/lib/apt/lists/*

# Install Google Cloud SDK (for gsutil)
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update && apt-get install -y google-cloud-cli && \
    rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e "install.packages(c('data.table', 'httr', 'jsonlite', 'remotes'), repos='https://cloud.r-project.org')"

# Install AnnotationGx from GitHub
RUN R -e "remotes::install_github('bhklab/AnnotationGx', dependencies=TRUE, build=FALSE, build_manual=FALSE, build_vignettes=FALSE)"

# Set up working directory
WORKDIR /app

# Copy R scripts
COPY pubchem_drug_extraction_from_input_ids.R .
COPY assay_extraction_input_ids.R .
COPY toxicity_extraction_input_ids.R .

# Copy entrypoint
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]