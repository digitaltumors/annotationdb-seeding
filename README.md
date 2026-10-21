# AnnotationDB Seeding

Repository for building and populating tables in the AnnotationDB MySQL database. Specifically, this database coordinates the creation and population of a drug database for PubChem, a ChEMBL drug database with FKs to pubchem CIDs, and a cell lines database using cellosaurus. This repository utilizes [Pixi](https://pixi.sh/latest/) for package management and [SQLAlchemy](https://www.sqlalchemy.org/) ORM for main database operations.

## Environment Variables
```Bash
DATABASE_IP=
PORT=
DATABASE_USER=
DATABASE_PASS=
SELECTED_DB=
```

## Compound Data Extraction Pipeline

The compound data extraction pipeline will take a CSV of compounds and map them to their respective Pubchem (properties, synonyms, bioassays, toxicity) and ChEMBL (mechanisms) metadata. The input file for compounds is to be structured as a csv with a column "input_id" for preserving the identity of the compound and another column that can take many different forms such as "cid", "pubchem", "inchikey", "smiles", depending on what your dataset is providing for the compounds in questions.

The compound data extraction pipeline can be run in three different modes:

### 1. Local Sequential (Rscript)
Best for small to medium lists of compounds (< 5,000) or local testing. It runs locally without Docker and synchronizes directly with Google Cloud Storage (GCS).

```bash
cd data_extraction/drugs/pubchem/kubernetes

Rscript input_id_coordinator.R \
  --in_csv "input/new_input_compounds.csv" \
  --gcs_bucket "gs://your_bucket_name_here" \
  --batch_size 50
```
This script handles its own chunking, automatically uploads intermediate batch results to GCS, and merges them into the master `gs://.../output/` files upon completion. It will automatically resume from previous runs.

### 2. Local Docker
Useful for validating the exact container environment that will be used in GKE, or for isolating dependencies on your local machine.

```bash
cd data_extraction/drugs/pubchem/kubernetes

# 1. Generate the batches in GCS
python coordinator.py \
  --in_csv "gs://annotationdb_data_retrieval/input/new_input_compounds.csv" \
  --batch-size 5000 \
  --max-batches 10

# 2. Build the Docker image
docker build -t annotationdb-scraper .

# 3. Run 10 containers in the background (simulating a Kubernetes Job for batches 0 through 9)
# Note: The -v flag mounts your local Google Cloud credentials into the container so that `gsutil` inside the container can authenticate and access the GCS bucket.
for i in {0..9}; do
  docker run -d -e BATCH_INDEX=$i \
    -v ~/.config/gcloud:/root/.config/gcloud \
    annotationdb-scraper
done

# 4. Merge the batch output back into the master files
python merge_outputs.py
```
*Note: Like GKE execution, running the container will only upload its specific batch to the `output_batches/` directory in GCS. You must manually run `merge_outputs.py` afterwards to consolidate those batches into the main `output/` directory.*

### 3. Cloud Execution (Google Kubernetes Engine)
Best for massive lists (100k+ compounds). This orchestrates dozens of Pods in parallel.

```bash
# Create cluster on GKE

gcloud container clusters create annotationdb-extraction \
  --zone northamerica-northeast2-a \
  --num-nodes 4 \
  --machine-type e2-standard-2 \
  --network default \
  --subnetwork default \
  --service-account annotationdb-scraping@annotationdb.iam.gserviceaccount.com


gcloud container clusters get-credentials annotationdb-extraction --location=northamerica-northeast2-a
```

```bash
cd data_extraction/drugs/pubchem/kubernetes

# 1. Generate the batches and the K8s Job YAML
python coordinator.py \
  --in_csv "gs://annotationdb_data_retrieval/input/new_input_compounds.csv" \
  --batch-size 5000 \
  --max-batches 10

# 2. Build and Push the Docker image (required if you changed the code/dependencies)
docker buildx build --platform linux/amd64 -t bhklabmattbocc/annotationdb-r:v9 --push .

# 3. Deploy the Job to GKE
kubectl apply -f job.yaml

# 3. alternative: deploy job locally

docker run --rm \
  -v ~/.config/gcloud:/root/.config/gcloud \
  bhklabmattbocc/annotationdb-r:v14 0


# 4. After the Job completes, merge all parallel batch outputs into the master files
python merge_outputs.py


# 5. Delete job once it has been completed
kubectl delete -f job.yaml

```

```bash

# Tracking all pod status
kubectl get pods -l app=annotationdb-extract -w

# Tracking a single pod progress
kubectl logs -f {pod_id}

# Tracking job progress across all pods
kubectl logs -f -l app=annotationdb-extract

```

## Database Seeding

Once the data is extracted and available in your master `output/` directory, you can push it to your MySQL database. If the latest files are in GCS, download them locally first, or just run the local extraction pipeline which caches them.

### Appending New Compounds Safely (Recommended)
To incrementally add new compounds without modifying or duplicating existing records, use the safe append script. This utilizes `INSERT IGNORE` and strictly respects all Foreign Key insert order requirements.

```bash
# 1. Dry-run to see what will be inserted without committing to the DB
pixi run python seeding/append_new_compounds.py \
  --data_dir data_extraction/drugs/pubchem/kubernetes/output \
  --dry-run

# 2. Execute the actual append
pixi run python seeding/append_new_compounds.py \
  --data_dir data_extraction/drugs/pubchem/kubernetes/output
```

### Initializing the DB (First-Time Setup)
To create the database tables from scratch for the first time, use the `seeding_coordinator.py` script. **Note: This script does NOT wipe existing data.** It runs a `create_all()` to build missing tables and attempts to naive-insert all rows. If the database already has data, it will crash on Primary Key conflicts, which is why `append_new_compounds.py` is recommended for subsequent updates.

```bash
pixi run python seeding/seeding_coordinator.py \
  --data_dir data_extraction/drugs/pubchem/kubernetes/output
```
