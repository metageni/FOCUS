# Docker / Singularity — focus

Two Dockerfiles are provided:

| File | Base | Best for |
|---|---|---|
| `Dockerfile` (root) | `continuumio/miniconda3` | Local development, quick runs |
| `Docker/Dockerfile` | `ubuntu:24.04` + miniforge | HPC, Singularity, `linux/amd64` |

---

## Local Docker (root Dockerfile)

### Build

```bash
cd focus/
docker build -t focus .
```

### Run

Mount your input and output directories into the container:

```bash
docker run --rm \
  -v $(pwd)/samples:/data/samples \
  -v $(pwd)/results:/data/results \
  focus -q /data/samples -o /data/results
```

With k-mer size 7 and 8 threads:

```bash
docker run --rm \
  -v $(pwd)/samples:/data/samples \
  -v $(pwd)/results:/data/results \
  focus -q /data/samples -o /data/results -k 7 -t 8
```

---

## HPC / Singularity (Docker/Dockerfile)

### Build and push

```bash
cd focus/Docker/
docker build -t focus_hpc .

# tag and push to your registry
TAG=$(docker images focus_hpc --format "{{.ID}}" | head -1)
docker tag focus_hpc $USER/focus:v2.0.0_${TAG}
docker push $USER/focus:v2.0.0_${TAG}
```

### Pull as a Singularity image

```bash
export TAG=<your_tag>
module load singularity
mkdir -p sif tmp
export SINGULARITY_TMPDIR=$PWD/tmp/
singularity pull --dir sif docker://$USER/focus:v2.0.0_${TAG}
```

### Run with Singularity

```bash
singularity exec \
  --bind /path/to/project:/data \
  sif/focus_v2.0.0_${TAG}.sif \
  focus -q /data/samples -o /data/results -k 7
```

### Example SLURM job

```bash
#!/bin/bash
#SBATCH --job-name=focus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

module load singularity

singularity exec \
  --bind /path/to/project:/data \
  sif/focus_v2.0.0_${TAG}.sif \
  focus -q /data/samples -o /data/results -k 7 -t 8
```

---

## Cleaning up Docker

```bash
# remove all stopped containers and dangling images
docker system prune -f

# remove a specific image
docker rmi focus
```
