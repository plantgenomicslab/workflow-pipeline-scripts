#!/bin/bash

snakemake --cluster-config cluster.json \
          --cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} --parsable -A {cluster.account} -p {cluster.partition} -t {cluster.time} -o {cluster.output} -e {cluster.error}" \
          --cluster-status ./scripts/status.py \
          --use-conda \
          --use-singularity \
          --cores 8 \
          "$@"

