#!/usr/bin/bash

# This script will use the example files to run the HARE pipeline
# This requires the user to have 1. a genome reference and 2. a VEP cache
# Users who wish to use the default ~$HOME/.vep directory for cache can remove
# the --cache_dir flag from the first command

# Example command: bash example_run.sh {REFERENCE_PATH} {CACHE_DIRECTORY}

hare intersect --gwas example_gwas.txt --eoi example_eoi.bed --ref $1 --cache_dir $2 --out test_example -n 100
hare sigtest --input test_example.intersections --out test_example
hare prerank --input example_gwas.txt --out test_prerank
