#!/bin/bash

## runs NLRtracker snakemake pipeline

function usage {
  cat <<EOM

usage: $(basename "$0") [OPTION]...

  -h        Display help
  dryrun    Do a dry run of the pipeline and get a report on what steps need to be done.
  unlock    Remove the lockfile preventing unauthorised run after failure of process.
  dag       Generate a graphviz dot file of the process DAG

Requires a file 'config.yaml' with three entries, e.g

 scratch: "/tsl/scratch/user_name/nlrsnake/"
 sample_fasta_file: "/hpc-home/user_name/samples_to_fasta.csv"
 seqs_per_file: 6000


  * 'scratch' is the path of a temporary working directory for all steps of the pipeline. A subdirectory for each sample will be created.
  * 'sample_fasta_file' is the name of a csv file linking input sequences to sample names
  * 'seqs_per_file' is the number of sequences per chunk for `interproscan` 6000 is a good number

EOM
  exit 2
}

if [ -z "$1" ]
then
   sbatch -J nlrtrack \
    -o nlrtrack.log \
    --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack.snakefile all --cluster 'sbatch --partition={params.queue} -c {threads} --mem={params.mem} --constraint="intel" ' -j 40 --latency-wait 60" \
    --constraint="intel"
elif [ $1 = 'unlock' ]
then
    sbatch -J unlock \
        -o nlrtrack.log \
        --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack.snakefile --unlock" \
        --constraint="intel" \
        --partition="tsl-short" \
        --mem="16G"
elif [ $1 = "dryrun" ]
then
    sbatch -J dryrun \
    -o nlrtrack.log \
    --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack.snakefile -n" \
    --constraint="intel" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = "dag" ]
then
    sbatch -J dag \
    -o nlrtrack.log \
    --wrap="source snakemake-5.5.3; snakemake --dag -s scripts/nlrtrack.snakefile  > nlrtrack.dot" \
    --constraint="intel" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = '-h' ]
then
  usage
fi