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
 species: "Arabidopsis_thaliana"
 helixer_options: "--subsequence-length 35000 --lineage fungi"


  * 'scratch' is the path of a temporary working directory for all steps of the pipeline. A subdirectory for each sample will be created.
  * 'sample_fasta_file' is the name of a csv file linking input sequences to sample names
  * 'seqs_per_file' is the number of sequences per chunk for `interproscan` 6000 is a good number.
  * 'species' is the species of the organism whos genome is to be annotated.  Should have an underscore (_) between genus and species names.
  * 'helixer_options' provides the opportunity to add additional usage options beyond the minimum.  Note default lineage is "land_plant"

EOM
  exit 2
}

if [ -z "$1" ]
then
   sbatch -J nlrtrack \
    -o nlrtrack.log \
    --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack_unannotated.snakefile all --cluster 'sbatch --partition={params.queue} -c {threads} --mem={params.mem} ' -j 40 --max-jobs-per-second 5 --max-status-checks-per-second 5 --latency-wait 120 --keep-going --rerun-incomplete" \
    --partition="tsl-long" \
    --mem="4G"
elif [ $1 = 'unlock' ]
then
    sbatch -J unlock \
        -o nlrtrack.log \
        --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack_unannotated.snakefile --unlock" \
        --partition="tsl-short" \
        --mem="16G"
elif [ $1 = "dryrun" ]
then
    sbatch -J dryrun \
    -o nlrtrack.log \
    --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack_unannotated.snakefile -n" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = "debug" ]
then
    sbatch -J debug \
    -o nlrtrack.log \
    # Edit below as required to debug
    --wrap="source snakemake-5.5.3; snakemake -s scripts/nlrtrack_unannotated.snakefile -n -p add_rule_to_debug" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = "dag" ]
then
    source package graphviz-2.38.0
    sbatch -J dag \
    -o nlrtrack.log \
    # Edit to have the command make a dag image directly
    --wrap="source snakemake-5.5.3; snakemake --dag -s scripts/nlrtrack_unannotated.snakefile | dot -Tsvg > nlrtrack.svg" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = '-h' ]
then
  usage
fi
