#!/bin/bash

source package /tsl/software/testing/bin/helixer-0.3.3

input_fasta=$1
lineage=$2
species=$3
output_gff=$4
model_path=$5
subsequence_length=$6
additional_options="${@:7}"

Helixer.py --lineage $lineage \
  --fasta-path $input_fasta \
  --species $species \
  --gff-output-path $output_gff \
  --model-filepath $model_path \
  --subsequence-length $subsequence_length \
  $additional_options
