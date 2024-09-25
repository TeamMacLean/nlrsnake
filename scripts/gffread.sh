#!/bin/bash

 source package gffread-0.12.7

input_fasta=$1
output_fasta=$2
input_gff=$3

gffread -g $input_fasta -y $output_fasta $input_gff