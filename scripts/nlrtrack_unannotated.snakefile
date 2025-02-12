import os
import csv
from snakemake.io import glob_wildcards

configfile: "config.yaml"

# Define the final output
FINAL_OUTPUT = "done.txt"

# Function to get samples and fasta paths
def get_fa_samples(file):
    samples = []
    fa = []
    with open(file, 'r') as input:
        for line in input:
            line = line.rstrip()
            items = line.split(",")
            # Adding check to ensure samples_to_Fasta file is formatted correctly
            if len(items) < 2:
                raise ValueError(f"Malformed line in {file}: {line}")
            samples.append(items[0])
            fa.append(items[1])
    return fa, samples

# Function to map sample to fasta path
def sample_to_fa(sample, fastas, samples):
    return fastas[samples.index(sample)]

# Function to read the new CSV generated by "create_new_samples"
def updated_samples(wildcards):
    # Get the output from the checkpoint "create_new_samples" which is required after running rule gffread
    checkpoint_output = checkpoints.create_new_samples.get(**wildcards).output.new_samples
    # Read the updated sample list
    return get_fa_samples(checkpoint_output)

FASTAS, SAMPLES = get_fa_samples(config['sample_fasta_file'])

# Main rule to ensure correct order of execution
rule all:
    input: FINAL_OUTPUT

# Rule to run Helixer
rule run_helixer:
    input:
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES)
    output:
        helixer=config['scratch'] + "/{sample}/helixer/{sample}_helixer.gff"
    params:
        lineage="land_plant",
        species=config['species'],
        model_path="/tsl/data/helixer/models/land_plant/land_plant_v0.3_a_0080.h5",
        subsequence_length=64152,
        additional_options=config.get('helixer_options', ''),
        mem="32G",
        queue="tsl-medium"
    threads: 16
    shell:
        """
        bash scripts/helixer.sh {input.fasta} {params.lineage} {params.species} \
            {output.helixer} {params.model_path} {params.subsequence_length} {params.additional_options}
        """

# Rule to run gffread
rule run_gffread:
    input:
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES),
        helixer=config['scratch'] + "/{sample}/helixer/{sample}_helixer.gff"
    output:
        gffread=config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta"
    params:
        mem="32G",
        queue="tsl-short",
        temp_dir=config['scratch'] + "/{sample}/gffread"
    threads: 16
    run:
        
        # Create temp directory if it doesn't exist
        shell("mkdir -p {params.temp_dir}")
        
        # Create a temporary uncompressed file
        temp_fasta = f"{params.temp_dir}/temp.fasta"
        
        # Check if input fasta is compressed without using Snakemake wildcard
        shell("""
            if [[ {input.fasta} == *.gz ]]; then 
                gunzip -c {input.fasta} > {params.temp_dir}/temp.fasta
            else
                cp {input.fasta} {params.temp_dir}/temp.fasta
        fi
        """)
        
        # Run gffread
        shell("bash scripts/gffread.sh {temp_fasta} {output.gffread} {input.helixer}")
        
        # Remove temporary file
        shell("rm {temp_fasta}")


# Checkpoint to create new samples file and ensure pipeline waits until complete
checkpoint create_new_samples:
    input:
        gffread=expand(config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta", sample=SAMPLES)
    output:
        new_samples=config['scratch'] + "/new_samples.csv"
    params:
        mem="32G",
        queue="tsl-short"
    threads: 16
    run:
        with open(output.new_samples, 'w') as outfile:
            for sample, gffread_file in zip(SAMPLES, input.gffread):
                outfile.write(f"{sample},{os.path.realpath(gffread_file)}\n")

# Assign new sample paths from files generated by gffread
rule assign_new_samples:
    input:
        new_csv=config['scratch'] + "/new_samples.csv"
    output:
        update_samples=config['scratch'] + "/{sample}_update_done.txt"
    params:
        mem="4G",
        queue="tsl-short"
    threads: 4
    run:
        # Load the updated samples
        FASTAS, SAMPLES = updated_samples({})
        shell("touch {output.update_samples}")

# Rule to strip asterisks
rule strip_asterisks:
    input:
        fasta=config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta"
    output:
        no_asterisk=config['scratch'] + "/{sample}/no_asterisk.fa"
    params:
        mem="4G",
        queue="tsl-short"
    shell:
        "cat {input.fasta} | sed 's/*//g' > {output.no_asterisk}"

# Checkpoint for interpro
checkpoint interpro:
    input:
        config['scratch'] + "/{sample}/no_asterisk.fa"
    output:
        fa=directory(config['scratch'] + "/{sample}/ipro_fa/")
    params:
        mem="4G",
        queue="tsl-short",
        seqs_per_file=config['seqs_per_file']
    threads: 1
    shell: "bash scripts/split_fa.sh {input} {params.seqs_per_file} {output}"

# Rule to run interpro on split files
rule run_interpro:
    input:
        fa=config['scratch'] + "/{sample}/ipro_fa/{n}.fa"
    output:
        gff=config['scratch'] + "/{sample}/ipro_gff/{n}.gff"
    params:
        mem="16G",
        queue="tsl-short",
        tmp=config['scratch'] + "/{sample}/ipro_tmp"
    threads: 16
    shell:
        "bash scripts/interpro.sh {input.fa} {params.tmp} {threads} {output.gff}"

# Rule to run FIMO
rule run_fimo:
    input:
        # Update to obtain new sample files
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, *updated_samples(wildcards)),
        update_samples=config['scratch'] + "/{sample}/{sample}_update_done.txt"
    output:
        fimo="{sample}/results/fimo_out/fimo.gff"
    params:
        mem="32G",
        queue="tsl-short"
    threads: 16
    shell:
    	"bash scripts/fimo.sh {wildcards.sample}/results/fimo_out lib/meme.xml {input.fasta}"

# Rule to run HMMER
rule run_hmmer:
    input:
        # Update to obtain new sample paths
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, *updated_samples(wildcards)),
        update_samples=config['scratch'] + "/{sample}/{sample}_update_done.txt"
    output:
        hmmer="{sample}/results/CJID.txt"
    params:
        mem="32G",
        queue="tsl-short"
    threads: 16
    shell:
        "bash scripts/hmmer.sh {output.hmmer} lib/abe3069_Data_S1.hmm {input.fasta}"

def aggregate_ipro(wildcards):
    checkpoint_output = checkpoints.interpro.get(**wildcards).output[0]
    return expand(config['scratch'] + "/{sample}/ipro_gff/{n}.gff", sample=wildcards.sample, n=glob_wildcards(os.path.join(checkpoint_output, "{n}.fa")).n)

rule aggregate:
    input:
        aggregate_ipro
    output:
        "{sample}/results/interpro_result.gff"
    params:
        mem="8G",
        queue="tsl-short",
    threads: 1
    shell:
        "cat {input} > {output}"

# Rule to run NLRtracker
rule run_nlrtracker:
    input:
        interpro="{sample}/results/interpro_result.gff",
        fimo="{sample}/results/fimo_out/fimo.gff",
        hmmer="{sample}/results/CJID.txt",
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES)
    output:
        nlrtracker="{sample}/nlrtracker_out/done.txt"
    params:
        mem="32G",
        queue="tsl-short"
    threads: 16
    shell:
        """
        bash scripts/run_tracker.sh lib/interproscan_v5.69-101.0_entry.list {input.interpro} {input.fimo} \
        {input.fasta} {wildcards.sample}/nlrtracker_out {wildcards.sample} p {input.hmmer} lib/iTOL_NLR_template.txt
        """

# Rule to finalize the pipeline
rule finalize:
    input:
        # Changed to delay expansion
        lambda wildcards: expand("{sample}/nlrtracker_out/done.txt", sample=SAMPLES)
    output:
        touch(FINAL_OUTPUT)
    params:
        mem="4G",
        queue="tsl-short"
    threads: 1
    shell:
        "echo 'Pipeline completed successfully' > {output}"
