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
            samples.append(items[0])
            fa.append(items[1])
    return fa, samples

# Function to map sample to fasta path
def sample_to_fa(sample, fastas, samples):
    return fastas[samples.index(sample)]

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
        queue="tsl-short"
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
        queue="tsl-short"
    threads: 16
    shell:
        """
        bash scripts/gffread.sh {input.fasta} {output.gffread} {input.helixer}
        """

# Rule to update samples file
rule update_samples:
    input:
       fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES),
       gffread=config['scratch'] + "/{sample}/gffread/{sample}_gffread.fasta"
    output:
        update_samples=config['scratch'] + "/{sample}/{sample}_update_done.txt"
    params:
        mem="4G",
        queue="tsl-short"
    threads: 4
    run:
        samples_file = config['scratch'] + "/samples_to_fasta.csv"
        temp_file = config['scratch'] + "/temp_samples_to_fasta.csv"

        with open(samples_file, 'r') as in_file, open(temp_file, 'w', newline='') as out_file:
            reader = csv.reader(in_file)
            writer = csv.writer(out_file)

            for row in reader:
                if row[0] == wildcards.sample:
                    new_path = os.path.realpath(config['scratch'] + f'/{wildcards.sample}/gffread/{wildcards.sample}_gffread.fasta')
                    writer.writerow([row[0], new_path])
                else:
                    writer.writerow(row)

        os.replace(temp_file, samples_file)
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
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES),
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
        fasta=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES),
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
        expand("{sample}/nlrtracker_out/done.txt", sample=SAMPLES)
    output:
        touch(FINAL_OUTPUT)
    params:
        mem="4G",
        queue="tsl-short"
    threads: 1
    shell:
        "echo 'Pipeline completed successfully' > {output}"