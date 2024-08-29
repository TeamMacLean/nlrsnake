'''Run NLRtracker pipeline on TSL HPC

Usage: bash scripts/do_nlrtrack.sh

'''



configfile: "config.yaml"

def get_fa_samples(file):
    print(file)
    samples = []
    fa = []
    with open(file, 'r') as input:
        for line in input:
            line = line.rstrip()
            items = line.split(",")
            samples.append(items[0])
            fa.append(items[1])
    return fa, samples

def sample_to_fa(sample, fastas, samples):
    return fastas[samples.index(sample)]

FASTAS, SAMPLES = get_fa_samples(config['sample_fasta_file'])

rule all:
    input: "done.txt"

rule strip_asterisks:
    input: lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES)
    output: config["scratch"] + "{sample}/no_asterisk.fa"
    params:
        mem="4G",
        queue="tsl-short"
    threads: 1
    shell: "cat {input} | sed 's/*//g' > {output}"


checkpoint interpro:
    input:
        config["scratch"] + "{sample}/no_asterisk.fa"
        #config["input_fasta"]
    output:
        fa=directory(config["scratch"] + "{sample}/ipro_fa/")
    params:
        mem="4G",
        queue="tsl-short",
        seqs_per_file=config['seqs_per_file']
    threads: 1
    shell: "bash scripts/split_fa.sh {input} {params.seqs_per_file} {output}"

rule run_interpro:
    input:
        config['scratch'] + "{sample}/ipro_fa/{n}.fa"
    output:
         config['scratch'] + "{sample}/ipro_gff/{n}.gff"
    params:
        tmp=config['scratch'] + "{sample}/ipro_tmp",
        mem="16G",
        queue="tsl-short",
    threads: 16
    #shell: "cp {input} {output}"
    shell: "bash scripts/interpro.sh {input} {params.tmp} {threads} {output}"

rule fimo:
    input: lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES)
    output: "{sample}/results/fimo_out/fimo.gff",
    params:
        xml="lib/meme.xml",
        mem="16G",
        queue="tsl-short",
        fimo_outdir="{sample}/results/fimo_out"
    threads: 8
    shell: "bash scripts/fimo.sh {params.fimo_outdir} {params.xml} {input}"
        #"source nlrtracker-1.0.0; "
        #"fimo -o {params.fimo_outdir} {params.xml} {input};"

rule hmmer:
    input: lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES)
    output: "{sample}/results/CJID.txt"
    params:
        hmm="lib/abe3069_Data_S1.hmm",
        mem="16G",
        queue="tsl-short"
    threads: 8
    shell: "bash scripts/hmmer.sh {output} {params.hmm} {input}"
        #"source nlrtracker-1.0.0; "
        #"hmmsearch --domtblout {output} {params.hmm} {input};"

rule nlrtracker:
    input:
        interpro="{sample}/results/interpro_result.gff",
        fimo="{sample}/results/fimo_out/fimo.gff",
        fa=lambda wildcards: sample_to_fa(wildcards.sample, FASTAS, SAMPLES),
        hmmer="{sample}/results/CJID.txt"
    output: "{sample}/nlrtracker_out/done.txt"
    params:
        itol="lib/iTOL_NLR_template.txt",
        ipro_list="lib/interproscan_v5.69-101.0_entry.list",
        run_dir="{sample}/nlrtracker_out",
        mem="16G",
        queue="tsl-short",
    threads: 8
    shell: "bash scripts/run_tracker.sh {params.ipro_list} {input.interpro} {input.fimo} {input.fa} {params.run_dir} {wildcards.sample} p {input.hmmer} {params.itol}"
        #"source nlrtracker-1.0.0; "
        #"Rscript scripts/NLRtracker.R {params.ipro_list} {input.interpro} {input.fimo} {input.fa} {output} p {input.hmmer} {params.itol};"


def aggregate_ipro(wildcards):
    checkpoint_output = checkpoints.interpro.get(**wildcards).output[0]
    return expand(config['scratch'] + "{sample}/ipro_gff/{n}.gff", sample=wildcards.sample, n=glob_wildcards(os.path.join(checkpoint_output, "{n}.fa")).n)

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


rule gather_dones:
    input: expand("{sample}/nlrtracker_out/done.txt", sample=SAMPLES)
    output: "done.txt"
    params:
        mem="8G",
        queue="tsl-short",
    threads: 1
    shell:
        "cat {input} > {output}"
