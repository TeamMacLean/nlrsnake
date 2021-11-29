'''Run NLRtracker pipeline on TSL HPC

Usage: bash scripts/do_nlrtrack.sh

'''



configfile: "config.yaml"


rule all:
    input: "nlrtracker_out/done.txt"

rule strip_asterisks:
    input: config["input_fasta"]
    output: config["ipro_tmp"] + "no_asterisk.fa"
    params:
        mem="4G",
        queue="tsl-short"
    threads: 1
    shell: "cat {input} | sed 's/*//g' > {output}"


checkpoint interpro:
    input:
        config["ipro_tmp"] + "no_asterisk.fa"
        #config["input_fasta"]
    output:
        fa=directory("ipro_fa/")
    params:
        mem="4G",
        queue="tsl-short",
        seqs_per_file=config['seqs_per_file']
    threads: 1
    shell: "bash scripts/split_fa.sh {input} {params.seqs_per_file} {output}"

rule run_interpro:
    input:
        "ipro_fa/{n}.fa"
    output:
        "ipro_gff/{n}.gff"
    params:
        tmp=config["ipro_tmp"],
        mem="16G",
        queue="tsl-short",
    threads: 8
    #shell: "cp {input} {output}"
    shell: "bash scripts/interpro.sh {input} {params.tmp} {threads} {output}"

rule fimo:
    input: config["input_fasta"]
    output: "results/fimo_out/fimo.gff",
    params:
        xml="lib/meme.xml",
        mem="16G",
        queue="tsl-short",
        fimo_outdir="results/fimo_out"
    threads: 8
    shell: "bash scripts/fimo.sh {params.fimo_outdir} {params.xml} {input}"
        #"source nlrtracker-1.0.0; "
        #"fimo -o {params.fimo_outdir} {params.xml} {input};"

rule hmmer:
    input: config['input_fasta']
    output: "results/CJID.txt"
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
        interpro="results/interpro_result.gff",
        fimo="results/fimo_out/fimo.gff",
        fa=config['input_fasta'],
        hmmer="results/CJID.txt"
    output: "nlrtracker_out/done.txt"
    params:
        itol="lib/iTOL_NLR_template.txt",
        ipro_list="lib/interproscan_5.51-85.0.list",
        run_dir="nlrtracker_out",
        mem="16G",
        queue="tsl-short",
    threads: 8
    shell: "bash scripts/run_tracker.sh {params.ipro_list} {input.interpro} {input.fimo} {input.fa} {params.run_dir} p {input.hmmer} {params.itol}"
        #"source nlrtracker-1.0.0; "
        #"Rscript scripts/NLRtracker.R {params.ipro_list} {input.interpro} {input.fimo} {input.fa} {output} p {input.hmmer} {params.itol};"


def aggregate_ipro(wildcards):
    checkpoint_output = checkpoints.interpro.get(**wildcards).output[0]
    return expand("ipro_gff/{n}.gff", n=glob_wildcards(os.path.join(checkpoint_output, "{n}.fa")).n)

rule aggregate:
    input:
        aggregate_ipro
    output:
        "results/interpro_result.gff"
    params:
        mem="8G",
        queue="tsl-short",
    threads: 1
    shell:
        "cat {input} > {output}"
