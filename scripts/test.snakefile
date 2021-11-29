rule all:
    input:
        "aggregated/a.txt",
        "aggregated/b.txt"


# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint clustering:
    input:
        "samples/{sample}.txt"
    output:
        clusters=directory("clustering/{sample}")
    shell:
        "mkdir clustering/{wildcards.sample}; "
        "for i in 1 2 3; do echo $i > clustering/{wildcards.sample}/$i.txt; done"


# an intermediate rule
rule intermediate:
    input:
        "clustering/{sample}/{i}.txt"
    output:
        "post/{sample}/{i}.txt"
    shell:
        "cp {input} {output}"


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
    return expand("post/{sample}/{i}.txt",
           sample=wildcards.sample,
           i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)


# an aggregation over all produced clusters
rule aggregate:
    input:
        aggregate_input
    output:
        "aggregated/{sample}.txt"
    shell:
        "cat {input} > {output}"