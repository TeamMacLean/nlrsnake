# nlrsnake - running NLRtracker.R on the HPC

This repo contains files that will run `NLRtracker.R` and its dependencies on the TSL HPC as a snakemake
pipeline. It replaces (so never runs) the previous wrapper `NLRtracker.sh`. 

## Why snakemake?

`Snakemake` is a workflow manager. Crucially, it can restart where it left off, so if one step of the pipeline fails, not all 
previous steps need to be re-run. It manages the steps that need to be run for you. This makes it very 
useful for bioinformatics pipelines especially this one, where we want to divide some steps into many similar small ones.

## The pipeline

The pipeline takes multiple sets of protein sequences then in parallel, for each set:

    1. Runs `interprocan.sh` 
    2. Runs `hmmer`
    3. Runs `fimo`    
    
Once these are done, `NLRtracker.R` is run and a results folder for each sample created

    
### Interproscan step 

 The pipeline begins by removing asterisks, splitting the protein file into chunks of specified size and runnning each chunk as  
 as a separate job on the HPC. Once all the jobs are complete the results are compiled into one large output.
 
### Other steps

All other steps are run individually. Steps run as soon as they are able (that is, the files they rely on are prepared). 
So `fimo` and `hmmer` run alongside `interproscan` jobs, but `NLRtracker.R` needs to wait until all its inputs are ready.
Snakemake manages this for you.

## Running the pipeline

The pipeline is just a collection of scripts, so just pull this repo to a folder in your home area and `cd` into
that folder.

### `config.yaml`

The file `config.yaml` contains job info and needs to be filled in. The pipeline looks for it automatically.
It is brief

```
 scratch: "/tsl/scratch/user_name/nlrsnake/"
 sample_fasta_file: "/hpc-home/user_name/samples_to_fasta.csv"
 seqs_per_file: 6000
 species: "Arabidopsis_thaliana"
 helixer_options: ""
```
    * `scratch` is the path of a temporary working directory for all steps of the pipeline. A subdirectory for each sample will be created.
    * `sample_fasta_file` is the name of a csv file linking input sequences to sample names
    * `seqs_per_file` is the number of sequences per chunk for `interproscan` 6000 is a good number
    * 'species' is the species of the organism whos genome is to be annotated.  Should have an underscore (_) between genus and species names.
    * 'helixer_options' provides the opportunity to add additional usage options beyond the minimum.  Note default lineage is "land_plant".
    
Try to use the absolute (full length) path in the `config.yaml` file

### `sample_fasta_file`

The format is `sample_name,path` , one per line. No blank lines at the end of the file. No spaces. Alphanumeric characters and underscores only. 

```
    sample_1,/hpc-home/user_name/seqs1.fa
```
    
### The submission script

Submitting the snakemake pipeline to the cluster is done from the script `scripts/do_nlrtrack.sh`.

    1. See the help - `bash scripts/do_nlrtrack.sh -h`
    2. Do a dryrun (check what remains of the pipeline to run and check everything is in place) - `bash scripts/do_nlrtrack.sh dryrun`
    3. Submit the pipeline - `bash scripts/do_nlrtrack.sh`
    
## Checking the pipeline status once submitted

The pipeline creates a master job called `nlrtrack` and that creates subjobs, this will persist for as long as the pipeline runs. Output from the main process goes
into a file `nlrtrack.log`. Other job output goes into job id specific `slurm**.out` files.  


## Output folder

In this pipeline all `NLRtracker.R` output goes into the folder `nlrtracker_out`. Other folders (`results`, `ipro_fa` and `ipro_gff`) are 
interim data files and can be discarded. 


## Restarting the pipeline

If the pipeline fails (usually because some file wasn't created somewhere), once a fix is made, the pipeline can be restarted from the point it failed. Just redo the submission step
`bash scripts/do_nlrtrack.sh` . The pipeline will go from as far as it got previously, so if you need to redo earlier steps, you need to delete their output files.



    
