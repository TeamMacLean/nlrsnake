
#source nlrtracker-1.0.0
mkdir -p $5
Rscript scripts/NLRtracker.R \
  $1 \
  $2 \
  $3 \
  $4 \
  $5 \
  $6 \
  $7 \
  $8 \
  $9

date > $5/done.txt
#vars =
 # {params.ipro_list} {input.interpro} {input.fimo} {input.fa} {params.rundir} p {input.hmmer} {params.itol}
