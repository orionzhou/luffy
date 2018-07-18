# Pipeline to annotate gene functions by InterProScan and Gene Ontology

This pipeline assumes you already have a multi-fasta file containing the 
representative (longest?) peptide sequences you are interested in.

In the current working directory create a sub-directiry and move the fasta
file there:  
```bash
mkdir 01.seqs
cd 01.seqs` 
ln -sf ../../52.rep.fas pro.fas
```

Split the fasta file into 50 pieces (so that jobs can be executed in parallel 
in servers):
```bash
pyfasta split -n 50 pro.fas
cd ..
mkdir 05.out
```

Create a job script (e.g., job) with the following content: 
```bash
#PBS -l nodes=1:ppn=8,walltime=10:00:00
#PBS -m ae
#PBS -M zhoux379@umn.edu
#PBS -q small

module load interproscan
cd $genome/Zmays_v4/61.interpro
let i=${PBS_ARRAYID}-1
printf -v pre 'pro.%02d' $i
interproscan.sh -i 01.seqs/${pre}.fas \
  -appl TIGRFAM,SFLD,ProDom,SMART,ProSiteProfiles,SUPERFAMILY,PRINTS,PANTHER,Pfam,Coils,MobiDBLite \
  -f tsv --goterms -o 05.out/${pre}.tsv
```

Submit the job by:
> qsub -t 1-50 job

Concatenate the interproscan outputs:
> cat 05.out/* > 06.ipr

Convert to tabular format:
> ipr2tsv.py 06.ipr 07.tsv

Use R script to generate `10_mrna.tsv` and `11_gene.tsv`

Convert to gene2go file:
> go.py 11_gene.tsv 15.tsv
