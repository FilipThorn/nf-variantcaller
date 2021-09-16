# nf-variantcaller
Nextflow pipeline that calls variants from reference in parallel per chromsome with bam and bam.bai as input files using freebayes. Merge chromosomal VCF to individual VCF. Produce summary plots to guide parameterisation for later VCF filtering using 10 % of variants randomly sampled from merged VCF.   

## Workflow

1) Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04) \
   Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10) 

2) Download git clone of this repository:
   ```bash
   git clone https://github.com/FilipThorn/nf-variantcaller
   ```
3) Edit nextflow.config file:
   ```bash
    ref = "/PATH/TO/INDEXED/REFERENCE.fa"
    chr = "/PATH/TO/CHROMOSOME.list"
   ```
4) Run vcfcaller workflow:
   ```bash
    nextflow run variantcaller.nf --bams /PATH/TO/'*.{bam,bam.bai}' --outdir PATH/TO/OUTDIR/
   ```
&nbsp;
&nbsp;
&nbsp;

#Plot examples 
![plot](./example_plots/Indivxxx_depth.png | width=100)

![plot](./example_plots/Indivxxx_quality.png | width=100)

## HPC enviroment
Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
