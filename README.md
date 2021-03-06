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
   &nbsp;
   *Example of CHROMOSOME.list*
   ```bash
   chr1
   chr2
   chr3
   chr6u
   ```
   
4) Run vcfcaller workflow:
   ```bash
    nextflow run ./nf-variantcaller/variantcaller.nf --bams /PATH/TO/'Indivxxx.{bam,bam.bai}' --outdir PATH/TO/OUTDIR/
   ```
&nbsp;
&nbsp;
&nbsp;

# Example plots 
&nbsp;

![plot](./example_plots/Indivxxx_depth.png )

&nbsp;
&nbsp;

![plot](./example_plots/Indivxxx_quality.png)

## HPC enviroment
Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
