#!/usr/bin/env nextflow

if (params.help) {
    log.info """\
            __________________
            |                |
            | |```````| |`````
            | |____   | |
            |     |   | |
            | |````   | |
            | |ilip   | |hörn     
            –––––––––––––––––––––––––––––––––––––––
            VariantCaller
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            'USAGE'
            
            nextflow run variantcaller.nf --bams /PATH/TO/'*.{bam,bam.bai}' --outdir /PATH/TO/OUTDIR/


            'Mandatory arguments:'
            --outdir                 Path to output directory
            --bams                   Path to bamfile and indexed bamfile pairs
            
            'OPTIONS'
            --help                   Outputs this help log      
            -resume                  Nextflow cmd to resume modified workflow


            'HPC'
            -profile       FILE      If intention to run workflow on HPC please provide a suitable profile 
                                     in the nextflow.config file 


            For freebayes see https://github.com/freebayes/freebayes
            'SUPPORT'
            Email Filip.Thorn@NRM.se for questions on script
            """
    exit 1
}


log.info """\
         –––––––––––––––––––––––––––––––––––––––
         VariantCaller 
         NEXTFLOW   P I P E L I N E                
         –––––––––––––––––––––––––––––––––––––––
         outdir       : ${params.outdir}
         bams         : ${params.bams}
         ref          : ${params.ref}
         """
         .stripIndent()


// Channels for bam
Channel.fromFilePairs( params.bams, flat:true )
       .set { bam_ch } 

// make chromosome list
chr_list = file(params.chr).readLines()

process Variantcall {

    tag "${chr};${sample_id}"

    publishDir "${params.outdir}/$sample_id", mode:'copy'

    input:
    tuple val(sample_id), file(bam), file(bai) from bam_ch
    each chr from chr_list

    output:
    tuple val(sample_id), file("${sample_id}_${chr}.vcf") into vcf_ch

    script:
    """
    freebayes -f $params.ref -p 2 -r $chr $bam > ${sample_id}_${chr}.vcf
    """
}


process merge_vcf{

    publishDir "${params.outdir}/$sample_id", mode:'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), file(chrs) from vcf_ch.groupTuple(by:0, size: chr_list.size() , sort:true)

    output:
    tuple val(sample_id), file("${sample_id}.vcf") into merge_ch

    script:
    """
    bcftools concat -o ${sample_id}.vcf $chrs

    """
}

merge_ch.into { merge_ch2; stats_ch }

process Stats {
	
	publishDir "${params.outdir}/$sample_id/stats/", mode:'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), file(vcf) from stats_ch

    output:
    file("*")

    script:
    """
    vcfstats $vcf > ${sample_id}_vcfstats.stats

    Variants_summary.sh $vcf ${sample_id}_variants.txt ${sample_id}_trans.txt
    """

}

process Subset {
    
    publishDir "${params.outdir}/$sample_id/plots/", mode:'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), file(vcf) from merge_ch

    output:
    tuple val(sample_id), file("${sample_id}_subset.lqual"), file("${sample_id}_subset.ldepth.mean") into subset_ch

    script:
    """
    vcfrandomsample -r 0.01  $vcf > ${sample_id}_subset.vcf
    bgzip ${sample_id}_subset.vcf
    bcftools index ${sample_id}_subset.vcf.gz

    vcftools --gzvcf ${sample_id}_subset.vcf.gz --site-mean-depth --out ${sample_id}_subset
    vcftools --gzvcf ${sample_id}_subset.vcf.gz --site-quality --out ${sample_id}_subset
    """
}

process plot {
    
    publishDir "${params.outdir}/$sample_id/plots/", mode:'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), file(qual), file(depth) from subset_ch

    output:
    file("*.png") 

    script:
    """
    #!/usr/bin/env Rscript 
    
    library(ggplot2)
    library(dplyr)
    library(readr)


    var_qual <- read_delim( "${qual}", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

    var_depth <- read_delim("${depth}", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)


    a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "pink", colour = "red", alpha = 0.3) + theme_light() + labs(title = "${sample_id} quality subset 10% variants") + xlim(0,500)

    b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "pink", colour = "red", alpha = 0.3) + theme_light() + labs(title = "${sample_id} depth subset 10% variants") + xlim(0,500)

    plot(a)   
   
    ggsave(
      "${sample_id}_quality.png",
      plot = last_plot(),
      device = "png",
      scale=1,
      width = 15,
      height = 15,
      unit= "in")

    plot(b)   
   
    ggsave(
      "${sample_id}_depth.png",
      plot = last_plot(),
      device = "png",
      scale=1,
      width = 15,
      height = 15,
      unit= "in")
    """
}
