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
            Mitomania
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            'USAGE'
            
            nextflow run variantcaller.nf --bams /home/filip/Documents/nextflowpipes/varientcalls/data/'*_mem_merged_sorted.{bam,bam.bai}' --outdir /home/filip/Documents/nextflowpipes/varientcalls/results


            'Mandatory arguments:'
            --outdir                 Path to output directory
            
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
         Mitomania 
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

    tag "$sample_id"

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
    file("${sample_id}.stats")    

    script:
    """
    bcftools concat -o ${sample_id}.vcf $chrs
    vcfstats ${sample_id}.vcf > ${sample_id}.stats
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
    bgzip _${sample_id}_subset.vcf
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
    file("${sample_id}_subset.*") 

    script:
    """
    #!/usr/bin/env Rscript 

    var_qual <- read_delim( "${qual}", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

    var_depth <- read_delim("${depth}", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)


    a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "pink", colour = "red", alpha = 0.3) + theme_light() + labs(title = "Indivxxx quality subset 10% variants")

    b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "pink", colour = "red", alpha = 0.3) + theme_light() + labs(title = "Indivxxx depth subset 10% variants")

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
