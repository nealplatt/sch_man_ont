import os
import pandas as pd
from snakemake.utils import min_version

#using
# - snakemake v6.9.1
# - using conda v22.9.0

# #ex:
# snakemake \
# --printshellcmds \
# --cluster 'qsub -V -cwd -j y -S /bin/bash -pe smp {threads} -q all.q -o {log} -e {log}' \
# --jobs 20 \
# --latency-wait 200 \
# --keep-going \
# --rerun-incomplete \
# --snake code/sch_man_ont.snk.py \
# --use-conda \
# --jobname snk.sv.{name}.jid{jobid}

##### set minimum snakemake version #####
min_version("6.9.1")

#set main project dir and work from there
proj_dir    = "/master/nplatt/sch_man_ont"
data_dir    = "{}/data".format(proj_dir)
results_dir = "{}/results".format(proj_dir)
envs_dir    = "{}/envs".format(proj_dir)
logs_dir    = "{}/logs".format(results_dir)
seq_dir     = "{}/raw_reads".format(results_dir)

#set ref genome
ref_genome = "{}/genome/SM_V10.fa".format(data_dir)

#insert list of populations
lab_pops =["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]
callers = ["cutesv", "sniffles", "debreak"]
#working on "nanovar", "mamnet" "svim",

localrules: 
    all,
    coverage_table,
    prep_survivor_r1,
    prep_survivor_r2,
    svjedi_to_tsv,
    tsvs_to_merged_table,
    extract_exons_from_gtf

rule all:
    input:      
        expand("{dir}/{lab_pop}.fastq",                         dir = seq_dir,     lab_pop = lab_pops),
        f"{seq_dir}/n50_raw_reads.csv",
        f"{results_dir}/nanoq/n50_filtered_reads.csv",
        f"{results_dir}/mosdepth/covs.cvs",
        expand("{dir}/nanoq/{lab_pop}.nanoq.fastq",             dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/minimap2/{lab_pop}.bam",                  dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/minimap2/{lab_pop}.flagstat.txt",         dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/mosdepth/{lab_pop}.mosdepth.summary.txt", dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/cutesv/{lab_pop}.cutesv.vcf",             dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/sniffles/{lab_pop}.sniffles.vcf",         dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/debreak/{lab_pop}.debreak.vcf",           dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/survivor_r1/{lab_pop}.survivor_r1.vcf",   dir = results_dir, lab_pop = lab_pops),
        f"{results_dir}/survivor_r2/survivor_r2.vcf",
        expand("{dir}/svjedi/{lab_pop}.jedi_genotype.vcf",   dir = results_dir, lab_pop = lab_pops),
        expand("{dir}/svjedi/{lab_pop}.jedi_genotype.tsv",   dir = results_dir, lab_pop = lab_pops),
        f"{results_dir}/svjedi/merged.jedi.tsv",
        expand("{dir}/svs_in_genes/{lab_pop}.{caller}.exons.bed", dir = results_dir, lab_pop = lab_pops, caller=callers)


rule n50_raw_reads:
    input:
       raw_bre_fq = seq_dir + "/smbre.fastq",
       raw_eg_fq  = seq_dir + "/smeg.fastq",
       raw_er_fq  = seq_dir + "/smle-pzq-er.fastq",
       raw_es_fq  = seq_dir + "/smle-pzq-es.fastq",
       raw_or_fq  = seq_dir + "/smor.fastq",
    output:
       n50_csv = seq_dir + "/n50_raw_reads.csv"
    threads:
        1
    log:
        logs_dir + "/n50_raw_reads.log"
    conda:
        envs_dir + "/n50.yml"
    shell:
        """
        n50 \
            --format csv \
            {input.raw_bre_fq} \
            {input.raw_eg_fq} \
            {input.raw_er_fq} \
            {input.raw_es_fq} \
            {input.raw_or_fq} \
            >{output.n50_csv}
        """

rule nanoq:
    input:
        raw_fq = seq_dir + "/{lab_pop}.fastq"
    output:
        nanoq_fq  = results_dir + "/nanoq/{lab_pop}.nanoq.fastq",
        nanoq_txt = results_dir + "/nanoq/{lab_pop}.nanoq.txt"
    threads:
        1
    params:
        min_len=500,
        min_qual=7,
        trim_start=25,
        trim_end=25
    log:
        logs_dir + "/nanoq.{lab_pop}.log"
    conda:
        envs_dir + "/nanoq.yml"
    shell:
        """
        nanoq \
            --input {input.raw_fq} \
            --output {output.nanoq_fq} \
            --min-len {params.min_len} \
            --min-qual {params.min_qual} \
            --trim-start {params.trim_start} \
            --trim-end {params.trim_start} \
            --report {output.nanoq_txt}
        """

rule n50_filtered_reads:
    input:
       filt_bre_fq = results_dir + "/nanoq/smbre.nanoq.fastq",
       filt_eg_fq  = results_dir + "/nanoq/smeg.nanoq.fastq",
       filt_er_fq  = results_dir + "/nanoq/smle-pzq-er.nanoq.fastq",
       filt_es_fq  = results_dir + "/nanoq/smle-pzq-es.nanoq.fastq",
       filt_or_fq  = results_dir + "/nanoq/smor.nanoq.fastq"
    output:
       n50_csv = results_dir + "/nanoq/n50_filtered_reads.csv"
    threads:
        1
    log:
        logs_dir + "/n50_filtered_reads.log"
    conda:
        envs_dir + "/n50.yml"
    shell:
        """
        n50 \
            --format csv \
            {input.filt_bre_fq} \
            {input.filt_eg_fq} \
            {input.filt_er_fq} \
            {input.filt_es_fq} \
            {input.filt_or_fq} \
            >{output.n50_csv}
        """

#map reads
rule minimap2:
    input:
        fq = results_dir + "/nanoq/{lab_pop}.nanoq.fastq",
        ref_fas = ref_genome
    output:
        sam  = results_dir + "/minimap2/{lab_pop}.sam",
    threads:
        48
    log:
        logs_dir + "/minimap2.{lab_pop}.log"
    conda:
        envs_dir + "/minimap2.yml"
    shell:
        """
        minimap2 \
            -a \
            -o {output.sam} \
            -t {threads} \
            -x map-ont \
            {input.ref_fas} \
            {input.fq}
        """
#map reads
rule flagstat:
    input:
        sam = results_dir + "/minimap2/{lab_pop}.sam"
    output:
        txt  = results_dir + "/minimap2/{lab_pop}.flagstat.txt"
    threads:
        1
    log:
        logs_dir + "/flagstat.{lab_pop}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools flagstat {input.sam} >{output.txt}
        """

rule sort_bam:
    input:
        sam = results_dir + "/minimap2/{lab_pop}.sam"
    output:
        bam  = results_dir + "/minimap2/{lab_pop}.bam"
    threads:
        12
    params:
        tmp_prefix = "tmp.{lab_pop}"
    log:
        logs_dir + "/sort_bam.{lab_pop}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools sort -@ {threads} -T {params.tmp_prefix} -O bam -o {output.bam} {input.sam}
        """

rule index_bam:
    input:
        bam = results_dir + "/minimap2/{lab_pop}.bam"
    output:
        bai  = results_dir + "/minimap2/{lab_pop}.bam.bai"
    threads:
        1
    log:
        logs_dir + "/index_bam.{lab_pop}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools index {input.bam}
        """

rule mosdepth:
    input:
        bam = results_dir + "/minimap2/{lab_pop}.bam",
        bai = results_dir + "/minimap2/{lab_pop}.bam.bai"
    output:
        dis = results_dir + "/mosdepth/{lab_pop}.mosdepth.global.dist.txt",
        sum = results_dir + "/mosdepth/{lab_pop}.mosdepth.summary.txt",
        bed = results_dir + "/mosdepth/{lab_pop}.per-base.bed.gz",
        csi = results_dir + "/mosdepth/{lab_pop}.per-base.bed.gz.csi"
    threads:
        4
    params:
        out_prefix=results_dir + "/mosdepth/{lab_pop}"
    log:
        logs_dir + "/mosdepth.{lab_pop}.log"
    conda:
        envs_dir + "/mosdepth.yml"
    shell:
        """
        mosdepth --threads {threads} {params.out_prefix} {input.bam}
        """

rule coverage_table:
    input:
        bre_sum = results_dir + "/mosdepth/smbre.mosdepth.summary.txt",
        eg_sum  = results_dir + "/mosdepth/smeg.mosdepth.summary.txt",
        er_sum  = results_dir + "/mosdepth/smle-pzq-er.mosdepth.summary.txt",
        es_sum  = results_dir + "/mosdepth/smle-pzq-es.mosdepth.summary.txt",
        or_sum  = results_dir + "/mosdepth/smor.mosdepth.summary.txt"
    output:
        cov_cvs = results_dir + "/mosdepth/covs.cvs",
    params:
        mosdepth_term = results_dir + "/mosdepth/*mosdepth.summary.txt"
    threads:
        1
    log:
        logs_dir + "/coverage_table.log"
    shell:
        """
        echo "sample,length,bases,mean,min,max">{output.cov_cvs}
        (grep "total" {params.mosdepth_term} | sed 's/:/\t/' | sed 's/\t/,/g' | cut -f1,3-7 -d"," | sed 's/.mosdepth.summary.txt//') >>{output.cov_cvs}
        """

rule cutesv:
    input:
        bam = results_dir + "/minimap2/{lab_pop}.bam",
        bai = results_dir + "/minimap2/{lab_pop}.bam.bai",
        ref_fas = ref_genome
    output:
        cutesv_vcf  = results_dir + "/cutesv/{lab_pop}.cutesv.vcf"
    threads:
        48
    params:
        min_mapq=20,
        min_read_len=500,
        md=100,
        mi=100,
        min_support=5,
        min_size=50,
        max_size=200000,
        outdir = results_dir + "/cutesv/{lab_pop}"
    log:
        logs_dir + "/cutesv.{lab_pop}.log"
    conda:
        envs_dir + "/cutesv.yml"
    shell:
        """

        if [ -d {params.outdir} ]; then
            rm -rf {params.outdir}
        fi
        mkdir -p {params.outdir}
        
        cuteSV \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads {threads} \
            --min_mapq {params.min_mapq} \
            --min_read_len {params.min_read_len} \
            -md {params.md} \
            -mi {params.mi} \
            --min_support {params.min_support} \
            --min_size {params.min_size} \
            --max_size {params.max_size} \
            --genotype \
            {input.bam} \
            {input.ref_fas} \
            {output.cutesv_vcf} \
            {params.outdir}        
    """

rule sniffles:
    input:
        bam = results_dir + "/minimap2/{lab_pop}.bam",
        bai = results_dir + "/minimap2/{lab_pop}.bam.bai",
        ref_fas = ref_genome
    output:
        sniffles_vcf  = results_dir + "/sniffles/{lab_pop}.sniffles.vcf",
        sniffles_snf  = results_dir + "/sniffles/{lab_pop}.sniffles.snf"
    threads:
        48
    params:
        sample_id="{wildcards.lab_pop}",
        mosaic_max=0.95,
        mosaic_min=0.5,
        minsvlen=50
    log:
        logs_dir + "/sniffles.{lab_pop}.log"
    conda:
        envs_dir + "/sniffles.yml"
    shell:
        """
        sniffles \
            --input {input.bam} \
            --vcf {output.sniffles_vcf} \
            --snf {output.sniffles_snf} \
            --threads {threads} \
            --reference {input.ref_fas} \
            --allow-overwrite \
            --sample-id {params.sample_id} \
            --minsvlen {params.minsvlen} \
            --mosaic \
            --mosaic-af-max {params.mosaic_max} \
            --mosaic-af-min {params.mosaic_min} \
            --qc-output-all \
            --output-rnames        
        """

def get_mean_cov(lab_pop):
    
    file_path=f"{results_dir}/mosdepth/{lab_pop}.mosdepth.summary.txt"
    
    df = pd.read_csv(file_path, sep="\t", header=0)
    mean_cov = df.loc[df["chrom"] == "total", "mean"].values[0].astype(float)
    
    return mean_cov

rule debreak:
    input:
        bam = results_dir + "/minimap2/{lab_pop}.bam",
        bai = results_dir + "/minimap2/{lab_pop}.bam.bai",
        ref_fas = ref_genome,
        cov = results_dir + "/mosdepth/{lab_pop}.mosdepth.summary.txt"
    output:
        debreak_vcf  = results_dir + "/debreak/{lab_pop}/debreak.vcf",
        renamed_vcf  = results_dir + "/debreak/{lab_pop}.debreak.vcf"
    threads:
        48
    params:
        outpath= results_dir + "/debreak/{lab_pop}",
        min_size=50,
        max_size=200000,
        depth=lambda wildcards: get_mean_cov(wildcards.lab_pop),
        aligner="minimap2",
        max_cov=lambda wildcards: int(get_mean_cov(wildcards.lab_pop) * 2),
        min_quality=20,
        min_support=lambda wildcards: int(get_mean_cov(wildcards.lab_pop) / 4)
    log:
        logs_dir + "/debreak.{lab_pop}.log"
    conda:
        envs_dir + "/debreak.yml"
    shell:
        """
        debreak \
            --bam {input.bam} \
            --poa \
            --outpath {params.outpath} \
            --min_size {params.min_size} \
            --max_size {params.max_size} \
            --depth {params.depth} \
            --aligner {params.aligner} \
            --thread {threads} \
            --rescue_dup \
            --rescue_large_ins \
            --ref {input.ref_fas} \
            --maxcov {params.max_cov} \
            --min_quality {params.min_quality} \
            --min_support {params.min_support}

        cp {output.debreak_vcf} {output.renamed_vcf}
        """

# rule filter_gt_calls:
#     input:
#         vcf   = results_dir + "/{caller}/{lab_pop}.{caller}.vcf",
#     output:
#         filtered_vcf   = results_dir + "/{caller}/{lab_pop}.{caller}.filtered.vcf",
#     threads:
#         1
#     log:
#         logs_dir + "/filter_gt_calls.{lab_pop}.{caller}.log"
#     conda:
#         envs_dir + "/vcftools.yml"
#     shell:
#         """
#         vcftools \
#             --vcf {input.vcf} \
#             --remove-filtered-all \
#             --recode \
#             --recode-INFO-all \
#             --stdout \
#             >{output.filtered_vcf}
#         """

rule prep_survivor_r1:
    input:
        cutesv_vcf   = results_dir + "/cutesv/{lab_pop}.cutesv.vcf",
        sniffles_vcf = results_dir + "/sniffles/{lab_pop}.sniffles.vcf",
        debreak_vcf  = results_dir + "/debreak/{lab_pop}.debreak.vcf",
    output:
        vcf_list  = results_dir + "/survivor_r1/{lab_pop}.list"
    threads:
        1
    log:
        logs_dir + "/prep_survivor_r1.{lab_pop}.log"
    shell:
        """
        echo {input.cutesv_vcf}    >{output.vcf_list}
        echo {input.sniffles_vcf} >>{output.vcf_list}
        echo {input.debreak_vcf} >>{output.vcf_list}
        """
rule survivor_r1:
    input:
        cutesv_vcf   = results_dir + "/cutesv/{lab_pop}.cutesv.filtered.vcf",
        sniffles_vcf = results_dir + "/sniffles/{lab_pop}.sniffles.filtered.vcf",
        debreak_vcf  = results_dir + "/debreak/{lab_pop}.debreak.filtered.vcf",
        vcf_list  = results_dir + "/survivor_r1/{lab_pop}.list"
    output:
        merged_vcf  = results_dir + "/survivor_r1/{lab_pop}.survivor_r1.vcf"
    threads:
        12
    params:
        max_distance_between_breakpoints=0.05,
        min_calls=3,
        sv_type=1,
        sv_strand=0,
        dist_estimate=0,
        min_size=50
    log:
        logs_dir + "/survivor_r1.{lab_pop}.log"
    conda:
        envs_dir + "/survivor.yml"
    shell:
        """
        SURVIVOR \
            merge \
            {input.vcf_list} \
            {params.max_distance_between_breakpoints} \
            {params.min_calls} \
            {params.sv_type} \
            {params.sv_strand} \
            {params.dist_estimate} \
            {params.min_size} \
            {output.merged_vcf}
        """

rule prep_survivor_r2:
    input:
        smbre_vcf        = results_dir + "/survivor_r1/smbre.survivor_r1.vcf",
        smor_vcf         = results_dir + "/survivor_r1/smor.survivor_r1.vcf",
        smle_pzq_er_vcf  = results_dir + "/survivor_r1/smle-pzq-er.survivor_r1.vcf",
        smle_pzq_es_vcf  = results_dir + "/survivor_r1/smle-pzq-es.survivor_r1.vcf",
        smeg_vcf         = results_dir + "/survivor_r1/smeg.survivor_r1.vcf"
    output:
        vcf_list  = results_dir + "/survivor_r2/pop_vcfs.list"
    threads:
        1
    log:
        logs_dir + "/prep_survivor_r2.log"
    shell:
        """
        echo {input.smbre_vcf}        >{output.vcf_list}
        echo {input.smor_vcf}        >>{output.vcf_list}
        echo {input.smle_pzq_er_vcf} >>{output.vcf_list}
        echo {input.smle_pzq_es_vcf} >>{output.vcf_list}
        echo {input.smeg_vcf}        >>{output.vcf_list}
        """

rule survivor_r2:
    input:
        smbre_vcf        = results_dir + "/survivor_r1/smbre.survivor_r1.vcf",
        smor_vcf         = results_dir + "/survivor_r1/smor.survivor_r1.vcf",
        smle_pzq_er_vcf  = results_dir + "/survivor_r1/smle-pzq-er.survivor_r1.vcf",
        smle_pzq_es_vcf  = results_dir + "/survivor_r1/smle-pzq-es.survivor_r1.vcf",
        smeg_vcf         = results_dir + "/survivor_r1/smeg.survivor_r1.vcf",
        vcf_list         = results_dir + "/survivor_r2/pop_vcfs.list"

    output:
        merged_vcf  = results_dir + "/survivor_r2/survivor_r2.vcf"
    threads:
        12
    params:
        max_distance_between_breakpoints=0.05,
        min_calls=1,
        sv_type=1,
        sv_strand=0,
        dist_estimate=0,
        min_size=50,
    log:
        logs_dir + "/survivor_r2.log"
    conda:
        envs_dir + "/survivor.yml"
    shell:
        """
        SURVIVOR \
            merge \
            {input.vcf_list} \
            {params.max_distance_between_breakpoints} \
            {params.min_calls} \
            {params.sv_type} \
            {params.sv_strand} \
            {params.dist_estimate} \
            {params.min_size} \
            {output.merged_vcf}
        """

#svjedi
rule svjedi:
    input:
        merged_vcf  = results_dir + "/survivor_r2/survivor_r2.vcf",
        nanoq_fq  = results_dir + "/nanoq/{lab_pop}.nanoq.fastq",
        ref_fas = ref_genome
    output:
        jedi_vcf  = results_dir + "/svjedi/{lab_pop}.jedi_genotype.vcf",
    threads:
        48
    params:
        prefix=results_dir + "/svjedi/{lab_pop}.jedi",
        min_support=5,
    log:
        logs_dir + "/svjedi.{lab_pop}.log"
    conda:
        envs_dir + "/svjedi.yml"
    shell:
        """
        svjedi-graph.py \
            --vcf {input.merged_vcf} \
            --ref {input.ref_fas} \
            --reads {input.nanoq_fq} \
            --prefix {params.prefix} \
            --threads {threads} \
            --minsupport {params.min_support} 
        """

    rule svjedi_to_tsv:
        input:
            jedi_vcf  = results_dir + "/svjedi/{lab_pop}.jedi_genotype.vcf",
        output:
            jedi_tsv  = results_dir + "/svjedi/{lab_pop}.jedi_genotype.tsv",
        params:
            sample_id = "{lab_pop}"
        threads:
            1
        log:
            logs_dir + "/svjedi_to_tsv.{lab_pop}.log"
        conda:
            envs_dir + "/bcftools.yml"
        shell:
            """
            echo -e "chr\tpos\tid\tref\talt\tqual\tfilter\tsvlen\tsvtype\tchr2\tend\tcipos\tgt.{params.sample_id}\tdp.{params.sample_id}\tad.{params.sample_id}\tpl.{params.sample_id}" \
            > {output.jedi_tsv}
            bcftools query --format "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%SVLEN\t%SVTYPE\t%CHR2\t%END\t%CIPOS\t[%GT\t%DP\t%AD\t%PL]\n" \
            {input.jedi_vcf} \
            >> {output.jedi_tsv}
            """
    
    rule tsvs_to_merged_table:
        input:
            bre_jedi_tsv = results_dir + "/svjedi/smbre.jedi_genotype.tsv",
            eg_jedi_tsv  = results_dir + "/svjedi/smeg.jedi_genotype.tsv",
            er_jedi_tsv  = results_dir + "/svjedi/smle-pzq-er.jedi_genotype.tsv",
            es_jedi_tsv  = results_dir + "/svjedi/smle-pzq-es.jedi_genotype.tsv",
            or_jedi_tsv  = results_dir + "/svjedi/smor.jedi_genotype.tsv"
        output:
            merged_tsv  = results_dir + "/svjedi/merged.jedi.tsv",
        threads:
            1
        log:
            logs_dir + "/tsvs_to_merged_table.log"
        shell:
            """
            cat {input.bre_jedi_tsv} >bre.tmp
            cut -f13-16 {input.eg_jedi_tsv} >eg.tmp
            cut -f13-16 {input.er_jedi_tsv} >er.tmp
            cut -f13-16 {input.es_jedi_tsv} >es.tmp
            cut -f13-16 {input.or_jedi_tsv} >or.tmp
    
            paste bre.tmp eg.tmp er.tmp es.tmp or.tmp >{output.merged_tsv}
    
            rm bre.tmp eg.tmp er.tmp es.tmp or.tmp
            """
###########################################################################################################################


rule debreak_to_bed:
    input:
        vcf  = results_dir + "/{caller}/{lab_pop}.{caller}.vcf",
    output:
        bed  = results_dir + "/svs_in_genes/{lab_pop}.{caller}.bed",
    threads:
        1
    log:
        logs_dir + "/debreak_to_bed.{lab_pop}.{caller}log"
    conda:
        envs_dir + "/bcftools.yml"
    shell:
        """
        bcftools query --format "%CHROM\t%POS\t%END\t%ID\n" {input.vcf}> {output.bed}
        """

rule extract_exons_from_gtf:
    input:
        gtf  = data_dir + "/genome/SM_V10.gtf"
    output:
        exon_gtf  = results_dir + "/svs_in_genes/exons.gtf"
    threads:
        1
    log:
        logs_dir + "/extract_exons_from_gtf.log"
    conda:
        envs_dir + "/bcftools.yml"
    shell:
        """
        grep -i exon {input.gtf} >{output.exon_gtf}
        """

rule intersect_svs_vs_genes:
    input:
        bed  = results_dir + "/svs_in_genes/{lab_pop}.{caller}.bed",
        exon_gtf  = results_dir + "/svs_in_genes/exons.gtf"
    output:
        bed  = results_dir + "/svs_in_genes/{lab_pop}.{caller}.exons.bed",
    threads:
        1
    log:
        logs_dir + "/intersect_svs_vs_genes.{lab_pop}.{caller}.log"
    conda:
        envs_dir + "/bedtools.yml"
    shell:
        """
         bedtools intersect -a {input.bed} -b {input.exon_gtf} -u >{output.bed}
        """







######################################## UNUSED SV CALLERS ##############################################################
# rule prep_nanosv_faidx:
#     input:
#         ref_fas = ref_genome
#     output:
#         ref_fai  = data_dir  + "/genome/SM_V10.fa.fai"
#     threads:
#         1
#     log:
#         logs_dir + "/prep_nanosv_faidx.log"
#     conda:
#         envs_dir + "/samtools.yml"
#     shell:
#         """
#         samtools faidx {input.ref_fas}
#         """

# rule prep_nanosv_bed:
#     input:
#         ref_fai = data_dir  + "/genome/SM_V10.fa.fai"
#     output:
#         bed         = results_dir + "/nanosv/random_sites.bed",
#         sorted_bed  = results_dir + "/nanosv/random_sites.sorted.bed"
#     threads:
#         1
#     params:
#         window_size = 1,
#         n_windows   = 1000000,
#         random_seed = 123456789
#     log:
#         logs_dir + "/prep_nanosv_bed.log"
#     conda:
#         envs_dir + "/bedtools.yml"
#     shell:
#         """
#         bedtools random -l {params.window_size} -n {params.n_windows} -seed {params.random_seed} -g {input.ref_fai} >{output.bed}
#         bedtools sort -i {output.bed} >{output.sorted_bed}
#         """

# rule nanosv:
#     input:
#         bam = results_dir + "/minimap2/{lab_pop}.bam",
#         bed  = results_dir + "/nanosv/random_sites.sorted.bed"
#     output:
#         nanosv_vcf   = results_dir + "/nanosv/{lab_pop}.nanosv.vcf"
#     threads:
#         48
#     params:
#         sambamba_path="/master/nplatt/anaconda3/envs/nanosv/bin/sambamba"
#     log:
#         logs_dir + "/nanosv.{lab_pop}.log"
#     conda:
#         envs_dir + "/nanosv.yml"
#     shell:
#         """
#         NanoSv \
#             --threads {threads} \
#             --sambamba {params.sambamba_path}\
#             --bed {input.bed} \
#             --output {output.nanosv_vcf} \
#             {input.bam} 
#         """




#         """
# #nanovar
# rule nanovar:
#     input:
#         bam = rule.sort_bam.output.bam,
#         bai = rule.sort_bam.output.bai,
#         ref_fas = ref_genome,
#     output:
#         nanovar_vcf   = f"{results_dir}/nanovar/{pop}.nanovar.pass.vcf",
#         nanovar_html  = f"{results_dir}/nanovar/{pop}.nanovar.pass.report.html"
#     threads:
#         48
#     params:
#         outpath=f"{results_dir}/nanovar"
#     log:
#         logs_dir + "/nanovar.{pop}.log"
#     conda:
#         envs_dir + "/nanovar.yml"
#     shell:
#         """
#         nanovar [Options] \
#             -t {threads} \
#             {input.bam} \
#             {input.ref} \
#             {params.outpath}     
#         """



# #svision
# rule svision:
#     input:
#         bam = rule.sort_bam.output.bam,
#         bai = rule.sort_bam.output.bai,
#         ref_fas = ref_genome,
#     output:
#         nanovar_vcf   = f"{results_dir}/svision/{pop}.svision.vcf
#     threads:
#         48
#     params:
#         outpath=f"{results_dir}/nanovar"
#     log:
#         logs_dir + "/nanovar.{pop}.log"
#     conda:
#         envs_dir + "/nanovar.yml"
#     shell:
#         """
#         nanovar [Options] \
#             -t {threads} \
#             {input.bam} \
#             {input.ref} \
#             {params.outpath}     
#         """

# rule sensv:
#     input:
#         fastq   = results_dir + "/nanoq/{lab_pop}.nanoq.fastq",
#         ref_fas = ref_genome,
#     output:
#         sensv_vcf   = results_dir + "/sensv/{lab_pop}/{lab_pop}_final.vcf",
#         vcf   = results_dir + "/sensv/{lab_pop}.sensv.vcf"
#     threads:
#         48
#     params:
#         sample_name = "{lab_pop}",
#         out_prefix  = results_dir + "/sensv/{lab_pop}/{lab_pop}",
#         min_sv_size =50,
#         max_sv_size =200000
#     log:
#         logs_dir + "/sensv.{lab_pop}.log"
#     conda:
#         envs_dir + "/sensv.yml"
#     shell:
#         """
#         sensv \
#             --sample_name {params.sample_name} \
#             --fastq {input.fastq}\
#             --ref {input.ref_fas} \
#             --output_prefix {params.out_prefix} \
#             --min_sv_size {params.min_sv_size} \
#             --max_sv_size {params.max_sv_size} \
#             --nprocs {threads}

#         cp {output.sens_vcf} {output.vcf} 
#         """

# rule duet:
#     input:
#         fastq   = results_dir + "/nanoq/{lab_pop}.nanoq.fastq",
#         ref_fas = ref_genome,
#     output:
#         duet_vcf = results_dir + "/sensv/{lab_pop}/{lab_pop}_final.vcf",
#         vcf      = results_dir + "/sensv/{lab_pop}.sensv.vcf"
#     threads:
#         48
#     params:
#         out_dir  = results_dir + "/duet/{lab_pop}",
#         include_all_ctgs = "{lab_pop}",
#         sv_min_size = 50,
#         min_allele_frequency =0.25,
#         min_support_read = 2,
#         sv_caller = "svim"
#     log:
#         logs_dir + "/duet.{lab_pop}.log"
#     conda:
#         envs_dir + "/duet.yml"
#     shell:
#         """
#         duet \
#             {input.bam} \
#             {input.ref_fa} \
#             {params.out_dir}
#             --include_all_ctgs \
#             --sv_min_size {params.sv_min_size} \
#             --min_allele_frequency {params.min_allele_frequency} \
#             --min_support_read {params.min_support_read} \
#             --sv_caller {params.sv_caller} \
#             --thread {threads}

#         cp {output.duet_vcf} {output.vcf} 
#         """