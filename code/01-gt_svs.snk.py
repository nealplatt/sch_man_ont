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
lab_pops = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]


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
        expand("{dir}/{lab_pop}.fastq",                         dir = seq_dir,     lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        f"{seq_dir}/n50_raw_reads.csv",
        f"{results_dir}/nanoq/n50_filtered_reads.csv",
        f"{results_dir}/mosdepth/covs.cvs",
        expand("{dir}/nanoq/{lab_pop}.nanoq.fastq",             dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        expand("{dir}/minimap2/{lab_pop}.bam",                  dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        expand("{dir}/minimap2/{lab_pop}.flagstats.txt",        dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        expand("{dir}/minimap2/{lab_pop}.stats.txt",            dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        expand("{dir}/mosdepth/{lab_pop}.mosdepth.summary.txt", dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        expand("{dir}/cutesv/{lab_pop}.cutesv.vcf",             dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        expand("{dir}/cutesv/{lab_pop}.cutesv.filtered.vcf",    dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        f"{results_dir}/survivor/survivor.vcf",
        expand("{dir}/svjedi/{lab_pop}.jedi_genotype.vcf",      dir = results_dir, lab_pop = ["smbre", "smeg", "smle-pzq-er", "smle-pzq-es", "smor"]),
        f"{results_dir}/svjedi/merged.jedi.tsv",


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
        nanoq_fq  = results_dir + "/nanoq/{lab_pop, sm.*\.}.nanoq.fastq",
        nanoq_txt = results_dir + "/nanoq/{lab_pop, sm.*\.}.nanoq.txt"
    threads:
        1
    params:
        min_len=500,
        min_qual=7,
        trim_start=25,
        trim_end=25
    log:
        logs_dir + "/nanoq.{lab_pop, sm.*\.}.log"
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
        sam  = results_dir + "/minimap2/{lab_pop, sm.*\.}.sam",
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
rule flagstats:
    input:
        sam = results_dir + "/minimap2/{lab_pop}.sam"
    output:
        txt  = results_dir + "/minimap2/{lab_pop, sm.*\.}.flagstats.txt"
    threads:
        1
    log:
        logs_dir + "/flagstats.{lab_pop}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools flagstat {input.sam} >{output.txt}
        """

rule stats:
    input:
        sam = results_dir + "/minimap2/{lab_pop}.sam"
    output:
        txt  = results_dir + "/minimap2/{lab_pop, sm.*\.}.stats.txt"
    threads:
        1
    log:
        logs_dir + "/stats.{lab_pop}.log"
    conda:
        envs_dir + "/samtools.yml"
    shell:
        """
        samtools stats {input.sam} >{output.txt}
        """

rule sort_bam:
    input:
        sam = results_dir + "/minimap2/{lab_pop}.sam"
    output:
        bam  = results_dir + "/minimap2/{lab_pop, sm.*\.}.bam"
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
        bai  = results_dir + "/minimap2/{lab_pop, sm.*\.}.bam.bai"
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
        dis = results_dir + "/mosdepth/{lab_pop, sm.*\.}.mosdepth.global.dist.txt",
        sum = results_dir + "/mosdepth/{lab_pop, sm.*\.}.mosdepth.summary.txt",
        bed = results_dir + "/mosdepth/{lab_pop, sm.*\.}.per-base.bed.gz",
        csi = results_dir + "/mosdepth/{lab_pop, sm.*\.}.per-base.bed.gz.csi"
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

        rm -r {params.outdir}
    """

rule filter_cutesv:
    input:
        cutesv_vcf = results_dir + "/cutesv/{lab_pop}.cutesv.vcf"
    output:
        filtered_vcf = results_dir + "/cutesv/{lab_pop}.cutesv.filtered.vcf",
    threads:
        4
    params:
        min_q  = 10,
        min_dp = 10,
    log:
        logs_dir + "/filter_cutesv.{lab_pop}.log"
    conda:
        envs_dir + "/vcftools.yml"
    shell:
        """
        vcftools --vcf {input.cutesv_vcf} --minQ {params.min_q} --minDP {params.min_dp} --remove-filtered-all --max-missing-count 0 --recode --recode-INFO-all --stdout >{output.filtered_vcf}
        """

rule survivor:
    input:
        vcfs = expand(results_dir + "/cutesv/{lab_pop}.cutesv.filtered.vcf", lab_pop = lab_pops)
    output:
        merged_vcf  = results_dir + "/survivor/survivor.vcf",
        vcf_list  = results_dir + "/survivor/pop_vcfs.list",
        survivor_dir = directory(results_dir + "/survivor")
    threads:
        12
    params:
        max_distance_between_breakpoints=100,
        min_calls=1,
        sv_type=1,
        sv_strand=0,
        dist_estimate=0,
        min_size=50,
    log:
        logs_dir + "/survivor.log"
    conda:
        envs_dir + "/survivor.yml"
    shell:
        """

        if [ ! -d {output.survivor_dir} ]; then
            mkdir -p {output.survivor_dir}
        fi
        
        echo {input.vcfs} | sed 's/ /\\n/g' >{output.vcf_list}
        
        SURVIVOR \
            merge \
            {output.vcf_list} \
            {params.max_distance_between_breakpoints} \
            {params.min_calls} \
            {params.sv_type} \
            {params.sv_strand} \
            {params.dist_estimate} \
            {params.min_size} \
            {output.merged_vcf}
        """

rule svjedi:
    input:
        merged_vcf  = results_dir + "/survivor/survivor.vcf",
        nanoq_fq  = results_dir + "/nanoq/{lab_pop}.nanoq.fastq",
        ref_fas = ref_genome
    output:
        jedi_vcf  = results_dir + "/svjedi/{lab_pop}.jedi_genotype.vcf",
    threads:
        48
    params:
        prefix=results_dir + "/svjedi/{lab_pop}.jedi",
        min_support=2,
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

